#%%
from rifl import analyze
import rifl.utils as utils

from sortedcontainers import SortedSet
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import scipy.stats as stats

import pathlib

MIN_INTENSITY_PER_CONTROL_CHANNEL = 5000
STANDARD_DEVIATION_FILTER_THRESHOLD = 0.5
MIN_NUM_UNIQUE_SEQUENCES = 2
MIN_NUM_DATASETS = 2
MIN_NUM_DATASETS_HYDROXYLAMINE = 2
DYNAMIC_FOLD_CHANGE_THRESHOLD = 3
CHASE_FOLD_CHANGE_THRESHOLD = 2
HYDROXYLAMINE_REDUCTION_THRESHOLD = 75
DATA_OUTPUT_PATH = pathlib.Path('output')
DATA_INPUT_PATH = pathlib.Path('input')
FIG_OUTPUT_PATH = pathlib.Path('plots')
UNIPROT_SEARCH_DB_PATH = DATA_INPUT_PATH.joinpath('20181024_human_reviewed_uniprot_GFP_NRAS.fasta')
CACHE_DIR = pathlib.Path('.cache')


#%%

def analysis_to_df(analysis, channel_layout: dict, control_channels: str):
    m = analysis.data_model

    query = (m
        .select(
            m.id,
            m.experiment,
            m.filename,
            m.scan_num,
            m.uniprot,
            m.symbol,
            m.description,
            m.unique_sequence,
            m.clean_sequence,
            m.channel_intensities,
        )
        .where(
            (m.experiment_id.in_(analysis.experiment_ids_included))
        )
    )

    df = pd.DataFrame.from_records(list(query.dicts()))
    df = df.set_index('id')

    for channel in channel_layout:
        channel_num, channel_name = next(iter(channel.items()))
        new_col_name = f'{channel_name}.intensity'
        # append number to column to distinguish duplicates
        new_col_name = f'{new_col_name}_{len(df.filter(regex=new_col_name).columns)}'
        df[new_col_name] = df.channel_intensities.str.get(channel_num - 1)

    df = df.drop(columns=['channel_intensities'])

    control_cols = df.filter(regex=f'{control_channels}.intensity')
    data_cols = df.filter(regex=f'.+\.intensity')

    percentages = data_cols.div(control_cols.mean(axis=1), axis=0).mul(100, axis=0)
    percentages = percentages.rename(columns={
        c: c.replace('intensity', 'percent_of_control')
        for c in df.filter(regex='\.intensity').columns
    })

    df = df.join(percentages)

    return df

def output_psm_level_report(df, output_stem: str):
    df = df[~df.uniprot.str.startswith('Reverse_')]
    df = df[~df.uniprot.str.startswith('contaminant')]
    psm_df = df.set_index(['uniprot', 'symbol', 'description', 'clean_sequence']).sort_index(level=1)
    psm_df = psm_df.drop(columns=['unique_sequence'])

    psm_df.to_excel(
        utils.get_timestamped_report_path(f'{output_stem}_{{}}.xlsx', DATA_OUTPUT_PATH),
    )
    return psm_df

def group_by_protein_and_filter(df, agg_fun, min_unique_sequences: int, min_num_datasets: int):
    result = df.groupby(['uniprot', 'experiment']).apply(agg_fun)
    # result = result[(result['num_unique_peptides'] >= min_unique_sequences) | (result.index.get_level_values('uniprot') == 'P01111')]
    result = result[result['num_unique_peptides'] >= min_unique_sequences]
    # result = result.groupby('uniprot').apply(lambda x: x if len(x) >= min_num_datasets else None)
    # result = result.unstack()

    return (result
        .groupby(level=0, group_keys=False)
        .apply(lambda x: x if len(x) >= min_num_datasets else None)
        .unstack()
    )

def add_meta(df, meta):
    meta = meta.drop_duplicates()
    meta = meta[meta.uniprot.isin(df.index)]
    meta = meta.set_index('uniprot')
    df = pd.merge(df, meta, left_index=True, right_index=True)
    return df

def filter_timecourse_data(
    output_name,
    params_filename,
    uniprots_passing_ha_filter,
    min_num_unique_sequences=1,
    min_num_datasets=1,
    output_psm_level=False,
    output_to_excel=True,
):
    analysis, params = analyze(params_filename, user='radu')
    dataset = analysis.datasets[0]
    df = analysis_to_df(analysis, dataset.channel_layout, dataset.control_channels)

    if output_psm_level:
        output_psm_level_report(df, f'unfiltered_timecourse_{output_name}_psm_level')

    df = df[~df.index.isin(analysis.filtered_out)]
    df = df[~df.uniprot.str.startswith('contaminant_')]
    df = df[df.uniprot.isin(uniprots_passing_ha_filter)]

    def agg(x):
        percentages = x.filter(regex='percent_of_control')
        cv = percentages.apply(stats.variation)
        have_cv_ge = percentages[cv[cv.ge(0.5)].index.values]

        if len(x) > 2 and not have_cv_ge.empty:
            to_filter_out = have_cv_ge.apply(stats.zscore).abs().idxmax().values
            x = x[~x.index.isin(to_filter_out)]

        result = x.filter(regex='percent_of_control').mean()
        result['num_unique_peptides'] = x['clean_sequence'].nunique()
        return result

    # group by protein
    result = group_by_protein_and_filter(df, agg, min_num_unique_sequences, min_num_datasets)

    # result.columns = pd.MultiIndex.from_tuples(
    #     (i.split('.percent')[0], j) for i, j in result.columns
    # )

    result.columns.set_levels([f'Replicate {i}' for i in range(1, len(dataset.experiments) + 1)], level=1, inplace=True)
    means = result.drop(columns=['num_unique_peptides']).groupby(level=0, axis=1).mean()
    means = means.groupby(by=lambda x: x.split('.percent')[0], axis='columns').mean()
    result = add_meta(result, meta=df[['uniprot', 'symbol', 'description']])
    result = pd.merge(result, means, left_index=True, right_index=True)


    # ordering for supp table
    cols = list(result.columns)
    unique_conditions = list(SortedSet([next(iter(x.values())) for x in dataset.channel_layout]))
    first_cols = ['symbol', 'description', *unique_conditions]
    replicate_cols = [c for c in cols if c not in first_cols and c[0] != 'num_unique_peptides']
    num_peptide_cols = [c for c in cols if c not in first_cols and c[0] == 'num_unique_peptides']
    result = result[first_cols + num_peptide_cols + replicate_cols]

    result = result.rename(columns={x: f'{x} (mean)' for x in unique_conditions})
    result = result.rename(columns=dict(zip(replicate_cols, ['{} ({})'.format(*c) for c in replicate_cols])))
    result = result.rename(columns=dict(zip(num_peptide_cols, ['{} ({})'.format(*c) for c in num_peptide_cols])))

    result.index.rename('uniprot', inplace=True)

    result.to_excel(
        utils.get_timestamped_report_path(f'filtered_timecourse_{output_name}_{{}}.xlsx', DATA_OUTPUT_PATH),
        freeze_panes=(1,1),
        # index=False
    )

    return result

def filter_hydroxylamine_data(
    params_filename: str,
    output_name: str,
    output_psm_level=False,
    output_to_excel=True,
    min_num_unique_sequences=1,
    min_num_datasets=1
):
    analysis, params = analyze(params_filename, user='radu')
    dataset = analysis.datasets[0]
    df = analysis_to_df(analysis, dataset.channel_layout, dataset.control_channels)

    if output_psm_level:
        output_psm_level_report(df, f'unfiltered_{output_name}_psm_level')

    df = df[~df.index.isin(analysis.filtered_out)]
    
    df = df[~df.uniprot.str.startswith('contaminant_')]

    # processing for filtered protein level table and plotting
    def agg(x):
        result = x.filter(regex='percent_of_control').mean()
        result['num_unique_peptides'] = x['clean_sequence'].nunique()
        return result

    result = group_by_protein_and_filter(df, agg, min_num_unique_sequences, min_num_datasets)

    result.columns = pd.MultiIndex.from_tuples(
        (i.replace('.percent_of_control_', ' '), j) for i, j in result.columns
    )

    result.columns.set_levels([f'Replicate {i}' for i in range(1, len(dataset.experiments) + 1)], level=1, inplace=True)
    means = result.drop(columns=['num_unique_peptides']).groupby(level=0, axis=1).mean()

    result = add_meta(result, meta=df[['uniprot', 'symbol', 'description']])
    result = pd.merge(result, means, left_index=True, right_index=True)

    result['mean_reduction'] = 100 - means.filter(regex='Hydroxylamine').mean(axis=1)

    # ordering for supp table
    cols = list(result.columns)
    first_cols = ['symbol', 'description', 'mean_reduction', 'PBS 0', 'PBS 1', 'PBS 2', 'Hydroxylamine 0', 'Hydroxylamine 1', 'Hydroxylamine 2']
    replicate_cols = [c for c in cols if c not in first_cols and c[0] != 'num_unique_peptides']
    num_peptide_cols = [c for c in cols if c not in first_cols and c[0] == 'num_unique_peptides']
    result = result[first_cols + num_peptide_cols + replicate_cols]

    result = result.rename(columns=dict(zip(replicate_cols, ['{} ({})'.format(*c) for c in replicate_cols])))
    result = result.rename(columns=dict(zip(num_peptide_cols, ['{} ({})'.format(*c) for c in num_peptide_cols])))

    result.index.rename('uniprot', inplace=True)

    if output_to_excel:
        result.to_excel(
            utils.get_timestamped_report_path(f'filtered_hydroxylamine_{output_name}_{{}}.xlsx', DATA_OUTPUT_PATH),
            freeze_panes=(1,1),
            # index=False
        )

    return result


# oci_o = filter_timecourse_data(
#     'oci_on',
#     'resubmission_ocio_17odya.yaml',
#     set(['O15533']),
#     min_num_unique_sequences=2,
#     min_num_datasets=2,
#     output_psm_level=False,
#     output_to_excel=True
# )

# %%

def hydroxylamine():
    OUTPUT_PSM_LEVEL = True
    OUTPUT_TO_EXCEL = True
    MIN_NUM_UNIQUE_SEQUENCES = 2

    ha_nb4 = filter_hydroxylamine_data(
        'resubmission_nb4_nh2oh.yaml',
        'nb4',
        output_psm_level=OUTPUT_PSM_LEVEL,
        output_to_excel=OUTPUT_TO_EXCEL,
        min_num_unique_sequences=MIN_NUM_UNIQUE_SEQUENCES,
        min_num_datasets=1
    )
    ha_nb4 = ha_nb4[ha_nb4['mean_reduction'] >= HYDROXYLAMINE_REDUCTION_THRESHOLD]

    ha_oci = filter_hydroxylamine_data(
        'resubmission_nh2oh.yaml',
        'oci_p',
        output_psm_level=OUTPUT_PSM_LEVEL,
        output_to_excel=OUTPUT_TO_EXCEL,
        min_num_unique_sequences=MIN_NUM_UNIQUE_SEQUENCES,
        min_num_datasets=2
    )
    ha_oci = ha_oci[ha_oci['mean_reduction'] >= HYDROXYLAMINE_REDUCTION_THRESHOLD]

    ha_on = filter_hydroxylamine_data(
        'resubmission_ocio_nh2oh.yaml',
        'oci_o',
        output_psm_level=OUTPUT_PSM_LEVEL,
        output_to_excel=OUTPUT_TO_EXCEL,
        min_num_unique_sequences=MIN_NUM_UNIQUE_SEQUENCES,
        min_num_datasets=1
    )
    ha_on = ha_on[ha_on['mean_reduction'] >= HYDROXYLAMINE_REDUCTION_THRESHOLD]

    ha_proteins = set([
        *ha_oci.index.unique().values,
        *ha_on.index.unique().values,
        *ha_nb4.index.unique().values,
    ])

    from matplotlib_venn import venn3

    fig, ax = plt.subplots()
    fig.suptitle('Hydroxylamine Sensitive Proteins')

    venn3(subsets=(
        set(ha_oci.index.unique().values),
        set(ha_on.index.unique().values),
        set(ha_nb4.index.unique().values),
    ), set_labels = ('OCI-AML3', 'ON (75%)', 'NB-4'), ax=ax)

    fig.savefig('plots/hydroxylamine_venn.pdf')

    return ha_proteins


#%%

def timecourse(ha_proteins):
    OUTPUT_PSM_LEVEL = True
    OUTPUT_TO_EXCEL = True
    MIN_NUM_UNIQUE_SEQUENCES = 2

    oci = filter_timecourse_data(
        'oci_p',
        'resubmission_17odya.yaml',
        ha_proteins,
        min_num_unique_sequences=MIN_NUM_UNIQUE_SEQUENCES,
        min_num_datasets=2,
        output_psm_level=OUTPUT_PSM_LEVEL,
        output_to_excel=OUTPUT_TO_EXCEL
    )

    oci_on = filter_timecourse_data(
        'oci_o',
        'resubmission_ocio_17odya.yaml',
        ha_proteins,
        min_num_unique_sequences=MIN_NUM_UNIQUE_SEQUENCES,
        min_num_datasets=2,
        output_psm_level=OUTPUT_PSM_LEVEL,
        output_to_excel=OUTPUT_TO_EXCEL
    )

    nb4_10plex = filter_timecourse_data(
        'nb4',
        'resubmission_nb4_17odya.yaml',
        ha_proteins,
        min_num_unique_sequences=MIN_NUM_UNIQUE_SEQUENCES,
        min_num_datasets=1,
        output_psm_level=OUTPUT_PSM_LEVEL,
        output_to_excel=OUTPUT_TO_EXCEL
    )

    nb4_6plex = filter_timecourse_data(
        '6plex_nb4',
        'resubmission_nb4_17odya_6plex.yaml',
        ha_proteins,
        min_num_unique_sequences=MIN_NUM_UNIQUE_SEQUENCES,
        min_num_datasets=1,
        output_psm_level=OUTPUT_PSM_LEVEL,
        output_to_excel=OUTPUT_TO_EXCEL
    )

    #%%

    both = nb4_10plex.join(nb4_6plex, how='outer', lsuffix='_10plex', rsuffix='_6plex')
    both = both.drop(columns=list(both.filter(regex='mean')))
    # both = both.drop(columns=both.filter(regex='Palm').columns)
    means = both.filter(regex='percent').groupby(lambda x: x.split('.')[0], axis=1).mean()
    means.columns = [f'{c} (mean)' for c in list(means.columns)]
    counts = both.filter(regex='percent').groupby(lambda x: x.split('.')[0], axis=1).count()
    counts.columns = [f'{c} (count)' for c in list(counts.columns)]

    by_condition = both.filter(regex='percent').groupby(lambda x: x.split('.')[0], axis=1)
    stdevs = by_condition.std(ddof=0)
    stdevs.columns = [f'{c} (stdev)' for c in list(stdevs.columns)]

    both = both.join(means)
    both = both.join(counts)
    both = both.join(stdevs)

    both['symbol'] = both.symbol_6plex.fillna(both.symbol_10plex)
    both['description'] = both.description_6plex.fillna(both.description_10plex)


    both = both[both['ABD957 t1 (count)'].ge(4)]
    both = both[both['DMSO t1 (count)'].ge(4)]
    both = both[both['ABD957 t0 (count)'].ge(2)]
    both = both[both['DMSO t0 (count)'].ge(2)]
    both = both[~stdevs.ge(100).any(axis=1)]


    both = both[[
        'symbol',
        'description',
        'DMSO t0 (mean)',
        'Palm M t0 (mean)',
        'ABD957 t0 (mean)',
        'DMSO t1 (mean)',
        'ABD957 t1 (mean)',
        'Palm M t1 (mean)',
        # 'ABD957 t0 (count)',
        # 'ABD957 t1 (count)',
        # 'DMSO t0 (count)',
        # 'DMSO t1 (count)',
        # 'Palm M t0 (count)',
        # 'Palm M t1 (count)',
        # 'ABD957 t0 (stdev)',
        # 'ABD957 t1 (stdev)',
        # 'DMSO t0 (stdev)',
        # 'DMSO t1 (stdev)',
        # 'Palm M t0 (stdev)',
        # 'Palm M t1 (stdev)',
        'num_unique_peptides (Replicate 1)_10plex',
        'num_unique_peptides (Replicate 2)_10plex',
        'num_unique_peptides (Replicate 1)_6plex',
        'num_unique_peptides (Replicate 2)_6plex',
        'DMSO t0.percent_of_control_0 (Replicate 1)_10plex',
        'DMSO t0.percent_of_control_1 (Replicate 1)',
        'DMSO t0.percent_of_control_0 (Replicate 2)_10plex',
        'DMSO t0.percent_of_control_1 (Replicate 2)',
        'DMSO t0.percent_of_control_0 (Replicate 1)_6plex',
        'DMSO t0.percent_of_control_0 (Replicate 2)_6plex',
        'Palm M t0.percent_of_control_0 (Replicate 1)',
        'Palm M t0.percent_of_control_0 (Replicate 2)',
        'ABD957 t0.percent_of_control_0 (Replicate 1)_10plex',
        'ABD957 t0.percent_of_control_1 (Replicate 1)',
        'ABD957 t0.percent_of_control_0 (Replicate 2)_10plex',
        'ABD957 t0.percent_of_control_1 (Replicate 2)',
        'ABD957 t0.percent_of_control_0 (Replicate 1)_6plex',
        'ABD957 t0.percent_of_control_0 (Replicate 2)_6plex',
        'DMSO t1.percent_of_control_0 (Replicate 1)_10plex',
        'DMSO t1.percent_of_control_1 (Replicate 1)',
        'DMSO t1.percent_of_control_2 (Replicate 1)',
        'DMSO t1.percent_of_control_0 (Replicate 2)_10plex',
        'DMSO t1.percent_of_control_1 (Replicate 2)',
        'DMSO t1.percent_of_control_2 (Replicate 2)',
        'DMSO t1.percent_of_control_0 (Replicate 1)_6plex',
        'DMSO t1.percent_of_control_0 (Replicate 2)_6plex',
        'Palm M t1.percent_of_control_0 (Replicate 1)',
        'Palm M t1.percent_of_control_0 (Replicate 2)',
        'ABD957 t1.percent_of_control_0 (Replicate 1)_10plex',
        'ABD957 t1.percent_of_control_1 (Replicate 1)',
        'ABD957 t1.percent_of_control_2 (Replicate 1)',
        'ABD957 t1.percent_of_control_0 (Replicate 2)_10plex',
        'ABD957 t1.percent_of_control_1 (Replicate 2)',
        'ABD957 t1.percent_of_control_2 (Replicate 2)',
        'ABD957 t1.percent_of_control_0 (Replicate 1)_6plex',
        'ABD957 t1.percent_of_control_0 (Replicate 2)_6plex',
    ]]



    both_output_path = utils.get_timestamped_report_path('nb4_6plex_10plex_{}.xlsx', DATA_OUTPUT_PATH)
    both.to_excel(both_output_path)

#%%

def tabulate_timecourse(glob, name):
    time_course_path = utils.get_newest_file(DATA_OUTPUT_PATH, glob)

    df = pd.read_excel(time_course_path, index_col='uniprot')

    chase_col = 'DMSO t0/t1'
    palm_protect_col = 'Palm M t1/DMSO t1'
    abd_protect_col = 'ABD957 t1/DMSO t1'

    df[chase_col] = df['DMSO t0 (mean)']/df['DMSO t1 (mean)']
    df[palm_protect_col] = df['Palm M t1 (mean)']/df['DMSO t1 (mean)']
    df[abd_protect_col] = df['ABD957 t1 (mean)']/df['DMSO t1 (mean)']

    dynamic = df[chase_col] >= DYNAMIC_FOLD_CHANGE_THRESHOLD
    palm_protected = df[palm_protect_col] >= CHASE_FOLD_CHANGE_THRESHOLD
    abd_protected = df[abd_protect_col] >= CHASE_FOLD_CHANGE_THRESHOLD

    dynamic_colname = 'Dynamically palmitoylated?'

    df[dynamic_colname] = 'no'
    df.loc[dynamic, dynamic_colname] = 'yes'

    desired_cols = [
        'symbol', 'description',
        'DMSO t0 (mean)', 'Palm M t0 (mean)', 'ABD957 t0 (mean)',
        'DMSO t1 (mean)', 'Palm M t1 (mean)', 'ABD957 t1 (mean)',
        dynamic_colname,
    ]

    output_path = utils.get_timestamped_report_path(
        f'timecourse_{name.replace("-", "").lower()}_supp_table_extra_sheets_{{}}.xlsx',
        DATA_OUTPUT_PATH
    )

    with pd.ExcelWriter(output_path) as w:
        df[dynamic][desired_cols].drop(columns=[dynamic_colname]).to_excel(w, sheet_name=f'Dynamic ({name})')
        df[palm_protected][desired_cols].to_excel(w, sheet_name=f'Palm M Regulated ({name})')
        df[abd_protected][desired_cols].to_excel(w, sheet_name=f'ABD957 M Regulated ({name})')

def make_timecourse_extra_tables():
    tabulate_timecourse('nb4_6plex_10plex_*', 'NB-4')
    tabulate_timecourse('filtered_timecourse_oci_o_*', 'ON')
    tabulate_timecourse('filtered_timecourse_oci_p_*', 'OCI-AML3')

#%%

def get_hydroxylamine_sensitive_df(path, name):
    df = pd.read_excel(path, index_col='uniprot')
    df = df[df.mean_reduction >= HYDROXYLAMINE_REDUCTION_THRESHOLD]
    df = df[['symbol', 'description', 'mean_reduction']]
    df = df.rename(columns={c: f'{c} ({name})' for c in df.columns})
    return df

def tabulate_hydroxylamine():
    nb4_path = utils.get_newest_file(DATA_OUTPUT_PATH, 'filtered_hydroxylamine_nb4_*.xlsx')
    oci_path = utils.get_newest_file(DATA_OUTPUT_PATH, 'filtered_hydroxylamine_oci_p_*.xlsx')
    on_path = utils.get_newest_file(DATA_OUTPUT_PATH, 'filtered_hydroxylamine_oci_o_*.xlsx')
    oci = get_hydroxylamine_sensitive_df(oci_path, name='OCI-AML3')
    nb4 = get_hydroxylamine_sensitive_df(nb4_path, name='NB-4')
    on = get_hydroxylamine_sensitive_df(on_path, name='ON')

    df = oci.join(nb4, how='outer')
    df = df.join(on, how='outer')

    df['symbol'] = df['symbol (OCI-AML3)'].fillna(df['symbol (NB-4)']).fillna(df['symbol (ON)'])
    df['description'] = df['description (OCI-AML3)'].fillna(df['description (NB-4)']).fillna(df['description (ON)'])

    df = df[[
        'symbol',
        'description',
        'mean_reduction (OCI-AML3)',
        'mean_reduction (NB-4)',
        'mean_reduction (ON)',
    ]]

    df = df.sort_values(by='mean_reduction (OCI-AML3)', ascending=False).fillna('-')

    swisspalm = pd.read_excel('input/swisspalm_search.xlsx')

    swisspalm = swisspalm.set_index('Query identifier')
    swisspalm = swisspalm[[   
        # 'Query identifier',
        # 'UniProt AC',
        # 'UniProt ID',
        # 'UniProt status',
        # 'Organism',
        # 'Gene names',
        # 'Description',
        'Number of palmitoyl-proteomics articles',
        'Number of palmitoyl-proteomics studies where the protein appears in a high confidence hit list',
        'Number of technique categories used in palmitoyl-proteomics studies',
        'Technique categories used in palmitoyl-proteomics studies',
        'Number of targeted studies',
        'Targeted studies (PMIDs)',
        # 'PATs',
        # 'APTs',
        'Number of sites',
        'Sites in main isoform',
        'Number of isoforms',
        'Max number of cysteines',
        'Max number of cysteines in TM or cytosolic domain',
        'Predicted to be S-palmitoylated?',
        'Predicted to be S-palmitoylated in cytosolic domains?',
        'Protein has hits in SwissPalm?',
        'Orthologs of this protein have hits in SwissPalm?'
    ]]

    df = df.join(swisspalm)
    df = df.fillna('-')

    output_path = utils.get_timestamped_report_path('hydroxylamine_supp_table_extra_sheets_{}.xlsx', DATA_OUTPUT_PATH)

    with pd.ExcelWriter(output_path) as w:
        df.to_excel(w, sheet_name='Hydroxylamine Sensitive')


def main():
    ha_proteins = hydroxylamine()
    timecourse(ha_proteins)
    make_timecourse_extra_tables()

if __name__ == '__main__':
    main()

# %%


