#%%

from rifl import analyze
import rifl.utils as utils

import pandas as pd

import pathlib
import json



#%%
def sh_analysis(params_file: str, report_output_prefix: str):
    with open('input/human.json') as f:
        whitelist = json.loads(f.read())

    blacklist = [
        'P35030', 'P07477', 'P07478', # trypsins 
        'O15427', # SLC16A3
        'P00734', # F2
        'Q14703', # MBTSP1 (serine protease)
        'Q8NBP7', # PCSK9
    ]

    analysis, params = analyze(params_file, user='radu')

    m = analysis.data_model

    query = (m
        .select(
            m.id,
            m.experiment,
            m.uniprot,
            m.symbol,
            m.ratio,
        )
        .where(
            (m.experiment_id.in_(analysis.experiment_ids_included)) &
            (m.id.not_in(analysis.filtered_out))
        )
    )

    df = pd.DataFrame.from_records(list(query.dicts()))
    df = analysis.dataset_class.generate_id(df)

    for dataset in analysis.datasets:
        df.loc[df.experiment.isin(dataset.experiment_ids_included), 'condition'] = dataset.name

    df = df.set_index('id')

    df = df[df.uniprot.isin(whitelist)]
    df = df[~df.uniprot.isin(blacklist)]

    result_df = (df
        .groupby(['uniprot', 'symbol', 'condition', 'experiment'])
        .agg(ratio=('ratio', 'median'), num_peptides=('ratio', len))
        .groupby(level=('uniprot', 'symbol', 'condition'))  
        .agg(ratio=('ratio', 'median'), num_peptides=('num_peptides', 'sum'), ratio_list=('ratio', list))
    )

    # result_df['num_peptides'] = result_df.num_peptides.astype(int)
    result_df = result_df.unstack(level='condition')
    result_df['num_peptides'] = result_df.num_peptides.fillna(0).astype(int)
    # result_df['ratio'] = result_df.ratio.fillna('-')

    def ratios_to_string(ratios, invert=True):
        if not isinstance(ratios, list):
            return
        if invert:
            ratios = [1/r for r in ratios]
        return ', '.join(map(str, ratios))

    result_df['ratio_list'] = result_df.ratio_list.transform(lambda x: x.apply(ratios_to_string))
    result_df['ratio'] = result_df.ratio.rdiv(1)
    
    result_df.columns = result_df.columns.swaplevel()

    sorted_col_multiindex, _ = result_df.columns.sortlevel()
    result_df = result_df[sorted_col_multiindex]

    # result_df.loc[:, pd.IndexSlice[:, 'ratio']] = (result_df
    #     .loc[:, pd.IndexSlice[:, 'ratio']]
    #     .rdiv(1)
    # )

    result_df = (result_df
        .reset_index()
        .sort_values(by=[(analysis.datasets[0].name, 'ratio'), 'symbol'], ascending=True)
        .fillna('-')
    )

    report_output_path = utils.get_timestamped_report_path(
        f'{report_output_prefix}_{{}}.csv',
        pathlib.Path(analysis.params.output_folder),
    )

    result_df.to_csv(report_output_path, index=False, encoding='utf-8-sig')
    return analysis, result_df

insitu_analysis, insitu_df = sh_analysis('resubmission_insitu_abpp.yaml', 'in_situ')
invitro_analysis, invitro_df = sh_analysis('resubmission_invitro_abpp.yaml', 'in_vitro')
crisprko_analysis, crisprko_df = sh_analysis('resubmission_crisprko_abpp.yaml', 'crisprko')
#%%

crisprko_analysis.filter_report()
insitu_analysis.filter_report()
invitro_analysis.filter_report()

#%%

def unfiltered_report(analysis, report_output_prefix: str):
    m = analysis.data_model

    query = (m
        .select(
            m.id,
            m.experiment,
            m.uniprot,
            m.symbol,
            m.description,
            m.sequence,
            m.mass,
            m.charge,
            m.rsquared,
            m.ratio,
        )
        .where(
            m.experiment_id.in_(analysis.experiment_ids_included)
        )
    )

    df = pd.DataFrame.from_records(list(query.dicts()))

    for dataset in analysis.datasets:
        df.loc[df.experiment.isin(dataset.experiment_ids_included), 'condition'] = dataset.name

    df = df[~df.uniprot.str.startswith('Reverse_')]
    df['description'] = df.description.str.split().str[1:].str.join(' ')

    df = df.set_index(['uniprot', 'symbol', 'description', 'condition', 'experiment']).sort_index(level=0)

    report_output_path = utils.get_timestamped_report_path(
        f'unfiltered_{report_output_prefix}_{{}}.xlsx',
        pathlib.Path(analysis.params.output_folder),
    )

    df.to_excel(report_output_path, encoding='utf-8-sig')
    return df

# unfiltered_insitu_df = unfiltered_report(insitu_analysis, 'in_situ')
# unfiltered_invitro_df = unfiltered_report(invitro_analysis, 'in_vitro')
# unfiltered_crisprko_df = unfiltered_report(crisprko_analysis, 'crisprko')

#%%
