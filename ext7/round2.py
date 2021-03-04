#%%
from rifl import analyze

import pandas as pd


df = pd.read_excel(
    'Supplementary Table 3.xlsx',
    sheet_name='Hydroxylamine sensitive',
    skiprows=2,
)

uniprots_passing_ha_filter = set(df['Uniprot'].unique())

#%%

from resubmission import filter_timecourse_data

#%%

filter_timecourse_data(
    'dko1',
    'abhd17-dko1.yaml',
    uniprots_passing_ha_filter=uniprots_passing_ha_filter,
    min_num_unique_sequences=2,
    min_num_datasets=1,
    output_psm_level=True,
    output_to_excel=True
)

#%%

filter_timecourse_data(
    'parental_vs_dko1',
    'parental_vs_dko1.yaml',
    uniprots_passing_ha_filter=uniprots_passing_ha_filter,
    min_num_unique_sequences=2,
    min_num_datasets=1,
    output_psm_level=True,
    output_to_excel=True
)

#%%

import pandas as pd

palm_m_protected = [
    'NRAS',
    'SCRIB',
    'MPP6',
    'GNA12',
    'KCT2',
    'CYRIB',
    'LAPTM5',
    'ITM2C',
    'ADAM10',
    # 'HLA',
    'ASGR2',
    'MMP14',
    'MPP1',
    'PLPPR3',
    'NUP210',
    'DESI2',
    'NDFIP1',
    'ITM2B',
]

dko_10plex = pd.read_excel('output/filtered_timecourse_dko1_202101281504.xlsx')
dko_10plex = (dko_10plex
    .set_index('symbol')
    .filter(regex='DMSO.+percent|ABD957 t1\.percent')
    .rename(columns=lambda col: col.split('.')[0])
    .loc[palm_m_protected]
)

dko_16plex = pd.read_excel('output/filtered_timecourse_parental_vs_dko1_202101281508.xlsx')
dko_16plex = (dko_16plex
    .set_index('symbol')
    .filter([
        *dko_16plex.filter(regex='Parental DMSO.+percent|Parental ABD957 t1\.percent'),
        *dko_16plex.filter(regex='DKO1 DMSO.+percent.+1$|DKO ABD957 t1\.percent.+1$'),
    ])
    .rename(columns=lambda col: col.split('.')[0])
    .loc[palm_m_protected]
)

oci = pd.read_excel('Supplementary Table 4.xlsx', sheet_name='Time Course (OCI-AML3)', skiprows=2)
oci = (oci
    .set_index('Symbol')
    .filter(regex='DMSO.+Experiment|ABD957 t1 \(Experiment')
    .rename(columns=lambda col: col.split('(')[0])
    .loc[palm_m_protected]
)


dko = pd.concat([
    dko_10plex,
    dko_16plex.filter(regex='DKO').rename(columns=lambda col: col.replace('DKO1', 'DKO').replace('DKO ', ''))
], axis=1)

oci = pd.concat([
    oci,
    dko_16plex.filter(regex='Parental').rename(columns=lambda col: col.replace('Parental ', ''))
], axis=1)


dko_dynamic = (dko
    .filter(regex='DMSO t1')
    .transform(lambda x: x.div(dko.filter(regex='DMSO t0').mean(axis=1)))
)
dko_protected = (dko
    .filter(regex='ABD957 t1')
    .transform(lambda x: x.div(dko.filter(regex='DMSO t1').mean(axis=1)))
)
oci_dynamic = (oci
    .filter(regex='DMSO t1')
    .transform(lambda x: x.div(oci.filter(regex='DMSO t0').mean(axis=1)))
)
oci_protected = (oci
    .filter(regex='ABD957 t1')
    .transform(lambda x: x.div(oci.filter(regex='DMSO t1').mean(axis=1)))
)

protected = pd.concat([oci_protected, dko_protected], axis=1)
dynamic = pd.concat([oci_dynamic, dko_dynamic], axis=1)

#%%

from rifl import utils
import pathlib

dko_10plex = pd.read_excel('output/filtered_timecourse_dko1_202101281504.xlsx', index_col='uniprot')
dko_16plex = pd.read_excel('output/filtered_timecourse_parental_vs_dko1_202101281508.xlsx', index_col='uniprot')

dko_10plex = dko_10plex.rename(columns=lambda x: f'DKO1 {x}' if 'percent' in x else x)
dko_10plex = dko_10plex.rename(columns=lambda x: f'{x}_10plex')

dko_16plex = dko_16plex[['symbol', 'description', 'DKO ABD957 t0 (mean)', 'DKO ABD957 t1 (mean)',
       'DKO1 DMSO t0 (mean)', 'DKO1 DMSO t1 (mean)',
       'Parental ABD957 t0 (mean)', 'Parental ABD957 t1 (mean)',
       'Parental DMSO t0 (mean)', 'Parental DMSO t1 (mean)',
       'num_unique_peptides (Replicate 1)',
       'Parental DMSO t0.percent_of_control_0 (Replicate 1)',
       'Parental DMSO t0.percent_of_control_1 (Replicate 1)',
       'DKO1 DMSO t0.percent_of_control_0 (Replicate 1).1',
       'DKO1 DMSO t0.percent_of_control_1 (Replicate 1).1',
       'Parental ABD957 t0.percent_of_control_0 (Replicate 1)',
       'Parental ABD957 t0.percent_of_control_1 (Replicate 1)',
       'DKO ABD957 t0.percent_of_control_0 (Replicate 1).1',
       'DKO ABD957 t0.percent_of_control_1 (Replicate 1).1',
       'Parental DMSO t1.percent_of_control_0 (Replicate 1)',
       'Parental DMSO t1.percent_of_control_1 (Replicate 1)',
       'DKO1 DMSO t1.percent_of_control_0 (Replicate 1).1',
       'DKO1 DMSO t1.percent_of_control_1 (Replicate 1).1',
       'Parental ABD957 t1.percent_of_control_0 (Replicate 1)',
       'Parental ABD957 t1.percent_of_control_1 (Replicate 1)',
       'DKO ABD957 t1.percent_of_control_0 (Replicate 1).1',
       'DKO ABD957 t1.percent_of_control_1 (Replicate 1).1',
]]

dko_16plex = dko_16plex.rename(columns=lambda x: x.replace(').1', ')') if x.endswith('.1') else x)
dko_16plex = dko_16plex.rename(columns=lambda x: x.replace('DKO', 'DKO1') if x.startswith('DKO ') else x)
dko_16plex = dko_16plex.rename(columns=lambda x: f'{x}_16plex')

both = dko_10plex.join(dko_16plex, how='outer')

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

both['symbol'] = both.symbol_16plex.fillna(both.symbol_10plex)
both['description'] = both.description_16plex.fillna(both.description_10plex)


# both = both[both['ABD957 t1 (count)'].ge(4)]
# both = both[both['DMSO t1 (count)'].ge(4)]
# both = both[both['ABD957 t0 (count)'].ge(2)]
# both = both[both['DMSO t0 (count)'].ge(2)]


# both = both[~stdevs.ge(100).any(axis=1)]
both = both[[
       'symbol',
       'description',
       'num_unique_peptides (Replicate 1)_10plex',
       'num_unique_peptides (Replicate 1)_16plex',
       'DKO1 DMSO t0 (mean)',
       'Parental DMSO t0 (mean)',
       'DKO1 ABD957 t0 (mean)',
       'Parental ABD957 t0 (mean)',
       'DKO1 DMSO t1 (mean)',
       'Parental DMSO t1 (mean)',
       'DKO1 ABD957 t1 (mean)',
       'Parental ABD957 t1 (mean)',
       'DKO1 DMSO t0.percent_of_control_0 (Replicate 1)_10plex',
       'DKO1 DMSO t0.percent_of_control_1 (Replicate 1)_10plex',
       'DKO1 DMSO t0.percent_of_control_0 (Replicate 1)_16plex',
       'DKO1 DMSO t0.percent_of_control_1 (Replicate 1)_16plex',
       'Parental DMSO t0.percent_of_control_0 (Replicate 1)_16plex',
       'Parental DMSO t0.percent_of_control_1 (Replicate 1)_16plex',
       'DKO1 ABD957 t0.percent_of_control_0 (Replicate 1)_10plex',
       'DKO1 ABD957 t0.percent_of_control_1 (Replicate 1)_10plex',
       'DKO1 ABD957 t0.percent_of_control_0 (Replicate 1)_16plex',
       'DKO1 ABD957 t0.percent_of_control_1 (Replicate 1)_16plex',
       'Parental ABD957 t0.percent_of_control_0 (Replicate 1)_16plex',
       'Parental ABD957 t0.percent_of_control_1 (Replicate 1)_16plex',
       'DKO1 DMSO t1.percent_of_control_0 (Replicate 1)_10plex',
       'DKO1 DMSO t1.percent_of_control_1 (Replicate 1)_10plex',
       'DKO1 DMSO t1.percent_of_control_2 (Replicate 1)_10plex',
       'DKO1 DMSO t1.percent_of_control_0 (Replicate 1)_16plex',
       'DKO1 DMSO t1.percent_of_control_1 (Replicate 1)_16plex',
       'Parental DMSO t1.percent_of_control_0 (Replicate 1)_16plex',
       'Parental DMSO t1.percent_of_control_1 (Replicate 1)_16plex',
       'DKO1 ABD957 t1.percent_of_control_0 (Replicate 1)_10plex',
       'DKO1 ABD957 t1.percent_of_control_1 (Replicate 1)_10plex',
       'DKO1 ABD957 t1.percent_of_control_2 (Replicate 1)_10plex',
       'DKO1 ABD957 t1.percent_of_control_0 (Replicate 1)_16plex',
       'DKO1 ABD957 t1.percent_of_control_1 (Replicate 1)_16plex',
       'Parental ABD957 t1.percent_of_control_0 (Replicate 1)_16plex',
       'Parental ABD957 t1.percent_of_control_1 (Replicate 1)_16plex',
]]

both = both.rename(columns=lambda x: x.replace('Replicate 1)_10plex', 'Experiment 1)_10plex'))
both = both.rename(columns=lambda x: x.replace('Replicate 1)_16plex', 'Experiment 2)_16plex'))


both_output_path = utils.get_timestamped_report_path('dko_16plex_10plex_{}.xlsx', pathlib.Path('output'))
both.to_excel(both_output_path)


# %%
