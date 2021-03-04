#%%
from scipy.linalg.special_matrices import dft
from rifl import analyze
import rifl.filters.isobaric as filters

from scipy.optimize import curve_fit
from scipy.stats.distributions import t as t_distribution
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np

from multiprocessing import Pool
import itertools
import inspect

RED_COLOR = sns.color_palette('tab10')[3]
BLUE_COLOR = sns.color_palette('tab10')[0]
GREY_COLOR = sns.color_palette('tab10')[7]



# pd.set_option('use_inf_as_na', True)

analysis, params = analyze('hyperplex2.yaml', user='radu_hyperplex2')
combined = analysis.datasets[0]

required_columns = ['id', 'experiment', 'uniprot', 'symbol', 'description', 'sequence', 'unique_sequence', 'channel_intensities', 'filename', 'scan_num']
required_columns = [getattr(combined.data_model, col) for col in required_columns]

combined_df = combined.as_df(required_columns)
combined_df = combined_df[combined_df.unique_sequence]
combined_df = combined_df[~combined_df.sequence.str.contains('-')]
combined_df = combined_df[~combined_df.uniprot.str.startswith('Reverse')]
combined_df = combined_df[~combined_df.uniprot.str.startswith('contaminant')]
combined_df = combined_df[combined_df.sequence.apply(filters.is_not_half_tryptic.fn)]

# TODO: reconsider intensity filters
# heavy_df = heavy_df[heavy_df['heavy.intensity_0'].ge(5000)]
# light_df = light_df[light_df['light.intensity_0'].ge(5000)]
heavy_df = combined_df[combined_df.filename.str.startswith('H')]
light_df = combined_df[~combined_df.filename.str.startswith('H')]
heavy_df = heavy_df[heavy_df['heavy.intensity_0'].gt(0)]
light_df = light_df[light_df['light.intensity_0'].gt(0)]

def get_normalization(light_df, heavy_df):
    heavy_df = heavy_df.copy().groupby(['sequence', 'experiment']).agg(
        uniprot=('uniprot', 'first'),
        symbol=('symbol', 'first'),
        description=('description', 'first'),
        **{col: (col, 'sum') for col in list(heavy_df.filter(regex='intensity').columns)}
    )
    light_df = light_df.groupby(['sequence', 'experiment']).agg(
        uniprot=('uniprot', 'first'),
        symbol=('symbol', 'first'),
        description=('description', 'first'),
        **{col: (col, 'sum') for col in list(light_df.filter(regex='intensity').columns)}
    )
    shared_sequences = heavy_df.index.intersection(light_df.index)

    heavy_df = heavy_df[heavy_df.index.isin(shared_sequences)]
    light_df = light_df[light_df.index.isin(shared_sequences)]

    light_df['nrow'] = (light_df['light.intensity_0'] + light_df['heavy.intensity_0'])/(2*light_df['light.intensity_0'])
    heavy_df['nrow'] = (heavy_df['light.intensity_0'] + heavy_df['heavy.intensity_0'])/(2*heavy_df['heavy.intensity_0'])

    light_df[light_df.filter(regex='.intensity').columns] = light_df.filter(regex='.intensity').mul(light_df.nrow, axis=0)
    heavy_df[heavy_df.filter(regex='.intensity').columns] = heavy_df.filter(regex='.intensity').mul(heavy_df.nrow, axis=0)

    combined_again = light_df.append(heavy_df)

    # TODO: reconsider the fenced zero filter
    # fenced_zeros_dmso = combined_df.filter(regex='dmso.+intensity').apply(
    #     lambda x: np.sum(k for k, _ in itertools.groupby(np.trim_zeros(x) == 0)) > 0,
    #     axis=1,
    #     raw=True
    # )
    # fenced_zeros_957 = combined_df.filter(regex='957.+intensity').apply(
    #     lambda x: np.sum(k for k, _ in itertools.groupby(np.trim_zeros(x) == 0)) > 0,
    #     axis=1,
    #     raw=True
    # )

    combined_sum = combined_again.filter(regex='intensity').sum(level=1, axis=0)
    median_channel_sum = combined_sum.median(axis=1)

    nsum = combined_sum.median(axis=1).div(combined_sum)

    return light_df.nrow, heavy_df.nrow, median_channel_sum.div(combined_sum.T)


# get the row and channel normalization values
light_nrow, heavy_nrow, nsum = get_normalization(light_df, heavy_df)

# throw out any sequence that is not in both heavy and light
# note this is done in the process of normalization
light_df = light_df.set_index(['sequence', 'experiment'])
light_df = light_df[light_df.index.isin(light_nrow.index)]
heavy_df = heavy_df.set_index(['sequence', 'experiment'])
heavy_df = heavy_df[heavy_df.index.isin(heavy_nrow.index)]


def plot_channel_sums(light_df, heavy_df, title):
    light = light_df.filter(regex='intensity').sum()
    heavy = heavy_df.filter(regex='intensity').sum()

    palette = sns.color_palette('tab10')
    ind = np.arange(len(light))

    fig, ax = plt.subplots()
    ax.bar(ind, heavy, color=palette[0])
    ax.bar(ind, light, bottom=heavy, color=palette[3])

    ax.set_ylabel('Summed TMT Intensity')
    ax.set_xlabel('TMT channel')
    ax.set_xticks(ind)
    ax.set_xticklabels(ind + 1)

    fig.suptitle(title)

plot_channel_sums(light_df, heavy_df, 'Raw data')

light_df = light_df.join(light_nrow)
heavy_df = heavy_df.join(heavy_nrow)

# apply the normalization back
light_df[light_df.filter(regex='.intensity').columns] = light_df.filter(regex='.intensity').mul(light_df.nrow, axis=0)
heavy_df[heavy_df.filter(regex='.intensity').columns] = heavy_df.filter(regex='.intensity').mul(heavy_df.nrow, axis=0)

light_df[light_df.filter(regex='.intensity').columns] = light_df.filter(regex='.intensity').mul(nsum.T, level=1)
heavy_df[heavy_df.filter(regex='.intensity').columns] = heavy_df.filter(regex='.intensity').mul(nsum.T, level=1)

plot_channel_sums(light_df, heavy_df, 'Normalized data')

light_df = light_df.reset_index()
heavy_df = heavy_df.reset_index()

light_perc = light_df.filter(regex='intensity').div(light_df['light.intensity_0'], axis=0)
heavy_perc = heavy_df.filter(regex='intensity').div(heavy_df['heavy.intensity_0'], axis=0)

light_perc = light_perc.rename(columns={
    c: c.replace('intensity', 'ratio')
    for c in light_df.filter(regex='\.intensity').columns
})

heavy_perc = heavy_perc.rename(columns={
    c: c.replace('intensity', 'ratio')
    for c in heavy_df.filter(regex='\.intensity').columns
})

light_df = light_df.join(light_perc)
heavy_df = heavy_df.join(heavy_perc)

heavy_df = heavy_df[heavy_df['heavy.intensity_0'].gt(2500)]
light_df = light_df[light_df['light.intensity_0'].gt(2500)]

# pd.set_option('use_inf_as_na', True)


#%%

def synthesis(t, a_syn, b_syn, k_syn):
    return (b_syn - a_syn) * np.exp(-k_syn*t) + a_syn

def degradation(t, a_deg, b_deg, k_deg):
    return (a_deg - b_deg) * np.exp(-k_deg*t) + b_deg

dmso_data_cols = list(light_df.filter(regex='light.+ratio|dmso.+ratio'))
abd_data_cols = list(light_df.filter(regex='light.+ratio|957.+ratio'))

timepoints = np.array([0, 3, 6, 24, 72])

fit_tasks = [(
        light_df,
        dmso_data_cols,
        degradation,
        ['dmso_a_deg', 'dmso_b_deg', 'dmso_k_deg', 'dmso_rsq_deg'],
        [0.9, 0.1, 0.04]
    ), (
        light_df,
        abd_data_cols,
        degradation,
        ['abd_a_deg', 'abd_b_deg', 'abd_k_deg', 'abd_rsq_deg'],
        [0.9, 0.1, 0.04]
    ), (
        heavy_df,
        dmso_data_cols,
        synthesis,
        ['dmso_a_syn', 'dmso_b_syn', 'dmso_k_syn', 'dmso_rsq_syn'],
        [0.9, 0.1, 0.04]
    ), (
        heavy_df,
        abd_data_cols,
        synthesis,
        ['abd_a_syn', 'abd_b_syn', 'abd_k_syn', 'abd_rsq_syn'],
        [0.9, 0.1, 0.04]
    )
]


# %% 

bounds = ((0, 0, 0), (1, 1, 10))

def fit_curves_psm(g, fn, arg_names, p0):
    vals = g.dropna()
    output = None

    try:
        x = timepoints
        y = vals.values.flatten()

        popt, pcov = curve_fit(fn, x, y, p0=p0, bounds=bounds)

        model_predictions = fn(x, *popt)
        abs_error = model_predictions - y

        r_squared = 1 - (np.var(abs_error) / np.var(y))
        output = [*popt, r_squared]

    except (ValueError, RuntimeError) as e:
        output = [np.nan] * len(arg_names)

    return pd.Series(dict(zip(arg_names, output)))


def get_psm_fit(df, data_cols, fn, arg_names, p0):
    return df[data_cols].apply(fit_curves_psm, axis=1, fn=fn, arg_names=arg_names, p0=p0)

#%%

with Pool(20) as p:
    res = p.starmap(get_psm_fit, fit_tasks)
    dmso_deg_opt_psm, abd_deg_opt_psm, dmso_syn_opt_psm, abd_syn_opt_psm = res


#%%
def get_psm_sheet(df):
    psm_df = df[[
        'uniprot',
        'symbol',
        'description',
        'sequence',
        'experiment',
        'filename',
        'scan_num',
        'light.intensity_0',
        'dmso_3h.intensity_0',
        'dmso_6h.intensity_0',
        'dmso_24h.intensity_0',
        'dmso_72h.intensity_0',
        '957_3h.intensity_0',
        '957_6h.intensity_0',
        '957_24h.intensity_0',
        '957_72h.intensity_0',
        'heavy.intensity_0',
        'light.ratio_0',
        'dmso_3h.ratio_0',
        'dmso_6h.ratio_0',
        'dmso_24h.ratio_0',
        'dmso_72h.ratio_0',
        '957_3h.ratio_0',
        '957_6h.ratio_0',
        '957_24h.ratio_0',
        '957_72h.ratio_0',
        'heavy.ratio_0'
    ]]

    psm_df = psm_df.rename(
        columns={col: col.replace('_0', '') for col in list(psm_df.columns)}
    )
    psm_df.index = psm_df.index.rename('id')
    return psm_df

light_psm_df = get_psm_sheet(light_df)
heavy_psm_df = get_psm_sheet(heavy_df)

light_psm_df.to_excel('output/hyperplex_light_psms.xlsx')
heavy_psm_df.to_excel('output/hyperplex_heavy_psms.xlsx')

# for debugging in Proturn
light_psm_df[light_psm_df.symbol == 'NRAS'].to_csv('output/hyperplex_light_nras_psms.tsv', sep='\t')
heavy_psm_df[heavy_psm_df.symbol == 'NRAS'].to_csv('output/hyperplex_heavy_nras_psms.tsv', sep='\t')

#%%

def krab_filter(df, a=(0.67, 1.5), b=(0, 0.3), k=(0, 5), rsq=0.7):
    df = df[df.filter(regex='_a_').iloc[:, 0].between(a[0], a[1], inclusive=True)]
    df = df[df.filter(regex='_b_').iloc[:, 0].between(b[0], b[1], inclusive=True)]
    df = df[df.filter(regex='_k_').iloc[:, 0].between(k[0], k[1], inclusive=True)]
    df = df[df.filter(regex='_rsq_').iloc[:, 0].ge(rsq)]
    return df

dmso_deg_opt_filtered = krab_filter(dmso_deg_opt_psm)
abd_deg_opt_filtered = krab_filter(abd_deg_opt_psm)
dmso_syn_opt_filtered = krab_filter(dmso_syn_opt_psm)
abd_syn_opt_filtered = krab_filter(abd_syn_opt_psm)

light_psm_df = light_psm_df.loc[dmso_deg_opt_filtered.index.intersection(abd_deg_opt_filtered.index)]
heavy_psm_df = heavy_psm_df.loc[dmso_syn_opt_filtered.index.intersection(abd_syn_opt_filtered.index)]

light_psm_df[light_psm_df.symbol == 'NRAS'].to_csv('output/hyperplex_light_nras_filtered_psms.tsv', sep='\t')
heavy_psm_df[heavy_psm_df.symbol == 'NRAS'].to_csv('output/hyperplex_heavy_nras_filtered_psms.tsv', sep='\t')

#%%

dmso_data_cols = list(light_psm_df.filter(regex='light.+ratio|dmso.+ratio'))
abd_data_cols = list(light_psm_df.filter(regex='light.+ratio|957.+ratio'))

fit_tasks = [(
        light_psm_df,
        dmso_data_cols,
        degradation,
        ['dmso_a', 'dmso_b', 'dmso_k', 'dmso_rsq', 'dmso_ci95'],
        [0.9, 0.1, 0.04]
    ), (
        light_psm_df,
        abd_data_cols,
        degradation,
        ['abd_a', 'abd_b', 'abd_k', 'abd_rsq', 'abd_ci95'],
        [0.9, 0.1, 0.04]
    ), (
        heavy_psm_df,
        dmso_data_cols,
        synthesis,
        ['dmso_a', 'dmso_b', 'dmso_k', 'dmso_rsq', 'dmso_ci95'],
        [0.9, 0.1, 0.04]
    ), (
        heavy_psm_df,
        abd_data_cols,
        synthesis,
        ['abd_a', 'abd_b', 'abd_k', 'abd_rsq', 'abd_ci95'],
        [0.9, 0.1, 0.04]
    )
]

bounds = ((0, 0, 0), (1, 1, 10))

def fit_curves(g, fn, arg_names, p0):
    vals = g.dropna()
    output = None

    try:
        x = np.tile(timepoints, len(vals))
        y = vals.values.flatten()

        popt, pcov = curve_fit(fn, x, y, p0=p0, bounds=bounds)
        # method of extracting parameter CI values taken from
        # http://kitchingroup.cheme.cmu.edu/blog/2013/02/12/Nonlinear-curve-fitting-with-parameter-confidence-intervals/
        # other ref:
        # https://stackoverflow.com/a/60412600/383744

        alpha = 0.05 # 95% confidence interval = 100*(1-alpha)
        n = len(y) # number of data points
        p = len(popt) # number of parameters
        dof = max(0, n - p) # number of degrees of freedom

        tval = t_distribution.ppf(1.0-alpha/2., dof)

        var = np.diag(pcov)
        sigma_k = np.sqrt(var[-1])
        ci_k = tval * sigma_k

        model_predictions = fn(x, *popt)
        abs_error = model_predictions - y

        r_squared = 1 - (np.var(abs_error) / np.var(y))
        output = [*popt, r_squared, ci_k]

    except (ValueError, RuntimeError) as e:
        output = [np.nan] * len(arg_names)

    return pd.Series(dict(zip(arg_names, output)))


def get_fit(df, data_cols, fn, arg_names, p0):
    return (df
        .groupby('uniprot')[data_cols]
        .apply(fit_curves, fn=fn, arg_names=arg_names, p0=p0)
    )


with Pool(8) as p:
    res = p.starmap(get_fit, fit_tasks)
    dmso_deg_opt, abd_deg_opt, dmso_syn_opt, abd_syn_opt = res


#%%

light_grouped = light_psm_df.groupby('uniprot').agg({
    'experiment': set,
    'symbol': 'first',
    'description': 'first',
    'sequence': set,
    **{ col: ['mean', list] for col in light_psm_df.filter(regex='ratio')}
}).join([dmso_deg_opt, abd_deg_opt])

heavy_grouped = heavy_psm_df.groupby('uniprot').agg({
    'experiment': set,
    'symbol': 'first',
    'description': 'first',
    'sequence': set,
    **{ col: ['mean', list] for col in heavy_psm_df.filter(regex='ratio')}
}).join([dmso_syn_opt, abd_syn_opt])

# light_grouped['curve'] = 'degradation'
# heavy_grouped['curve'] = 'synthesis'

cell_counts = pd.read_excel('time_course_experiment_growth_data.xlsx', sheet_name='analysis')
initial_counts = pd.read_excel('time_course_experiment_growth_data.xlsx', sheet_name='initial_counts')


initial = initial_counts.merge(cell_counts[['Sample #', 'Treatment']], left_on='Unnamed: 0', right_on='Sample #')
initial_dmso = initial[initial['Treatment'] == 'DMSO'].filter(regex='Cell count').mean().mean()
initial_abd = initial[initial['Treatment'] == 'ABD957'].filter(regex='Cell count').mean().mean()

times, dmso_counts = (cell_counts[cell_counts['Treatment'] == 'DMSO']
    .filter(regex='Time|Cell count')
    .set_index('Time (hrs)')
    .unstack()
    .reset_index()
    .drop(columns=['level_0'])
    .T
    .values
)

dmso_opt, _ = curve_fit(lambda t, k: initial_dmso*np.exp(k*t), times, dmso_counts)
dmso_doubling_time = np.log(2)/dmso_opt

times, abd_counts = (cell_counts[cell_counts['Treatment'] == 'ABD957']
    .filter(regex='Time|Cell count')
    .set_index('Time (hrs)')
    .unstack()
    .reset_index()
    .drop(columns=['level_0'])
    .T
    .values
)

abd_opt, _ = curve_fit(lambda t, k: initial_abd*np.exp(k*t), times, abd_counts)
abd_doubling_time = np.log(2)/abd_opt


light_grouped['dmso_k_corr'] = light_grouped['dmso_k'].sub(dmso_opt[0])
light_grouped['abd_k_corr'] = light_grouped['abd_k'].sub(abd_opt[0])
heavy_grouped['dmso_k_corr'] = heavy_grouped['dmso_k'].sub(dmso_opt[0])
heavy_grouped['abd_k_corr'] = heavy_grouped['abd_k'].sub(abd_opt[0])

light_grouped['dmso_halflife'] = np.log(2)/light_grouped['dmso_k_corr']
light_grouped['abd_halflife'] = np.log(2)/light_grouped['abd_k_corr']
heavy_grouped['dmso_halflife'] = np.log(2)/heavy_grouped['dmso_k_corr']
heavy_grouped['abd_halflife'] = np.log(2)/heavy_grouped['abd_k_corr']

# lol
light_grouped['dmso_halflife_ci95'] = (light_grouped['dmso_ci95']
    .mul(np.log(2))
    .div(
        light_grouped['dmso_k_corr'].sub(light_grouped['dmso_ci95'])
        .mul(light_grouped['dmso_k_corr'].add(light_grouped['dmso_ci95']))
    )
)

light_grouped['abd_halflife_ci95'] = (light_grouped['abd_ci95']
    .mul(np.log(2))
    .div(
        light_grouped['abd_k_corr'].sub(light_grouped['abd_ci95'])
        .mul(light_grouped['abd_k_corr'].add(light_grouped['abd_ci95']))
    )
)

heavy_grouped['dmso_halflife_ci95'] = (heavy_grouped['dmso_ci95']
    .mul(np.log(2))
    .div(
        heavy_grouped['dmso_k_corr'].sub(heavy_grouped['dmso_ci95'])
        .mul(heavy_grouped['dmso_k_corr'].add(heavy_grouped['dmso_ci95']))
    )
)

heavy_grouped['abd_halflife_ci95'] = (heavy_grouped['abd_ci95']
    .mul(np.log(2))
    .div(
        heavy_grouped['abd_k_corr'].sub(heavy_grouped['abd_ci95'])
        .mul(heavy_grouped['abd_k_corr'].add(heavy_grouped['abd_ci95']))
    )
)

#%%

def krab_filter(df, a=(0.67, 1.5), b=(0, 0.3), k=(0, 5), rsq=0.7):
    df = df[df.filter(regex='_a').transform(lambda x: x.between(a[0], a[1], inclusive=True)).all(axis=1)]
    df = df[df.filter(regex='_b').transform(lambda x: x.between(b[0], b[1], inclusive=True)).all(axis=1)]
    df = df[df.filter(regex='_k').transform(lambda x: x.between(k[0], k[1], inclusive=True)).all(axis=1)]
    df = df[df.filter(regex='_rsq').transform(lambda x: x.ge(rsq)).all(axis=1)]
    return df

light_grouped = krab_filter(light_grouped)
heavy_grouped = krab_filter(heavy_grouped)


#%%

def get_protein_df(df):
    protein_df = df[[
        ('symbol', 'first'),
        ('description', 'first'),
        ('experiment', 'set'),
        ('sequence', 'set'),
        ('light.ratio', 'mean'),
        ('dmso_3h.ratio', 'mean'),
        ('dmso_6h.ratio', 'mean'),
        ('dmso_24h.ratio', 'mean'),
        ('dmso_72h.ratio', 'mean'),
        ('957_3h.ratio', 'mean'),
        ('957_6h.ratio', 'mean'),
        ('957_24h.ratio', 'mean'),
        ('957_72h.ratio', 'mean'),
        ('heavy.ratio', 'mean'),
        ('light.ratio', 'list'),
        ('dmso_3h.ratio', 'list'),
        ('dmso_6h.ratio', 'list'),
        ('dmso_24h.ratio', 'list'),
        ('dmso_72h.ratio', 'list'),
        ('957_3h.ratio', 'list'),
        ('957_6h.ratio', 'list'),
        ('957_24h.ratio', 'list'),
        ('957_72h.ratio', 'list'),
        ('heavy.ratio', 'list'),
        'dmso_a',
        'dmso_b',
        'dmso_k',
        'dmso_k_corr',
        'dmso_ci95',
        'dmso_rsq',
        'abd_a',
        'abd_b',
        'abd_k',
        'abd_k_corr',
        'abd_ci95',
        'abd_rsq',
        # 'curve',
        'dmso_halflife',
        'dmso_halflife_ci95',
        'abd_halflife',
        'abd_halflife_ci95'
    ]]

    protein_df = protein_df.rename(columns={
        ('experiment', 'set'): 'experiment',
        ('symbol', 'first'): 'symbol',
        ('description', 'first'): 'description',
        ('sequence', 'set'): 'sequences',
        ('light.ratio', 'mean'): 'light.ratio (mean)',
        ('dmso_3h.ratio', 'mean'): 'dmso_3h.ratio (mean)',
        ('dmso_6h.ratio', 'mean'): 'dmso_6h.ratio (mean)',
        ('dmso_24h.ratio', 'mean'): 'dmso_24h.ratio (mean)',
        ('dmso_72h.ratio', 'mean'): 'dmso_72h.ratio (mean)',
        ('957_3h.ratio', 'mean'): '957_3h.ratio (mean)',
        ('957_6h.ratio', 'mean'): '957_6h.ratio (mean)',
        ('957_24h.ratio', 'mean'): '957_24h.ratio (mean)',
        ('957_72h.ratio', 'mean'): '957_72h.ratio (mean)',
        ('heavy.ratio', 'mean'): 'heavy.ratio (mean)',
        ('light.ratio', 'list'): 'light.ratio (list)',
        ('dmso_3h.ratio', 'list'): 'dmso_3h.ratio (list)',
        ('dmso_6h.ratio', 'list'): 'dmso_6h.ratio (list)',
        ('dmso_24h.ratio', 'list'): 'dmso_24h.ratio (list)',
        ('dmso_72h.ratio', 'list'): 'dmso_72h.ratio (list)',
        ('957_3h.ratio', 'list'): '957_3h.ratio (list)',
        ('957_6h.ratio', 'list'): '957_6h.ratio (list)',
        ('957_24h.ratio', 'list'): '957_24h.ratio (list)',
        ('957_72h.ratio', 'list'): '957_72h.ratio (list)',
        ('heavy.ratio', 'list'): 'heavy.ratio (list)',
        'dmso_a': 'A (DMSO)',
        'dmso_b': 'B (DMSO)',
        'dmso_k': 'k (DMSO)',
        'dmso_k_corr': 'k (DMSO, cell doubling corrected)',
        'dmso_ci95': 'k ci95 (DMSO)',
        'dmso_rsq': 'R_squared (DMSO)',
        'abd_a': 'A (ABD957)',
        'abd_b': 'B (ABD957)',
        'abd_k': 'k (ABD957)',
        'abd_ci95': 'k ci95 (ABD957)',
        'abd_k_corr': 'k (ABD957, cell doubling corrected)',
        'abd_rsq': 'R_squared (ABD957)',
        # 'curve': '',
        'dmso_halflife': 'Halflife (DMSO, hours)',
        'dmso_halflife_ci95': 'Halflife margin of error at 95% CI (DMSO, hours)',
        'abd_halflife': 'Halflife (ABD957, hours)',
        'abd_halflife_ci95': 'Halflife margin of error at 95% CI (ABD957, hours)',
    })

    def list_to_comma_delimited(lst):
        return ', '.join(map(str, lst))

    protein_df['experiment'] = protein_df.experiment.apply(list_to_comma_delimited)
    protein_df['sequences'] = protein_df.sequences.apply(list_to_comma_delimited)
    protein_df[protein_df.filter(regex='list').columns] =  protein_df.filter(regex='list').apply(lambda x: x.apply(list_to_comma_delimited))

    return protein_df

light_protein_df = get_protein_df(light_grouped)
heavy_protein_df = get_protein_df(heavy_grouped)

light_protein_df.to_excel('output/hyperplex_light_proteins.xlsx')
heavy_protein_df.to_excel('output/hyperplex_heavy_proteins.xlsx')



#%%

plt.rcParams.update({
    # 'font.weight': 'normal',
    # 'axes.labelweight': 'normal',
    # 'axes.titleweight': 'normal',
    'font.sans-serif': 'Arial',
    'font.family': 'sans-serif',
    # # 'mathtext.fontset': 'custom',
    # 'mathtext.bf': 'Arial',
    # # 'mathtext.rm': 'Arial',
    'font.size': 14,
    'svg.fonttype': 'none',
    # 'xtick.minor.size': 4,
    # 'ytick.minor.size': 4,
    # 'xtick.minor.width': 1.5,
    # 'ytick.minor.width': 1.5,
    # 'xtick.major.size': 6,
    # 'ytick.major.size': 6,
    # 'xtick.major.width': 2,
    # 'ytick.major.width': 2,
    # 'axes.linewidth': 2,
    # 'axes.labelpad': 14,
    # 'axes.spines.top': False,
    # 'axes.spines.right': False,
    # 'legend.frameon': False,
    # 'hatch.linewidth': 1
})


def get_ratios_for_protein(df, uniprot, condition):
    protein = df.loc[uniprot]

    dmso = protein.filter(regex=condition)
    ratios = dmso.filter(regex='list')
    ratios.index = pd.Index([3, 6, 24, 72], name='time')
    ratios = ratios.explode().reset_index()
    
    return ratios

def get_fit_params(df, uniprot, condition):
    protein = df.loc[uniprot]
    params = ['a', 'b', 'k']
    cols = [f'{condition}_{x}' for x in params]

    return dict(zip(params, protein[cols].values))

def get_deg_params(uniprot, condition):
    protein = light_grouped.loc[uniprot]
    params = ['a', 'b', 'k']
    cols = [f'{condition}_{x}' for x in params]
    deg_params = ['a_deg', 'b_deg', 'k_deg']

    return dict(zip(deg_params, protein[cols].values))

def get_syn_params(uniprot, condition):
    protein = heavy_grouped.loc[uniprot]
    params = ['a', 'b', 'k']
    cols = [f'{condition}_{x}' for x in params]
    syn_params = ['a_syn', 'b_syn', 'k_syn']
    
    return dict(zip(syn_params, protein[cols].values))

def plot(uniprot):
    palette = sns.color_palette('Set2')

    light_protein = light_grouped.loc[uniprot]
    heavy_protein = heavy_grouped.loc[uniprot]

    dmso_light_ratios = get_ratios_for_protein(light_grouped, uniprot, 'dmso')
    dmso_heavy_ratios = get_ratios_for_protein(heavy_grouped, uniprot, 'dmso')
    abd_light_ratios = get_ratios_for_protein(light_grouped, uniprot, '957')
    abd_heavy_ratios = get_ratios_for_protein(heavy_grouped, uniprot, '957')

    dmso_light_params = get_deg_params(uniprot, 'dmso')
    abd_light_params = get_deg_params(uniprot, 'abd')
    dmso_heavy_params = get_syn_params(uniprot, 'dmso')
    abd_heavy_params = get_syn_params(uniprot, 'abd')

    x2 = np.linspace(1,72,250)

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharey='all', sharex='all')
    plt.subplots_adjust(wspace=0.1, hspace=0.1)

    ax1.plot(x2, degradation(x2, **dmso_light_params), label='fit', color=palette[1])
    ax1.scatter(dmso_light_ratios['time'].values, dmso_light_ratios[uniprot], color=palette[-1], alpha=0.5)
    t = f"A={dmso_light_params.get('a_deg'):.2f} B={dmso_light_params.get('b_deg'):.2f} K={dmso_light_params.get('k_deg'):.2f} R2={light_protein.get('dmso_rsq_deg', 0.0):.2f}\nT1/2={light_protein.get('dmso_deg_halflife'):.2f}" 
    ax1.annotate(t, (0.025, 0.975), xycoords='axes fraction', size='x-small', va='top')

    ax2.plot(x2, degradation(x2, **abd_light_params), label='fit', color=palette[1])
    ax2.scatter(abd_light_ratios['time'].values, abd_light_ratios[uniprot], color=palette[2], alpha=0.5)
    t = f"A={abd_light_params.get('a_deg'):.2f} B={abd_light_params.get('b_deg'):.2f} K={abd_light_params.get('k_deg'):.2f} R2={light_protein.get('abd_rsq_deg', 0.0):.2f}\nT1/2={light_protein.get('abd_deg_halflife'):.2f}" 
    ax2.annotate(t, (0.025, 0.975), xycoords='axes fraction', size='x-small', va='top')

    ax3.plot(x2, synthesis(x2, **dmso_heavy_params), label='fit', color=palette[0])
    ax3.scatter(dmso_heavy_ratios['time'].values, dmso_heavy_ratios[uniprot], color=palette[-1], alpha=0.5)
    t = f"A={dmso_heavy_params.get('a_syn'):.2f} B={dmso_heavy_params.get('b_syn'):.2f} K={dmso_heavy_params.get('k_syn'):.2f} R2={heavy_protein.get('dmso_rsq_syn', 0.0):.2f}\nT1/2={heavy_protein.get('dmso_syn_halflife'):.2f}" 
    ax3.annotate(t, (0.025, 0.975), xycoords='axes fraction', size='x-small', va='top')


    ax4.plot(x2, synthesis(x2, **abd_heavy_params), label='fit', color=palette[0])
    ax4.scatter(abd_heavy_ratios['time'].values, abd_heavy_ratios[uniprot], color=palette[2], alpha=0.5)
    t = f"A={abd_heavy_params.get('a_syn'):.2f} B={abd_heavy_params.get('b_syn'):.2f} K={abd_heavy_params.get('k_syn'):.2f} R2={heavy_protein.get('abd_rsq_syn', 0.0):.2f}\nT1/2={heavy_protein.get('abd_syn_halflife'):.2f}" 
    ax4.annotate(t, (0.025, 0.975), xycoords='axes fraction', size='x-small', va='top')


    ax1.set_ylabel('Degradation')
    ax3.set_ylabel('Synthesis')
    ax3.set_xlabel('DMSO')
    ax4.set_xlabel('ABD957')

    
    fig.suptitle(light_protein[('symbol', 'first')], fontweight='bold')

#%%
from matplotlib.path import Path
from matplotlib.axes import Axes
from matplotlib.axes._axes import _make_inset_locator
from matplotlib.transforms import Bbox, Transform, IdentityTransform, Affine2D
from matplotlib.backend_bases import RendererBase
import matplotlib._image as _image
import numpy as np


#%%
def plot_overlay(uniprot):
    palette = sns.color_palette('Set2')

    light_protein = light_grouped.loc[uniprot]
    heavy_protein = heavy_grouped.loc[uniprot]

    dmso_light_ratios = get_ratios_for_protein(light_grouped, uniprot, 'dmso')
    dmso_heavy_ratios = get_ratios_for_protein(heavy_grouped, uniprot, 'dmso')
    abd_light_ratios = get_ratios_for_protein(light_grouped, uniprot, '957')
    abd_heavy_ratios = get_ratios_for_protein(heavy_grouped, uniprot, '957')

    dmso_light_params = get_deg_params(uniprot, 'dmso')
    abd_light_params = get_deg_params(uniprot, 'abd')
    dmso_heavy_params = get_syn_params(uniprot, 'dmso')
    abd_heavy_params = get_syn_params(uniprot, 'abd')

    print(dmso_light_params)
    print(abd_light_params)
    print(dmso_heavy_params)
    print(abd_heavy_params)

    x2 = np.linspace(1,72,250)

    fig, (ax1, ax2) = plt.subplots(1, 2, sharey='all', sharex='all', figsize=(12,4))
    plt.subplots_adjust(wspace=0.1, top=0.8)

    ax1.plot(x2, degradation(x2, **dmso_light_params), color=BLUE_COLOR)
    ax1.scatter(dmso_light_ratios['time'].values, dmso_light_ratios[uniprot], label='DMSO', color=BLUE_COLOR, alpha=0.5)
    # t = f"A={dmso_light_params.get('a_deg'):.2f} B={dmso_light_params.get('b_deg'):.2f} K={dmso_light_params.get('k_deg'):.2f} R2={light_protein.get('dmso_rsq_deg', 0.0):.2f}\nT1/2={light_protein.get('dmso_deg_halflife'):.2f}" 
    # ax1.annotate(t, (0.025, 0.975), xycoords='axes fraction', size='x-small', va='top')

    ax1.plot(x2, degradation(x2, **abd_light_params), color=RED_COLOR)
    ax1.scatter(abd_light_ratios['time'].values, abd_light_ratios[uniprot], label='ABD957', color=RED_COLOR, alpha=0.5)
    # t = f"A={abd_light_params.get('a_deg'):.2f} B={abd_light_params.get('b_deg'):.2f} K={abd_light_params.get('k_deg'):.2f} R2={light_protein.get('abd_rsq_deg', 0.0):.2f}\nT1/2={light_protein.get('abd_deg_halflife'):.2f}" 
    # ax1.annotate(t, (0.025, 0.975), xycoords='axes fraction', size='x-small', va='top')

    ax2.plot(x2, synthesis(x2, **dmso_heavy_params), color=BLUE_COLOR)
    ax2.scatter(dmso_heavy_ratios['time'].values, dmso_heavy_ratios[uniprot], label='DMSO', color=BLUE_COLOR, alpha=0.5)
    # t = f"A={dmso_heavy_params.get('a_syn'):.2f} B={dmso_heavy_params.get('b_syn'):.2f} K={dmso_heavy_params.get('k_syn'):.2f} R2={heavy_protein.get('dmso_rsq_syn', 0.0):.2f}\nT1/2={heavy_protein.get('dmso_syn_halflife'):.2f}" 
    # ax2.annotate(t, (0.025, 0.975), xycoords='axes fraction', size='x-small', va='top')


    ax2.plot(x2, synthesis(x2, **abd_heavy_params), color=RED_COLOR)
    ax2.scatter(abd_heavy_ratios['time'].values, abd_heavy_ratios[uniprot], label='ABD957', color=RED_COLOR, alpha=0.5)
    # t = f"A={abd_heavy_params.get('a_syn'):.2f} B={abd_heavy_params.get('b_syn'):.2f} K={abd_heavy_params.get('k_syn'):.2f} R2={heavy_protein.get('abd_rsq_syn', 0.0):.2f}\nT1/2={heavy_protein.get('abd_syn_halflife'):.2f}" 
    # ax2.annotate(t, (0.025, 0.975), xycoords='axes fraction', size='x-small', va='top')


    # axins = ZoomViewAxes(ax1, Bbox.from_bounds(0.6, 0.6, 0.35, 0.35), ax1.transAxes)  # Use the new zoom axes...
    axins = ax1.inset_axes([0.6, 0.6, 0.35, 0.35])
    axins.plot(x2, degradation(x2, **dmso_light_params), color=BLUE_COLOR)
    axins.scatter(dmso_light_ratios['time'].values, dmso_light_ratios[uniprot], label='DMSO', color=BLUE_COLOR, alpha=0.5)
    # t = f"A={dmso_light_params.get('a_deg'):.2f} B={dmso_light_params.get('b_deg'):.2f} K={dmso_light_params.get('k_deg'):.2f} R2={light_protein.get('dmso_rsq_deg', 0.0):.2f}\nT1/2={light_protein.get('dmso_deg_halflife'):.2f}" 
    # axins.annotate(t, (0.025, 0.975), xycoords='axes fraction', size='x-small', va='top')

    axins.plot(x2, degradation(x2, **abd_light_params), color=RED_COLOR)
    axins.scatter(abd_light_ratios['time'].values, abd_light_ratios[uniprot], label='ABD957', color=RED_COLOR, alpha=0.5)


    axins.set_xlim(1, 9)
    axins.set_ylim(0.5, 1.2)
    axins.set_xticks([3, 6])
    axins.set_xlim(left=2, right=7)
    axins.set_xticklabels(['3 h', '6 h'])

    ax1.indicate_inset_zoom(axins)

    ax1.set_xticks([3, 6, 12, 24, 48, 72])
    ax1.set_ylabel('Fraction of t${_0}$')
    ax1.set_xlabel('Time (hours)')
    ax2.set_xlabel('Time (hours)')
    ax1.set_title('Degradation')
    ax2.set_title('Synthesis')
    ax1.set_xlim(left=0)
    ax1.set_ylim(bottom=0)

    ax2.legend()

    # ax1.set_ylabel('Degradation')
    # ax3.set_ylabel('Synthesis')
    # ax3.set_xlabel('DMSO')
    # ax4.set_xlabel('ABD957')

    symbol = light_protein[('symbol', 'first')]
    fig.suptitle(symbol, fontweight='bold')
    fig.savefig(f'{symbol}_turnover.svg', transparent=True)
    fig.savefig(f'{symbol}_turnover.png')

plot_overlay('P01111')
# plot_overlay('Q9NZW5') # MPP6
plot_overlay('O75608') # LYPLA1
plot_overlay('Q14160') # SCRIB


#%%
# for uniprot in maybe_fast:
#     try:
#         plot_overlay(uniprot)
#     except:
#         continue


# %%


#%%

def plot_overlap(uniprot):
    palette = sns.color_palette('Set2')
    light_protein = light_grouped.loc[uniprot]
    heavy_protein = heavy_grouped.loc[uniprot]

    dmso_light_ratios = get_ratios_for_protein(light_grouped, uniprot, 'dmso')
    dmso_heavy_ratios = get_ratios_for_protein(heavy_grouped, uniprot, 'dmso')
    abd_light_ratios = get_ratios_for_protein(light_grouped, uniprot, '957')
    abd_heavy_ratios = get_ratios_for_protein(heavy_grouped, uniprot, '957')

    light_ratios = dmso_light_ratios.assign(Treatment='DMSO').append(abd_light_ratios.assign(Treatment='ABD957')).reset_index()
    heavy_ratios = dmso_heavy_ratios.assign(Treatment='DMSO').append(abd_heavy_ratios.assign(Treatment='ABD957')).reset_index()

    dmso_light_params = get_deg_params(uniprot, 'dmso')
    abd_light_params = get_deg_params(uniprot, 'abd')
    dmso_heavy_params = get_syn_params(uniprot, 'dmso')
    abd_heavy_params = get_syn_params(uniprot, 'abd')

    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, sharey=True, sharex=True, figsize=(10,4))
    plt.subplots_adjust(wspace=0.1, top=0.8)

    # ax1.scatter(dmso_light_ratios['time'].values, dmso_light_ratios[uniprot], color=palette[-1], alpha=0.5)
    # ax1.scatter(abd_light_ratios['time'].values, abd_light_ratios[uniprot], color=palette[2], alpha=0.5)
    # ax2.scatter(dmso_heavy_ratios['time'].values, dmso_heavy_ratios[uniprot], color=palette[-1], alpha=0.5)
    # ax2.scatter(abd_heavy_ratios['time'].values, abd_heavy_ratios[uniprot], color=palette[2], alpha=0.5)

    sns.boxplot(x='time', y=uniprot, hue='Treatment', data=light_ratios, ax=ax1, fliersize=0)
    sns.swarmplot(x='time', y=uniprot, hue='Treatment', data=light_ratios, dodge=True, edgecolor='.25', alpha=1, color='.25', ax=ax1)
    sns.boxplot(x='time', y=uniprot, hue='Treatment', data=heavy_ratios, ax=ax2, fliersize=0)
    sns.swarmplot(x='time', y=uniprot, hue='Treatment', data=heavy_ratios, dodge=True, edgecolor='.25', alpha=1, color='.25', ax=ax2)

    # ax1.set_xticks([0, 3, 6, 24, 72])
    ax1.get_legend().remove()
    ax2.set_ylabel('')
    ax1.set_title('Degradation')
    ax2.set_title('Synthesis')
    ax1.set_xlabel('Time (hours)')
    ax2.set_xlabel('Time (hours)')
    ax1.set_ylabel('Fraction of ${t_0}$')
    fig.suptitle(light_protein[('symbol', 'first')], fontweight='bold')

# plot_overlap('Q9NZW5') # MPP6
# plot_overlap('Q03113') # GNA12
# plot_overlap('Q5VST6') # ABHD17B
# plot_overlap('O75608') # LYPLA1
# plot_overlap('P01112') # HRAS
plot_overlap('P01111') # NRAS
plot('P01111') # NRAS
