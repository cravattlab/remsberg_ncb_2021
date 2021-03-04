import sys
sys.path.insert(0, '/mnt/d/OneDrive - The Scripps Research Institute/ABHD17/NCB_Reviews/figures/figure3/data')
import rifl.utils as utils

from matplotlib.pyplot import cm as color_map
from matplotlib.ticker import FuncFormatter, AutoMinorLocator
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import pandas as pd
import seaborn as sns

import pathlib
import statistics

import os
os.chdir('/mnt/d/OneDrive - The Scripps Research Institute/ABHD17/NCB_Reviews/figures/figure1')


DATA_INPUT_PATH = pathlib.Path('input').absolute()
DATA_OUTPUT_PATH = pathlib.Path('output').absolute()

blues_palette = sns.color_palette('Blues', n_colors=3) 


matplotlib.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': 'Arial',
    'font.cursive': 'Arial',
    # 'svg.fonttype': 'none'
    'mathtext.fontset': 'custom',
    # 'mathtext.bf': 'Arial',
    'mathtext.rm': 'Arial',
    'svg.fonttype': 'none',
})

FIGURE1_RED = '#ee2822ff'

symbol_highlight_colormap = {
    'ABHD17A': blues_palette[-1],
    'ABHD17B': blues_palette[-1],
    'ABHD17C': blues_palette[-1],
    'LYPLA1': FIGURE1_RED,
    'LYPLA2': FIGURE1_RED,
    'ABHD10': FIGURE1_RED,
}

def heatmap(df, output_path, title=None):
    fig, ax = plt.subplots()
    fig.set_size_inches(14, 3, forward=True)

    cmap = color_map.get_cmap('Blues_r')
    cmap.set_over('white')
    cmap.set_under('yellow') # there shouldn't be any under..
    cmap.set_bad('grey')

    ax.tick_params(axis='both', which='both', length=0)
    ax.xaxis.set_ticks_position('top')

    sns.heatmap(df, cmap=cmap, linewidths=0.3, linecolor='k', square=True, ax=ax, vmin=0.05, vmax=1, cbar_kws={
        'orientation': 'horizontal',
        'label': 'Activity (% of DMSO)',
        'shrink': 0.15,
        'pad': 0.015,
        'anchor': (1, 1),
        # 'panchor': False,
        'aspect': 12,
        'use_gridspec': False,
        'format': FuncFormatter(lambda x, _: f'{x*100:0.0f}'),
        'ticks': [0.05, 0.25, 0.5, 0.75, 1],
    })

    x_texts = ax.set_xticklabels(df.columns, size=11.24, rotation='vertical')

    for label in x_texts:
        text = label.get_text()
        if symbol_highlight_colormap.get(text):
            label.set_color(symbol_highlight_colormap[text])
            label.set_fontweight('bold')

    cbar = ax.collections[0].colorbar
    cbar.ax.set_frame_on(True)
    cbar.ax.tick_params(labelsize=11.24)

    # rotate ticks back to normal
    y_texts = ax.set_yticklabels(df.index, size=11.24, rotation=0)
    ax.set_xlabel('')

    for label in y_texts:
        if label.get_text().isdigit():
            label.set_fontweight('bold')

    # box around the heatmap
    for _, spine in ax.spines.items():
        spine.set_visible(True)

    # tight layout ruins everything
    # plt.tight_layout()

    if title:
        fig.suptitle(title)

    plt.savefig(DATA_OUTPUT_PATH.joinpath(output_path), format='svg', transparent=True)
    plt.close()

def bar_plot(df, output_path, figsize=(10, 2.5)):
    fig, ax = plt.subplots(figsize=figsize)

    N = len(df)
    ind = np.arange(N)
    width = 0.7

    rects = ax.bar(ind, df.ratio.values, width, color=blues_palette[-1])

    replicate_values = df.reset_index().ratio_list.explode()

    scatter_args = dict(facecolor='white', zorder=2, edgecolors='.1', alpha=0.9, s=12)
    ax.scatter(replicate_values.index.values, replicate_values.values, **scatter_args)

    stdevs = df.ratio_list.apply(statistics.pstdev)

    error_args = dict(elinewidth=1, lolims=True, ecolor='0.2', zorder=1, capsize=0, linestyle='none')
    _, caplines_a, _ = ax.errorbar(ind, df.ratio.values, yerr=stdevs.values, **error_args)

    for c in caplines_a:
        c.set_marker('_')
        c.set_markersize(6)


    ax.axhline(y=0.5, zorder=-1, linestyle='--', linewidth=0.75, color='lightgrey')
    ax.axhline(y=1, zorder=-1, linestyle='--', linewidth=0.75, color='lightgrey')

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # texts[10].set_color('red')

    ax.set_ylabel('Activity (% of DMSO)', fontsize=12)
    ax.set_xticks(ind)
    # ax.set_yticks(np.linspace(0, 1.50, 7))
    ax.tick_params(direction='out')
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    # ax.yaxis.set_minor_locator(AutoMinorLocator(2))

    ax.yaxis.set_major_formatter(FuncFormatter(lambda x, _: f'{x*100:0.0f}'))

    plt.xlim(left=-width, right=N-width/2)
    plt.yticks(np.arange(0, float(replicate_values.max()), step=0.25))

    texts = ax.set_xticklabels(df.index.values, rotation='vertical', family='sans-serif')

    for label in texts:
        text = label.get_text()
        if symbol_highlight_colormap.get(text):
            label.set_color(symbol_highlight_colormap[text])
            label.set_fontweight('bold')

    plt.tight_layout()
    plt.savefig(DATA_OUTPUT_PATH.joinpath(output_path), format='svg')
    plt.close()


def prep_data_og(df, subset=None):
    df = df.drop(columns='num_peptides', level=1)
    df.columns = df.columns.get_level_values(0)
    df = df.drop(columns='uniprot')
    df = df.set_index('symbol')
    df = df.apply(pd.to_numeric, errors='coerce')
    df = df.dropna(how='all')
    if subset:
        df = df.filter(regex=subset)
    return df

def prep_data(df, subset=None, include_ratio_list=False):
    ratio_list = df.loc[:, pd.IndexSlice[:, 'ratio_list']]
    ratio_list = ratio_list.set_index(df.loc[:, 'symbol'].values.flatten())

    def ratio_list_to_float(ratios):
        if isinstance(ratios, str):
            return list(map(float, ratios.split(', ')))
        return '-'

    ratio_list = ratio_list.replace('-', np.nan).transform(lambda x: x.apply(ratio_list_to_float))
    df = df.drop(columns=['num_peptides', 'ratio_list'], level=1)

    df.columns = df.columns.get_level_values(0)
    df = df.drop(columns='uniprot')
    df = df.set_index('symbol')
    df = df.apply(pd.to_numeric, errors='coerce')
    df = df.dropna(how='all')
    if include_ratio_list:
        ratio_list.columns = ratio_list.columns.get_level_values(0)
        df = df.join(ratio_list, rsuffix=' (ratio list)')
    if subset:
        df = df.filter(regex=subset)
    return df

def prep_invitro(df, concentration=None):
    df = prep_data(df, concentration)
    df = df[sorted(df.columns)]
    df = df.T
    df = df[sorted(df.columns)]
    df = df.rename(index={'5': 'ABD957'})
    return df

def prep_all_invitro(df):
    # df = prep_data(df, subset='[1-5]')
    df = prep_data(df)
    df = df[sorted(df.columns)]
    df = df.T
    df = df[sorted(df.columns)]
    df.index = df.index.str.replace('Compound\s(\d).+', lambda m: m.group(1), regex=True)
    df = df.rename(index={'5': 'ABD957'})
    return df

def prep_insitu(df, inhibitor):
    df = prep_data(df, inhibitor, include_ratio_list=True)
    df.columns = ['ratio', 'ratio_list']
    df = df.sort_values(by=['ratio', 'symbol'], ascending=[True, False])
    return df

def prep_all_insitu(df):
    df = prep_data(df)
    df = df[sorted(df.columns)]
    df = df.T
    df = df[sorted(df.columns)]
    df = df.rename(index={'5': 'ABD957'})
    return df

def main():
    latest_insitu_path = utils.get_newest_file(DATA_OUTPUT_PATH, 'in_situ_*.csv')
    insitu = pd.read_csv(latest_insitu_path, header=[0,1])
    # insitu_all = prep_all_insitu(insitu)
    # heatmap(1/insitu_all, 'insitu_all.svg')

    latest_invitro_path = utils.get_newest_file(DATA_OUTPUT_PATH, 'in_vitro_*.csv')
    invitro = pd.read_csv(latest_invitro_path, header=[0,1])
    # invitro_1um = invitro.drop(columns=invitro.filter(regex='10 µM'))
    # invitro_10um = invitro.drop(columns=invitro.filter(regex='1 µM'))

    invitro_1um = prep_all_invitro(invitro.drop(columns=invitro.filter(regex='10 µM')))
    invitro_10um = prep_all_invitro(invitro.drop(columns=invitro.filter(regex='1 µM')))
    # # invitro_all = prep_all_invitro(invitro)
    # # invitro_all.loc['Compound 4 (1 µM)', 'FAM135B'] = 0.79
    invitro_1um.loc['4', 'FAM135B'] = 0.79
    # invitro_1um = invitro_1um.drop(columns='CTSA')

    # plt.rcParams.update({
    #     'font.weight': 'normal',
    #     'axes.labelweight': 'normal',
    #     'axes.titleweight': 'normal',
    #     'font.sans-serif': 'Arial',
    #     'font.family': 'sans-serif',
    #     'font.cursive': 'Arial',
    #     'mathtext.fontset': 'custom',
    #     # 'mathtext.bf': 'Arial',
    #     'mathtext.rm': 'Arial',
    #     # 'font.size': 14,
    #     'xtick.minor.size': 2,
    #     'ytick.minor.size': 2,
    #     'xtick.minor.width': 1.25,
    #     'ytick.minor.width': 1.25,
    #     'xtick.major.size': 4,
    #     'ytick.major.size': 4,
    #     'xtick.major.width': 1.25,
    #     'ytick.major.width': 1.25,
    #     'axes.linewidth': 1.25,
    #     # 'axes.labelpad': 8,
    #     'axes.spines.top': False,
    #     'axes.spines.right': False,
    #     'legend.frameon': False,
    #     'hatch.linewidth': 1,
    #     'svg.fonttype': 'none',
    # })

    heatmap(invitro_1um, 'invitro_1um.svg', r'OCI-AML3 ($\it{in\ vitro}$, 1 µM)')
    heatmap(invitro_10um, 'invitro_10um.svg')

    # latest_crispr_path = utils.get_newest_file(DATA_OUTPUT_PATH, 'crisprko_*.csv')
    # crispr = pd.read_csv(latest_crispr_path, header=[0,1])
    # crispr_all = prep_all_insitu(crispr)
    # heatmap(1/crispr_all, 'crispr_all.svg')

    plt.rcParams.update({
        'font.weight': 'normal',
        'axes.labelweight': 'normal',
        'axes.titleweight': 'normal',
        'font.sans-serif': 'Arial',
        'font.family': 'sans-serif',
        'font.cursive': 'Arial',
        'mathtext.fontset': 'custom',
        # 'mathtext.bf': 'Arial',
        'mathtext.rm': 'Arial',
        # 'font.size': 14,
        'xtick.minor.size': 2,
        'ytick.minor.size': 2,
        'xtick.minor.width': 1.25,
        'ytick.minor.width': 1.25,
        'xtick.major.size': 4,
        'ytick.major.size': 4,
        'xtick.major.width': 1.25,
        'ytick.major.width': 1.25,
        'axes.linewidth': 1.25,
        # 'axes.labelpad': 8,
        'axes.spines.top': False,
        'axes.spines.right': False,
        'legend.frameon': False,
        'hatch.linewidth': 1,
        'svg.fonttype': 'none',
    })

    abd957 = prep_insitu(insitu, '957').dropna()
    palm_m = prep_insitu(insitu, 'Palm M').dropna()
    hdfp = prep_insitu(insitu, 'HDFP').dropna()
    jjh254 = prep_insitu(insitu, 'JJH254').dropna()
    abhd13i = prep_insitu(insitu, 'ABD298').dropna()

    # hdfp.loc['SEC11A', 'ratio_list'] = [1/0.66, 0.8, 0.8620689655172414]

    bar_plot(abd957, 'abd957.svg')
    bar_plot(palm_m, 'palm_m.svg')
    bar_plot(hdfp, 'hdfp.svg')
    bar_plot(jjh254, 'jjh254.svg')
    bar_plot(abhd13i, 'abhd13i.svg')


if __name__ == '__main__':
    main()

