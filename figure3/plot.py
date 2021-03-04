from constants import (
    CACHE_DIR,
    DYNAMIC_FOLD_CHANGE_THRESHOLD,
    CHASE_FOLD_CHANGE_THRESHOLD,
    DATA_OUTPUT_PATH
)
from filter import filter_time_course_data
from utils import get_path_for_fig, get_newest_file
from db import UniprotFlatfile
import constants as k

from matplotlib.ticker import AutoMinorLocator
from adjustText import adjust_text
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

import joblib
import math


mem = joblib.Memory(CACHE_DIR)
read_csv = mem.cache(pd.read_csv)

plt.rcParams.update({
    'font.weight': 'normal',
    'axes.labelweight': 'normal',
    'axes.titleweight': 'normal',
    'font.sans-serif': 'Arial',
    'font.family': 'sans-serif',
    'mathtext.fontset': 'custom',
    'mathtext.default': 'regular',
    'mathtext.bf': 'Arial',
    'mathtext.it': 'Arial',
    'mathtext.cal': 'Arial',
    'mathtext.rm': 'Arial',
    'svg.fonttype': 'none',
    'font.size': 12,
    'xtick.major.size': 6,
    'ytick.major.size': 6,
    'ytick.major.width': 1.25,
    'xtick.major.width': 1.25,
    'axes.linewidth': 1.25,
    'axes.labelpad': 12,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'legend.frameon': False,
    'hatch.linewidth': 1
})

get_filtered_time_course_data = mem.cache(filter_time_course_data)

def process_timecourse_data(file_pattern):
    time_course_path = get_newest_file(file_pattern, DATA_OUTPUT_PATH)
    df = pd.read_excel(time_course_path, index_col='symbol')

    df = pd.concat([
        df['DMSO t0 (mean)']/df['DMSO t1 (mean)'],
        df['ABD957 t1 (mean)']/df['DMSO t1 (mean)'],
        df['Palm M t1 (mean)']/df['DMSO t1 (mean)'],
    ], axis=1)

    df.columns=[
        'DMSO t0/t1',
        'ABD957 t1/DMSO t1',
        'Palm M t1/DMSO t1'
    ]

    df = df.groupby(level=0, axis=1).mean()

    df = np.log2(df)

    return df.reset_index()



def repel_labels(ax, x, y, labels, k=0.8, fontsize=10):
    G = nx.DiGraph()
    data_nodes = []
    init_pos = {}
    for xi, yi, label in zip(x, y, labels):
        data_str = 'data_{0}'.format(label)
        G.add_node(data_str)
        G.add_node(label)
        G.add_edge(label, data_str)
        data_nodes.append(data_str)
        init_pos[data_str] = (xi, yi)
        init_pos[label] = (xi, yi)

    pos = nx.spring_layout(G, pos=init_pos, fixed=data_nodes, k=k)

    # undo spring_layout's rescaling
    pos_after = np.vstack([pos[d] for d in data_nodes])
    pos_before = np.vstack([init_pos[d] for d in data_nodes])
    scale, shift_x = np.polyfit(pos_after[:,0], pos_before[:,0], 1)
    scale, shift_y = np.polyfit(pos_after[:,1], pos_before[:,1], 1)
    shift = np.array([shift_x, shift_y])
    for key, val in pos.items():
        pos[key] = (val*scale) + shift

    for label, data_str in G.edges():
        ax.annotate(label,
                    xy=pos[data_str], xycoords='data',
                    xytext=pos[label], textcoords='data',
                    ha='right',
                    va='center_baseline',
                    fontsize=fontsize,
                    arrowprops=dict(arrowstyle="->",
                                    shrinkA=0, shrinkB=2,
                                    connectionstyle="arc3", 
                                    ), )
    # expand limits
    all_pos = np.vstack(pos.values())
    x_span, y_span = np.ptp(all_pos, axis=0)
    mins = np.min(all_pos-x_span*0.15, 0)
    maxs = np.max(all_pos+y_span*0.15, 0)
    ax.set_xlim([mins[0], maxs[0]])
    ax.set_ylim([mins[1], maxs[1]])

def scatter(file_pattern, title, maxval=None):
    df = process_timecourse_data(file_pattern)

    MAX = df.max(numeric_only=True).max() + 0.1
    MIN = df.min(numeric_only=True).min() - 0.1
    # select largest absolute value and round up to next decimal

    MAX = math.ceil(max(MAX, abs(MIN)) * 10)/10

    if maxval is not None:
        MAX = maxval

    MIN = -MAX

    MIN_X = -1
    MIN_Y = -3

    _tick_range = [x for x in np.arange(math.ceil(MIN), MAX, 2) if x != 0]
    xtick_range = [x for x in np.arange(math.ceil(MIN_X), MAX, 2) if x != 0]
    ytick_range = [x for x in np.arange(math.ceil(MIN_Y), MAX, 2) if x != 0]

    dynamic_threshold = np.log2(DYNAMIC_FOLD_CHANGE_THRESHOLD)
    chase_threshold = np.log2(CHASE_FOLD_CHANGE_THRESHOLD)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))
    # fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 10))

    def plot_scatter(ax, data, ylabel_text, fontsize=10):
        data.columns = ['symbol', 'x', 'y']
        data['color'] = 'black'

        data.loc[data['x'] >= dynamic_threshold, 'color'] = 'green'
        data.loc[data['y'] >= chase_threshold, 'color'] = 'blue'
        data.loc[(data['x'] >= dynamic_threshold) & (data['y'] >= chase_threshold), 'color'] = 'red'

        data = data[data.symbol != 'HLA']
        
        ax.set(
            xlabel=r'$log_2(DMSO\ t_0/t_1)$',
            ylabel=ylabel_text,
            xticks=_tick_range,
            yticks=_tick_range,
            # xticks=xtick_range,
            # yticks=ytick_range,
        )

        ax.spines['left'].set_position('center')
        ax.spines['bottom'].set_position('center')

        ax.xaxis.set_label_coords(0.5, -0.03)
        ax.yaxis.set_label_coords(-0.03, 0.5)

        ax.xaxis.set_minor_locator(AutoMinorLocator(2))
        ax.yaxis.set_minor_locator(AutoMinorLocator(2))

        ax.tick_params(direction='inout', length=8)
        ax.tick_params(which='minor', direction='inout', length=6, width=1.25)

        ax.set_xlim(left=MIN, right=MAX)
        ax.set_ylim(bottom=MIN, top=MAX)
        ax.spines['left'].set_position(('data', 0))
        ax.spines['bottom'].set_position(('data', 0))

        ax.scatter(
            data['x'].values,
            data['y'].values,
            s=20, color=data['color'].values, alpha=0.5, zorder=1
        )

        to_label = data[(data.color == 'red')]
        # to_label = data[(data.color.isin(['red']))]
        # to_label = data[(data.symbol.isin(['HLA']))]

        for uniprot, row in to_label.iterrows():
            ax.annotate(
                row.symbol,
                xy=(row.x, row.y),
                # xytext=(5, 0),
                # textcoords='offset points',
                size=fontsize,
                ha='left',
            )

        ax.axvline(dynamic_threshold, color='red', linewidth=1, linestyle=':', zorder=-1)
        ax.axhline(chase_threshold, color='red', linewidth=1, linestyle=':', zorder=-1)
        ax.set_aspect('equal')

    palm_m = df[['symbol', 'DMSO t0/t1', 'Palm M t1/DMSO t1']]
    abd_957 = df[['symbol', 'DMSO t0/t1', 'ABD957 t1/DMSO t1']]

    plot_scatter(ax1, palm_m, r'$log_2(Palm\ M\ t_1/DMSO\ t_1)$', fontsize=6)
    plot_scatter(ax2, abd_957, r'$log_2(ABD957\ t_1/DMSO\ t_1)$')
    # plot_scatter(ax1, palm_m, r'$log_2(Palm\ M\ t_1/DMSO\ t_1)$', fontsize=4)
    # plot_scatter(ax2, abd_957, r'$log_2(ABD957\ t_1/DMSO\ t_1)$', fontsize=4)

    # plt.tight_layout()
    fig_output_path = get_path_for_fig(f'scatter_plots_{title}', datestamp='draft')
    fig.savefig(fig_output_path, format='svg', transparent=True)
    df.to_excel(fig_output_path.with_suffix('.xlsx'))


def scatter_10(time_course_path, title):
    df = pd.read_excel(time_course_path, index_col='symbol')

    df = pd.concat([
        df['DMSO t0 (mean)']/df['DMSO t1 (mean)'],
        df['ABD957 t1 (mean)']/df['DMSO t1 (mean)'],
    ], axis=1)

    df.columns=[
        'DMSO t0/t1',
        'ABD957 t1/DMSO t1',
    ]

    df = df.groupby(level=0, axis=1).mean()

    df = np.log2(df)

    df = df.reset_index()


    MAX = df.max(numeric_only=True).max() + 0.1
    MIN = df.min(numeric_only=True).min() - 0.1
    # select largest absolute value and round up to next decimal
    MAX = math.ceil(max(MAX, abs(MIN)) * 10)/10
    MIN = -MAX

    _tick_range = [x for x in np.arange(math.ceil(MIN), MAX, 2) if x != 0]

    dynamic_threshold = np.log2(DYNAMIC_FOLD_CHANGE_THRESHOLD)
    chase_threshold = np.log2(CHASE_FOLD_CHANGE_THRESHOLD)

    fig, ax = plt.subplots(figsize=(5, 5))

    def plot_scatter(ax, data, ylabel_text, fontsize=10):
        data = data.copy()
        data.columns = ['symbol', 'x', 'y']
        data['color'] = 'black'

        data.loc[data['x'] >= dynamic_threshold, 'color'] = 'green'
        data.loc[data['y'] >= chase_threshold, 'color'] = 'blue'
        data.loc[(data['x'] >= dynamic_threshold) & (data['y'] >= chase_threshold), 'color'] = 'red'
        
        ax.set(
            xlabel=r'$log_2(DMSO\ t_0/t_1)$',
            ylabel=ylabel_text,
            xticks=_tick_range,
            yticks=_tick_range,
        )

        ax.spines['left'].set_position('center')
        ax.spines['bottom'].set_position('center')

        ax.xaxis.set_label_coords(0.5, -0.03)
        ax.yaxis.set_label_coords(-0.03, 0.5)

        ax.xaxis.set_minor_locator(AutoMinorLocator(2))
        ax.yaxis.set_minor_locator(AutoMinorLocator(2))

        ax.tick_params(direction='inout', length=8)
        ax.tick_params(which='minor', direction='inout', length=6, width=1.25)

        ax.set_xlim(left=MIN, right=MAX)
        ax.set_ylim(bottom=MIN, top=MAX)

        ax.scatter(
            data['x'].values,
            data['y'].values,
            s=20, color=data['color'].values, alpha=0.5, zorder=1
        )

        to_label = data[
            (data.color.isin(['red', 'blue'])) |
            (data.symbol.isin(['NRAS', 'GNA12', 'MPP6', 'TAPBP']))
        ]


        for uniprot, row in to_label.iterrows():
            ax.annotate(
                row.symbol,
                xy=(row.x, row.y),
                # xytext=(5, 0),
                # textcoords='offset points',
                size=fontsize,
                ha='left',
            )

        ax.axvline(1, color='red', linewidth=1, linestyle=':', zorder=-1)
        ax.axhline(1, color='red', linewidth=1, linestyle=':', zorder=-1)
        ax.set_aspect('equal')

    abd_957 = df[['symbol', 'DMSO t0/t1', 'ABD957 t1/DMSO t1']]

    plot_scatter(ax, abd_957, r'$log_2(ABD957\ t_1/DMSO\ t_1)$')

    # plt.tight_layout()
    fig_output_path = get_path_for_fig(f'scatter_plots_{title}', datestamp='draft')
    fig.savefig(fig_output_path, format='svg', transparent=True)
    df.to_excel(fig_output_path.with_suffix('.xlsx'))


def main():
    scatter('filtered_timecourse_oci_p*.xlsx', 'oci_p')
    scatter('filtered_timecourse_oci_o*.xlsx', 'oci_o')
    scatter('nb4_6plex_10plex_*.xlsx', 'nb4')
    # scatter('filtered_timecourse_oci_on*.xlsx', 'oci_on')
    # scatter('filtered_timecourse_oci_combinedha_on*.xlsx', 'oci_combinedha_on')
    # scatter('filtered_timecourse_6plex_nb4*.xlsx')
    # scatter_10(get_newest_file('nb4_6plex_10plex.xlsx', DATA_OUTPUT_PATH), 'nb4_6plex_10_plex')


if __name__ == '__main__':
    main()
