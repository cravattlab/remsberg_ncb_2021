#%%
import sys
sys.path.insert(0, '/mnt/d/OneDrive - The Scripps Research Institute/ABHD17/NCB_Reviews/figures')


from matplotlib.patches import Patch
import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
import numpy as np
import seaborn as sns

import pathlib


#%%

prm_folder = pathlib.Path('prm')
prm_files = list(prm_folder.glob('*.csv'))

df = pd.DataFrame()

for file in prm_files:
    enzyme = file.stem.split()[-1]
    df = df.append(pd.read_csv(file).assign(enzyme=enzyme))

df = df.rename(columns={'Unnamed: 0': 'time'})

abd957 = df[['time', 'enzyme'] + list(df.filter(regex='957'))]
abd957['time'] = abd957.time.apply(lambda x: f'{x}_957')
abd957.columns = ['time', 'enzyme', 'rep1', 'rep2', 'rep3']
abd362 = df[['time', 'enzyme'] + list(df.filter(regex='362'))]
abd362['time'] = abd362.time.apply(lambda x: f'{x}_362')
abd362.columns = ['time', 'enzyme', 'rep1', 'rep2', 'rep3']
# df = abd957.append(abd362)
df = abd362.append(abd957)

main_df = df[df.enzyme.str.startswith('ABHD17')]

abhd17a = df[df.enzyme.str.startswith('ABHD17A')]

main_df = (main_df
    .set_index(['time', 'enzyme'])
    .stack()
    .reset_index()
    .drop(columns='level_2')
    .rename(columns={0: 'activity'})
)


plt.rcParams.update({
    'font.weight': 'normal',
    'axes.labelweight': 'normal',
    'axes.titleweight': 'normal',
    'font.sans-serif': 'Arial',
    'font.family': 'sans-serif',
    'mathtext.fontset': 'custom',
    # 'mathtext.bf': 'Arial',
    'mathtext.rm': 'Arial',
    'font.size': 12,
    'xtick.minor.size': 4,
    'ytick.minor.size': 4,
    'xtick.minor.width': 1.25,
    'ytick.minor.width': 1.25,
    'xtick.major.size': 6,
    'ytick.major.size': 6,
    'xtick.major.width': 1.25,
    'ytick.major.width': 1.25,
    'axes.linewidth': 1.25,
    'axes.labelpad': 8,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'legend.frameon': False,
    'hatch.linewidth': 1.25,
    'svg.fonttype': 'none',
})

#%%
def old_plot():
    fig, ax = plt.subplots()

    palette = sns.color_palette('Blues', n_colors=3)
    enzymes = list(main_df.enzyme.unique())

    common_args = dict(x='time', y='activity', hue='enzyme', data=main_df, ax=ax, hue_order=enzymes)
    sns.barplot(**common_args, palette=palette, errwidth=1, capsize=0.1, errcolor='0.2', ci='sd')
    sns.stripplot(**common_args, dodge=True, edgecolor='0.1', linewidth=1.2, jitter=False, alpha=0.8, size=4)
    # sns.swarmplot(**common_args, dodge=True, edgecolor='0.1', linewidth=1, alpha=0.8, size=4)
    for c in ax.get_children():
        if isinstance(c, matplotlib.collections.PathCollection):
            c.set_facecolor('none')

    ax.legend().remove()
    patches = [Patch(color=palette[i], label=enzyme) for i, enzyme in enumerate(enzymes)]
    ax.legend(handles=patches, loc='upper right')

    ax.set_ylabel('Activity (% of DMSO)')
    ax.set_xlabel(None)
    ax.set_xticklabels(['2h', '24h', '48h', '72h'] * 2)
    ax.set_ylim(top=100)
    ax.set_yticks(np.linspace(0, 100, 5))

    fig.suptitle(r'$\it{In\ situ}$ Time-course')

    plt.subplots_adjust(bottom=0.2)
    annotation_args = dict(
        xycoords='axes fraction',
        ha='center',
        arrowprops={
            'arrowstyle': '-[, widthB=5.8, lengthB=0.3',
        }
    )
    plt.annotate(r'$\bf{4}$ (1 µM)', (0.25, -0.14), xytext=(0.25, -0.24), **annotation_args)
    plt.annotate('ABD957 (1 µM)', (0.75, -0.14), xytext=(0.75, -0.24), **annotation_args)
    # fig

# %%

def k(x):
    split = x.str.split('h_').str
    return split.get(1) + split.get(0)

lol = df[(df.enzyme.str.startswith('ABHD17'))]
lol['means'] = lol.filter(regex='rep').mean(axis=1)
lol['stdev'] = lol.filter(regex='rep').std(axis=1, ddof=0)
lol = lol.sort_values('time', key=k)

a = lol[lol.enzyme == 'ABHD17A']
b = lol[lol.enzyme == 'ABHD17B']
c = lol[lol.enzyme == 'ABHD17C']

N = len(a)
ind = np.arange(N)
width = 0.27

palette = sns.color_palette('Blues', n_colors=3) 

fig, ax = plt.subplots(figsize=(4, 2.5))

bar_args = dict(zorder=1)
ax.bar(ind, a.means.values, width, label='ABHD17A', color=palette[0], **bar_args)
ax.bar(ind + width, b.means.values, width, label='ABHD17B', color=palette[1], **bar_args)
ax.bar(ind + width * 2, c.means.values, width, label='ABHD17C', color=palette[2], **bar_args)

scatter_args = dict(facecolor='white', linewidth=1.1, zorder=2, edgecolors='.1', alpha=0.9, s=12)
n = ind.repeat(3).astype(float)
# n[2::3] += width/4
# n[::3] += -width/4

ax.scatter(n, a.filter(regex='rep').values.flat, **scatter_args)

n = (ind + width).repeat(3).astype(float)
# n[2::3] += width/4
# n[::3] += -width/4

ax.scatter(n, b.filter(regex='rep').values.flat, **scatter_args)

n = (ind + width * 2).repeat(3).astype(float)
# n[2::3] += width/4
# n[::3] += -width/4
ax.scatter(n, c.filter(regex='rep').values.flat, **scatter_args)

# error_args = dict(elinewidth=1, lolims=True, ecolor='0.2', zorder=1, capsize=0, linestyle='none')
error_args = dict(elinewidth=1.1, lolims=True, ecolor='0.1', zorder=1, capsize=0, linestyle='none', alpha=0.9)
_, caplines_a, _ = ax.errorbar(ind, a.means.values, yerr=a.stdev.values, **error_args)
_, caplines_b, _ = ax.errorbar(ind + width, b.means.values, yerr=b.stdev.values, **error_args)
_, caplines_c, _ = ax.errorbar(ind + width * 2, c.means.values, yerr=c.stdev.values, **error_args)

for c in caplines_a + caplines_b + caplines_c:
    c.set_marker('_')
    c.set_markersize(6)

ax.set_ylabel('Activity (% of DMSO)')
ax.set_xlabel(None)
ax.set_xticks(ind + width)
ax.set_xticklabels(['2 ', '24h', '48h', '72h'] * 2)
ax.set_xlim(-width, N-width/1.5)
ax.set_ylim(top=100)
ax.set_yticks(np.linspace(0, 100, 5))
ax.legend(loc='upper right')

fig.suptitle(r'$\it{In\ situ}$ Time-course')

plt.subplots_adjust(bottom=0.2)
annotation_args = dict(
    xycoords='axes fraction',
    ha='center',
    arrowprops={
        'arrowstyle': '-[, widthB=5.8, lengthB=0.3',
    }
)
plt.annotate(r'$\bf{4}$ (1 µM)', (0.25, -0.14), xytext=(0.25, -0.24), **annotation_args)
plt.annotate('ABD957 (1 µM)', (0.75, -0.14), xytext=(0.75, -0.24), **annotation_args)

plt.tight_layout()

fig.savefig('output/prm.svg', transparent=True)

# %%
