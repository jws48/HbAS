import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from tslearn.utils import to_time_series_dataset
from tslearn.preprocessing import TimeSeriesScalerMeanVariance
from tslearn import metrics
from tslearn.metrics import dtw
from scipy.spatial.distance import cdist
import statsmodels.api as sa
import statsmodels.formula.api as sfa
import scikit_posthocs as sp
from numpy.random import seed
from numpy.random import randn
from numpy import mean
from numpy import std
from matplotlib import pyplot
from scipy import stats
from scipy.stats import ks_2samp
from scipy.stats import shapiro
from scipy.stats import normaltest
from scipy.stats import mannwhitneyu
    
def calc_dtw(ref, exp, timepoints):
    column_names = ["ORF", "dtw_score", "path"]
    df = pd.DataFrame(columns = column_names)
    res = []
    for index, row in ref.iterrows():
        x = ref.loc[[index]]
        x.columns = timepoints
        y = exp.loc[[index]]
        y.columns = timepoints
        dataset = x.append(y, ignore_index=True)
        dataset = dataset.to_numpy()
        dataset = to_time_series_dataset(dataset)
        scaler = TimeSeriesScalerMeanVariance(mu=0., std=1.)  # Rescale time series
        dataset_scaled = scaler.fit_transform(dataset)
        sz = dataset_scaled.shape[1]
        path, sim = metrics.dtw_path(dataset_scaled[0], dataset_scaled[1])
        res.append((index, sim, path))
    df = pd.DataFrame(data = res, columns=column_names)
    return df

## function to create dynamic time warping heatmap and optimal path trace
def pathmap(ref, exp, orf, timepoints):
    sns.set_style("white", {'xtick.bottom': True})
    sns.set_context("talk")
    x = ref.loc[[orf]]
    x.columns = timepoints
    y = exp.loc[[orf]]
    y.columns = timepoints
    dataset = x.append(y, ignore_index=True)
    dataset = dataset.to_numpy()
    dataset = to_time_series_dataset(dataset)
    scaler = TimeSeriesScalerMeanVariance(mu=0., std=1.)  # Rescale time series
    dataset_scaled = scaler.fit_transform(dataset)
    sz = dataset_scaled.shape[1]
    path, sim = metrics.dtw_path(dataset_scaled[0], dataset_scaled[1])

    plt.figure(1, figsize=(8, 8))

    # definitions for the axes
    left, bottom = 0.01, 0.1
    w_ts = h_ts = 0.2
    left_h = left + w_ts + 0.02
    left_cb = left + w_ts + 0.2
    width = height = 0.65
    bottom_h = bottom + height + 0.02

    rect_s_y = [left, bottom, w_ts, height]
    rect_gram = [left_h, bottom, width, height]
    rect_s_x = [left_h, bottom_h, width, h_ts]
    rect_cb = [left_cb, bottom+.05, width, .5]

    ax = plt.axes(rect_cb)
    ax_gram = plt.axes(rect_gram)
    ax_s_x = plt.axes(rect_s_x)
    ax_s_y = plt.axes(rect_s_y)

    mat = cdist(dataset_scaled[0], dataset_scaled[1])

    labels = list(range(3,51,3))

    ax_gram.set_yticks(np.arange(mat.shape[0]), minor=False)
    ax_gram.set_xticks(np.arange(mat.shape[1]), minor=False)

    ax_gram.invert_yaxis()
    ax_gram.set_yticklabels(labels, minor=False)
    ax_gram.set_xticklabels(labels, minor=False)
    ax_gram.yaxis.tick_right()

    heatmap = ax_gram.imshow(mat, cmap=plt.cm.Purples_r)
    # ax_gram.axis("off")
    ax_gram.autoscale(False)
    ax_gram.plot([j for (i, j) in path], [i for (i, j) in path], "w-",
             linewidth=3.)

    ax.axis("off")
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size = "5%", pad= 1)
#     cax = divider.append_axes("right", size = "5%", pad= .8)

    ax_s_x.plot(np.arange(sz), dataset_scaled[1], "b-", linewidth=3., color = '#d3160a')
    ax_s_x.axis("off")
#     ax_s_x.set_xticklabels([])
    ax_s_x.set_xlim((0, sz - 1))
    ax_s_y.plot(- dataset_scaled[0], np.arange(sz)[::-1], "b-", linewidth=3., color = '#054fb4')
    ax_s_y.axis("off")
#     ax_s_y.set_yticklabels([])
#     ax_s_y.invert_xaxis()
#     ax_s_y.invert_yaxis()
    ax_s_y.set_ylim((0, sz - 1))

    plt.colorbar(heatmap, cax = cax)
    plt.title(orf, loc='center', y=1.7, x=-18)
    plt.show()

    
def isoform_pathmap(ts1, ts2, tx1, tx2, timepoints):
    sns.set_style("white", {'xtick.bottom': True})
    sns.set_context("talk")
    x = ts1.loc[[tx1]]
    x.columns = timepoints
    y = ts2.loc[[tx2]]
    y.columns = timepoints
    dataset = x.append(y, ignore_index=True)
    dataset = dataset.to_numpy()
    dataset = to_time_series_dataset(dataset)
    scaler = TimeSeriesScalerMeanVariance(mu=0., std=1.)  # Rescale time series
    dataset_scaled = scaler.fit_transform(dataset)
    sz = dataset_scaled.shape[1]
    path, sim = metrics.dtw_path(dataset_scaled[0], dataset_scaled[1])

    plt.figure(1, figsize=(8, 8))

    # definitions for the axes
    left, bottom = 0.01, 0.1
    w_ts = h_ts = 0.2
    left_h = left + w_ts + 0.02
    left_cb = left + w_ts + 0.2
    width = height = 0.65
    bottom_h = bottom + height + 0.02

    rect_s_y = [left, bottom, w_ts, height]
    rect_gram = [left_h, bottom, width, height]
    rect_s_x = [left_h, bottom_h, width, h_ts]
    rect_cb = [left_cb, bottom+.05, width, .5]

    ax = plt.axes(rect_cb)
    ax_gram = plt.axes(rect_gram)
    ax_s_x = plt.axes(rect_s_x)
    ax_s_y = plt.axes(rect_s_y)

    mat = cdist(dataset_scaled[0], dataset_scaled[1])

    labels = list(range(3,51,3))

    ax_gram.set_yticks(np.arange(mat.shape[0]), minor=False)
    ax_gram.set_xticks(np.arange(mat.shape[1]), minor=False)

    ax_gram.invert_yaxis()
    ax_gram.set_yticklabels(labels, minor=False)
    ax_gram.set_xticklabels(labels, minor=False)
    ax_gram.yaxis.tick_right()

    heatmap = ax_gram.imshow(mat, cmap=plt.cm.Purples_r)
    # ax_gram.axis("off")
    ax_gram.autoscale(False)
    ax_gram.plot([j for (i, j) in path], [i for (i, j) in path], "w-",
             linewidth=3.)

    ax.axis("off")
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size = "5%", pad= 1)

    ax_s_x.plot(np.arange(sz), dataset_scaled[1], "b-", linewidth=3., color = '#d3160a')
    ax_s_x.axis("off")
    ax_s_x.set_xlim((0, sz - 1))
    ax_s_y.plot(- dataset_scaled[0], np.arange(sz)[::-1], "b-", linewidth=3., color = '#054fb4')
    ax_s_y.axis("off")
    ax_s_y.set_ylim((0, sz - 1))

    plt.colorbar(heatmap, cax = cax)
    plt.title(tx1 + "v" + tx2, loc='center', y=1.7, x=-18)
    plt.show()
    
    
def calc_dtw_by_k(df, ref, exp, time_points, k_range):
    column_names = ["ORF", "dtw_score", "path", "k"]
    mrg = pd.DataFrame(columns = column_names)
    for i in k_range:
        res = []
        data = df[df['k'] == i]
        ts1 = data.filter(regex=ref)
        ts2 = data.filter(regex=exp)
        for index, row in ts1.iterrows():
            x = ts1.loc[[index]]
            x.columns = time_points
            y = ts2.loc[[index]]
            y.columns = time_points
            dataset = x.append(y, ignore_index=True)
            dataset = dataset.to_numpy()
            dataset = to_time_series_dataset(dataset)
            scaler = TimeSeriesScalerMeanVariance(mu=0., std=1.)  # Rescale time series
            dataset_scaled = scaler.fit_transform(dataset)
            sz = dataset_scaled.shape[1]
            path, sim = metrics.dtw_path(dataset_scaled[0], dataset_scaled[1])
            res.append((index, sim, path))
    
        scores = pd.DataFrame(data = res, columns=('ORF', 'dtw_score', 'path'))
        scores_k = scores.copy()
        scores_k['k'] = i
        mrg = mrg.append(scores_k)
        mrg = mrg.reset_index(drop=True)
    return mrg

def calc_dtw_sakoe_chiba(ref, exp, timepoints, radius):
    column_names = ["ORF", "dtw_score", "path"]
    df = pd.DataFrame(columns = column_names)
    res = []
    for index, row in ref.iterrows():
        x = ref.loc[[index]]
        x.columns = timepoints
        y = exp.loc[[index]]
        y.columns = timepoints
        dataset = x.append(y, ignore_index=True)
        dataset = dataset.to_numpy()
        dataset = to_time_series_dataset(dataset)
        scaler = TimeSeriesScalerMeanVariance(mu=0., std=1.)  # Rescale time series
        dataset_scaled = scaler.fit_transform(dataset)
        sz = dataset_scaled.shape[1]
        path, sim = metrics.dtw_path(dataset_scaled[0], dataset_scaled[1], global_constraint="sakoe_chiba", sakoe_chiba_radius=radius)
        res.append((index, sim, path))
    df = pd.DataFrame(data = res, columns=column_names)
    return df

## function to create dynamic time warping heatmap and optimal path trace
def pathmap_sakoe_chiba(ref, exp, orf, timepoints, radius):
    sns.set_style("white", {'xtick.bottom': True})
    sns.set_context("talk")
    x = ref.loc[[orf]]
    x.columns = timepoints
    y = exp.loc[[orf]]
    y.columns = timepoints
    dataset = x.append(y, ignore_index=True)
    dataset = dataset.to_numpy()
    dataset = to_time_series_dataset(dataset)
    scaler = TimeSeriesScalerMeanVariance(mu=0., std=1.)  # Rescale time series
    dataset_scaled = scaler.fit_transform(dataset)
    sz = dataset_scaled.shape[1]
    path, sim = metrics.dtw_path(dataset_scaled[0], dataset_scaled[1], global_constraint="sakoe_chiba", sakoe_chiba_radius=radius)

    plt.figure(1, figsize=(8, 8))

    # definitions for the axes
    left, bottom = 0.01, 0.1
    w_ts = h_ts = 0.2
    left_h = left + w_ts + 0.02
    left_cb = left + w_ts + 0.2
    width = height = 0.65
    bottom_h = bottom + height + 0.02

    rect_s_y = [left, bottom, w_ts, height]
    rect_gram = [left_h, bottom, width, height]
    rect_s_x = [left_h, bottom_h, width, h_ts]
    rect_cb = [left_cb, bottom+.05, width, .5]

    ax = plt.axes(rect_cb)
    ax_gram = plt.axes(rect_gram)
    ax_s_x = plt.axes(rect_s_x)
    ax_s_y = plt.axes(rect_s_y)

    mat = cdist(dataset_scaled[0], dataset_scaled[1])

    labels = list(range(3,51,3))

    ax_gram.set_yticks(np.arange(mat.shape[0]), minor=False)
    ax_gram.set_xticks(np.arange(mat.shape[1]), minor=False)

    ax_gram.invert_yaxis()
    ax_gram.set_yticklabels(labels, minor=False)
    ax_gram.set_xticklabels(labels, minor=False)
    ax_gram.yaxis.tick_right()

    heatmap = ax_gram.imshow(mat, cmap=plt.cm.Purples_r)
    # ax_gram.axis("off")
    ax_gram.autoscale(False)
    ax_gram.plot([j for (i, j) in path], [i for (i, j) in path], "w-",
             linewidth=3.)

    ax.axis("off")
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size = "5%", pad= 1)

    ax_s_x.plot(np.arange(sz), dataset_scaled[1], "b-", linewidth=3., color = '#d3160a')
    ax_s_x.axis("off")
    ax_s_x.set_xlim((0, sz - 1))
    ax_s_y.plot(- dataset_scaled[0], np.arange(sz)[::-1], "b-", linewidth=3., color = '#054fb4')
    ax_s_y.axis("off")
    ax_s_y.set_ylim((0, sz - 1))

    plt.colorbar(heatmap, cax = cax)
    plt.title(orf, loc='center', y=1.7, x=-18)
    plt.show()
    
def ridgeplot(df, title):
    sns.set(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})
    # Initialize the FacetGrid object
    pal = sns.color_palette("colorblind")
    g = sns.FacetGrid(df, row="Comparison", hue="Comparison", aspect=15, height=.5, palette=pal)

    g.map(sns.kdeplot, "DTW", clip_on=False, shade=True, alpha=.8, lw=1.5, bw=.2)
    g.map(sns.kdeplot, "DTW", clip_on=False, color="w", lw=2, bw=.2)
    g.map(plt.axhline, y=0, lw=2, clip_on=False)

    # Define and use a simple function to label the plot in axes coordinates
    def label(x, color, label):
        ax = plt.gca()
        ax.text(0, .2, label, fontweight="bold", color=color,
            ha="right", va="center", transform=ax.transAxes)
    
    g.map(label, "DTW")

    # Set the subplots to overlap
    g.fig.subplots_adjust(hspace=-.35)

    # Remove axes details that don't play well with overlap
    g.set_titles("")
    g.set(yticks=[])
    g.despine(bottom=True, left=True)
    plt.gcf().subplots_adjust(left=0.20, bottom=.15)
    plt.xlabel("Dynamic time warping score")
    plt.suptitle(title, fontsize=15)
    
def plot(df, gene):
    sns.set_style("white", {'xtick.bottom': True, 'ytick.left': True})
    sns.set_context("paper")
    data = df[df['ORF'] == gene]
    g = sns.relplot(
        data = data,
        x='HPI', 
        y='TPM', 
        col='ID',
        col_wrap = 4,
        height=3,
        aspect = 1.5/1,
        hue = 'ID',
        kind = 'line',
        palette = 'bright'
        )
    plt.suptitle(gene)
    plt.gcf().subplots_adjust(bottom=0.15, top=0.9)
    g.set(xlim=(3,48),xticks=[3,6,9,12,15,18,21,24,27,30,33,36,39,42,45,48])
    plt.show()

def dtw_df(ts1, ts2, time_points, comp):
    df = calc_dtw(ts1, ts2, time_points)
    df = df.rename(columns = {"dtw_score" : comp})
    df = df.set_index('ORF')
    df = df.drop(labels = 'path', axis = 1)
    return df

def dtw_df_sakoe_chiba(ts1, ts2, time_points, comp, radius):
    df = calc_dtw_sakoe_chiba(ts1, ts2, time_points, radius)
    df = df.rename(columns = {"dtw_score" : comp})
    df = df.set_index('ORF')
    df = df.drop(labels = 'path', axis = 1)
    return df

def getGenes(time_series, rgx, geneList):
    df = time_series.filter(regex=rgx)
    df = df.loc[geneList,:]
    return df

def kruskal_wallis(df):
    kw_data_df = df.drop(columns = "AAvAS_mean", axis=1).reset_index()
    kw_data_df = kw_data_df.melt(id_vars=['ORF'], var_name='Sample', value_name='Shift').set_index('ORF')
    
    kw_data = [kw_data_df.loc[ids, 'Shift'].values for ids in kw_data_df.groupby('Sample').groups.values()]
    H, p = stats.kruskal(*kw_data)
    return kw_data_df, H, p

def calc_dtw_unscaled(ref, exp, timepoints):
    column_names = ["ORF", "dtw_score", "path"]
    df = pd.DataFrame(columns = column_names)
    res = []
    for index, row in ref.iterrows():
        x = ref.loc[[index]]
        x.columns = timepoints
        y = exp.loc[[index]]
        y.columns = timepoints
        dataset = x.append(y, ignore_index=True)
        dataset = dataset.to_numpy()
        dataset = to_time_series_dataset(dataset)
        sz = dataset.shape[1]
        path, sim = metrics.dtw_path(dataset[0], dataset[1])
        res.append((index, sim, path))
    df = pd.DataFrame(data = res, columns=column_names)
    return df

def dtw_df_unscaled(ts1, ts2, time_points, comp):
    df = calc_dtw_unscaled(ts1, ts2, time_points)
    df = df.rename(columns = {"dtw_score" : comp})
    df = df.set_index('ORF')
    df = df.drop(labels = 'path', axis = 1)
    return df