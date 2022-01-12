import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt
import scipy
import seaborn as sns
import matplotlib

##### Taken from Haase Lab git ####

#make colormap haase lab color scheme
norm = matplotlib.colors.Normalize(-1.5,1.5)
colors = [[norm(-1.5), "cyan"],
          [norm(0), "black"],
         [norm(1.5), "yellow"]]
haase = matplotlib.colors.LinearSegmentedColormap.from_list("", colors)


## function to plot data in a heat map, sorted based off of the maximum in the first period, and normalized by z-score
def heatmap_max(data, gene_list, first_period, yticks = False, cbar_bool = True, axis = None):
    # data = data frame with gene names as the index
    # gene list = list of gene names to plot
    # first period = integer value of the last column of the first period

    # yticks: default False doesn't plot ytick labels
    # cbar_bool: default True does include color bar
    # axis: default None doesn't specify an axis
    data = data.loc[gene_list]
    z_pyjtk = scipy.stats.zscore(data, axis=1)
    z_pyjtk_df = pd.DataFrame(z_pyjtk, index=data.index, columns=data.columns)
    z_pyjtk_1stperiod = z_pyjtk_df.iloc[:, 0:first_period]
    max_time = z_pyjtk_1stperiod.idxmax(axis=1)
    z_pyjtk_df["max"] = max_time
    z_pyjtk_df["max"] = pd.to_numeric(z_pyjtk_df["max"])
    z_pyjtk_df = z_pyjtk_df.sort_values(by="max", axis=0)
    z_pyjtk_df = z_pyjtk_df.drop(columns=['max'])
    order = z_pyjtk_df.index
    s = sns.heatmap(z_pyjtk_df, cmap=haase, vmin=-1.5, vmax=1.5, yticklabels = yticks, cbar = cbar_bool, ax = axis)
    return order


## function to plot data in a heat map, sorted based off of the supplied order, and normalized by z-score
def heatmap_order(data, order, yticks = False, cbar_bool = True, axis = None):
    # data = data frame with gene names as the index
    #order = list of genes to plot in desired order

    # yticks: default False doesn't plot ytick labels
    # cbar_bool: default True does include color bar
    # axis: default None doesn't specify an axis
    data = data.loc[order]
    z_pyjtk = scipy.stats.zscore(data, axis=1)
    z_pyjtk_df = pd.DataFrame(z_pyjtk, index=data.index, columns=data.columns)
    order2 = z_pyjtk_df.index
    s = sns.heatmap(z_pyjtk_df, cmap=haase, vmin=-1.5, vmax=1.5, yticklabels = yticks, cbar= cbar_bool,ax = axis)
    return order2

## function to plot data in a heat map, sorted based off of the left edge, and normalized by z-score
def heatmap_LE(data,gene_list, first_period, yticks = False, cbar_bool = True, axis = None):
    #adapted from R code orderHeatmapsLE (from R function autoorder)

    # data = data frame with gene names as the index
    # gene list = list of gene names to plot
    # first period = integer value of the last column of the first period

    # yticks: default False doesn't plot ytick labels
    # cbar_bool: default True does include color bar
    # axis: default None doesn't specify an axis
    data = data.loc[gene_list]
    z_pyjtk = scipy.stats.zscore(data, axis=1)
    z_pyjtk_df = pd.DataFrame(z_pyjtk, index=data.index, columns=data.columns)
    z_pyjtk_1stperiod = z_pyjtk_df.iloc[:, 0:first_period]
    LE = pd.DataFrame(index=z_pyjtk_1stperiod.index)
    LE["order"] = 0
    for i in range(0, len(z_pyjtk_1stperiod.index)):
        temp = z_pyjtk_1stperiod.iloc[i]
        sPos = -1
        maxLength = 0
        for j in range(0, len(temp)):
            if temp.iloc[j] > 1:
                if sPos == -1:
                    sPos = j
                    curLength = 0
                curLength = curLength + 1
            if temp.iloc[j] < 1:
                if sPos != -1:
                    if curLength > maxLength:
                        maxLength = curLength
                        LE.iloc[i] = sPos
        if sPos != -1:
            if curLength > maxLength:
                maxLength = curLength
                LE.iloc[i] = sPos
    z_pyjtk_df["order"] = LE["order"]
    z_pyjtk_df = z_pyjtk_df.sort_values(by="order", axis=0)
    z_pyjtk_df = z_pyjtk_df.drop(columns=['order'])
    order = z_pyjtk_df.index
    s = sns.heatmap(z_pyjtk_df, cmap=haase, vmin=-1.5, vmax=1.5, yticklabels = yticks, cbar = cbar_bool, ax = axis )

    return order

#########

def get_results(results_dir):
    # reads in pydl and pyjtk results into dfs for one ts
    dir_contents = os.listdir(results_dir)
    node_dir = [item for item in dir_contents if item.startswith("node")][0]
    node_dir_contents = os.listdir(os.path.join(results_dir, node_dir))
    pydl_file = [item for item in node_dir_contents if item.startswith("pydl")][0]
    pyjtk_file = [item for item in node_dir_contents if item.startswith("pyjtk")][0]
    pydl_df = pd.read_csv(os.path.join(results_dir, node_dir, pydl_file), sep="\t", comment="#", index_col=0)
    pyjtk_df = pd.read_csv(os.path.join(results_dir, node_dir, pyjtk_file), sep="\t", comment="#", index_col=0)
    return pydl_df, pyjtk_df

def get_ts(filename):
    # reads ts into df
    ts_df = pd.read_csv(os.path.join("Data", filename), sep="\t", index_col=0)
    return(ts_df)

def pydl_filter_ts(ts_df, pydl_df, threshold, thres_rank=False):
    # filters pydl df on dl_score column, or if thres_rank == True, filter by top N genes
    pydl_sorted = pydl_df.sort_values(by="dl_score")
    if thres_rank:
        pydl_sorted_thres = pydl_sorted[:threshold]
    else:
        pydl_sorted_thres = pydl_sorted[pydl_sorted["dl_score"] <= threshold]
    ts_filtered_df = ts_df[ts_df.index.isin(pydl_sorted_thres.index)]
    return ts_filtered_df

def pyjtk_filter_ts(ts_df, pyjtk_df, threshold, thres_rank=False):
    # filters pyjtk df on p-value column, or if thres_rank == True, filter by top N genes
    pyjtk_sorted = pyjtk_df.sort_values(by="p-value")
    if thres_rank:
        pyjtk_sorted_thres = pyjtk_sorted[:threshold]
    else:  
        pyjtk_sorted_thres = pyjtk_sorted[pyjtk_sorted["p-value"] <= threshold]
    ts_filtered_df = ts_df[ts_df.index.isin(pyjtk_sorted_thres.index)]
    return ts_filtered_df

def plot_heatmap_max(ts_df, results_df, method, thres, thres_rank=False, axis=None):
    # plot heatmap of a ts that has been ordered by phase of max gene expression in first period.
    # ts has also been filtered by a given threshold used on a given method of ranking periodicity
#     f, axs = plt.subplots(figsize=(8,12))
    if method == "pydl":
        genelist = pydl_filter_ts(ts_df, results_df, thres, thres_rank).index.tolist()
    elif method == "pyjtk":
        genelist = pyjtk_filter_ts(ts_df, results_df, thres, thres_rank).index.tolist()
        
    heat = heatmap_max(ts_df, genelist, 16, axis=axis) # Remove hard coded int
    return(heat)