import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc
import matplotlib.ticker as ticker
from scipy.signal import find_peaks
from numpy.random import seed
from numpy.random import randn
from numpy import mean
from numpy import std
from scipy import stats
from statsmodels.distributions.empirical_distribution import ECDF
from scipy.stats import ks_2samp
from scipy.stats import normaltest
from scipy.stats import mannwhitneyu
import statsmodels.api as sa
import statsmodels.formula.api as sfa
import scikit_posthocs as sp
from venn import venn

def kruskal_wallis(df):
    kw_data_df = df.drop(columns = "AAvAS_mean", axis=1).reset_index()
    kw_data_df = kw_data_df.melt(id_vars=['ORF'], var_name='Sample', value_name='Shift').set_index('ORF')
    
    kw_data = [kw_data_df.loc[ids, 'Shift'].values for ids in kw_data_df.groupby('Sample').groups.values()]
    H, p = stats.kruskal(*kw_data)
    return kw_data_df, H, p

def getIsolate(time_series, rgx):
    df = time_series.filter(regex=rgx)
    time_points = list(range(3, 51, 3))
    df.columns = time_points
    return df

def get_ts(directory, filename):
    # reads ts into df
    ts_df = pd.read_csv(os.path.join(directory, filename), sep="\t", index_col=0)
    return(ts_df)

def onePeakGenes(df, threshold, distance, peak_width, peak_prominence):
    orf_list = []
    orfs = df.index.tolist()
    rows = df.to_numpy()
    peaksData = dict(zip(orfs, rows)) 
    
    for key, value in peaksData.items():
        peaks, properties = find_peaks(value, threshold = threshold, distance = distance, width=peak_width, prominence=peak_prominence)
        if (len(peaks) == 1):
            orf_list.append(key)
    return orf_list

def make_peak_df(strain_ts, ts1_regex, ts1_id, ts2_regex, ts2_id, timepoints):
    ts1 = getIsolate(strain_ts, ts1_regex)
    ts2 = getIsolate(strain_ts, ts2_regex)
    ts1.columns = timepoints
    ts2.columns = timepoints
    new_df = pd.DataFrame(columns = [ts1_id, ts2_id]) 
    new_df[ts1_id] = ts1.idxmax(axis=1)
    new_df[ts2_id] = ts2.idxmax(axis=1)
    return new_df

def calc_peak_change(dataframe, col1, col2):
    counts = pd.DataFrame()
    counts_agg = pd.DataFrame()
    counts_trans = pd.DataFrame()
    dataframe_norm = pd.DataFrame()
    
    dataframe['peak_time_change'] = dataframe[col1] - dataframe[col2]
    dataframe = dataframe[abs(dataframe['peak_time_change']) < 45]
    
    counts = dataframe.groupby([col1, 'peak_time_change']).peak_time_change.agg('count').to_frame('Percent of transcripts').reset_index()
    counts_agg = counts.groupby([col1, 'peak_time_change']).agg({'Percent of transcripts': 'sum'})
    counts_trans = counts.groupby([col1]).agg({'Percent of transcripts': 'sum'})
    dataframe_norm = counts_agg.div(counts_trans, level = col1) * 100
    dataframe_norm = dataframe_norm.reset_index()
    return dataframe, dataframe_norm

def plot_peakShifts(df, ref_col, title):

    y_range = list(range(-48,48,3))
    x_range = list(range(3,51,3))
    sns.set_context("talk")
    sns.set_style("white", {"axes.facecolor": "1", 'axes.edgecolor': '.1', 'xtick.bottom': True, 'ytick.left': True, 'grid.linestyle': '--'} )
    g = sns.JointGrid(x=ref_col, y="peak_time_change", data=df, height=10, ratio=4, ylim = (-48,48))

    g = g.plot_joint(sns.scatterplot, palette = 'coolwarm', hue = df['peak_time_change'], size = df['Percent of transcripts'], sizes=(100, 500), legend='full') #color="black",

    _ = g.ax_marg_x.hist(df[ref_col], color="black", alpha=.6,
                     bins=np.arange(1.5, 52, 3))

    #18 hour max
    # _ = g.ax_marg_y.hist(max_df_FUP_pct["peak_time_change"], color="black", alpha=.6, align='right',
    #                      orientation="horizontal", bins=np.arange(-20.6, 21, 3))

    #24 hour max
#     _ = g.ax_marg_y.hist(df["peak_time_change"], color="black", alpha=.6, align='right',
#                          orientation="horizontal", bins=np.arange(-26.6, 25, 3))
    # No max
    _ = g.ax_marg_y.hist(df["peak_time_change"], color="black", alpha=.6, align='right',
                         orientation="horizontal", bins=np.arange(-50.6, 49, 3))
    
    g.set_axis_labels("Hours post-invasion", "Hours between peak expression time in \n {}".format(title))
    g.ax_joint.yaxis.labelpad = 45
    plt.setp(g.ax_marg_x.get_yticklabels(), visible=True)
    plt.setp(g.ax_marg_y.get_xticklabels(), visible=True)
    g.ax_joint.set(yticks=y_range , xticks=x_range)
    g.ax_joint.yaxis.grid()
    g.ax_joint.yaxis.set_major_locator(ticker.MultipleLocator(6))
    g.ax_joint.set(xlim=(1.5,50))
    g.ax_joint.legend().remove()
    # plt.legend(bbox_to_anchor=(1.4, 1), loc=2, borderaxespad=0., prop={'size': 30})
    plt.subplots_adjust(top=0.9)
    g.fig.suptitle('Shifts in peak transcription time between \n {}'.format(title), x=.5)
    plt.gcf().subplots_adjust(left = .2)
    # plt.savefig('Figures/PeakExpression/FUP_all_tx.png', dpi=300)

def plotkde(df):
    data = df['peak_time_change'].to_numpy()
    sns.kdeplot(data, cut=0, shade=True, bw=1.5)
    
def getStageData(data, ref, ring_end, troph_end):
    ring = data[data[ref] < ring_end]
    troph = data[(data[ref] >= ring_end) & (data[ref] <= troph_end)]
    schiz = data[data[ref] > troph_end]
    return ring, troph, schiz
    
    
def calcCorr(data1, data2, description):
    x = data1.to_numpy()
    y = data2.to_numpy()
    print ("Spearman and Pearson correlations for {}".format(description))
    print(scipy.stats.spearmanr(x, y))
    pearson = scipy.stats.pearsonr(x, y)
    print("Pearson: {}".format(pearson)) 

def testNormality(data, description):
    stat, p = normaltest(data)
    print('Data: {}'.format(description))
    print('Statistics=%.3f, p=%.3f' % (stat, p))
    # interpret
    alpha = 0.05
    if p > alpha:
        print('{} look Gaussian (fail to reject H0)'.format(description))
    else:
        print('{} do not look Gaussian (reject H0)\n'.format(description))
        
def kolmogorov_smirnov(data_list, names_list, alt_hyp):
    results_dict = []
    results_df = pd.DataFrame()
    for i in range(len(data_list)):
        for j in range(len(data_list)):
            if (j > i):
                ks, p = ks_2samp(data_list[i], data_list[j], alternative = alt_hyp)
                results = {'F1': names_list[i], 'F2': names_list[j], 'KS': ks, 'p-value': p}
                results_dict.append(results)
    results_df = pd.DataFrame.from_dict(results_dict)
    return results_df        
        
def testMannWhitney(test, x, y, alt):
    print(test)
    stat, p = mannwhitneyu(x, y, alternative=alt)
    print('Statistics=%.3f, p=%.5f' % (stat, p))
    # interpret
    alpha = 0.05
    if p > alpha:
        print('Same distribution (fail to reject H0)\n\n')
    else:
        print('Different distribution (reject H0)\n\n')
        
def plotECDF(data, names, comparison):
    sns.set_context("paper")
    colors = ["black", "red", "green", "blue", "orange", "gray"]
    fig, axs = plt.subplots(nrows=2, ncols=3, constrained_layout=True)
    fig.suptitle('{}'.format(comparison))
    k = 0
    for i in range(0,2):
        for j in range(0,3):
            axs[i, j].set_title('{}'.format(names[k]))
            axs[i, j].set_xlim([-24, 24])
            k+=1
    i = 0
    for row in axs:
        for col in row:
            col.plot(data[i].x, data[i].y, color = colors[i])
            i+=1

def writeList(dataList, col, path):
    import csv
    mylist= dataList
    with open(path, 'w', newline='') as myfile:
        wr = csv.writer(myfile, delimiter=',')
        header = [col] 
        wr.writerow(header) 
        for val in mylist:
            wr.writerow([val])

def ridgeplot(df, title):
    sns.set(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})
    # Initialize the FacetGrid object
    pal = sns.color_palette("colorblind")
    g = sns.FacetGrid(df, row="Comparison", hue="Comparison", aspect=15, height=.5, palette=pal)

    g.map(sns.kdeplot, "peakShift", clip_on=False, shade=True, alpha=.8, lw=1.5, bw=2)
    g.map(sns.kdeplot, "peakShift", clip_on=False, color="w", lw=2, bw=2)
    g.map(plt.axhline, y=0, lw=2, clip_on=False)

    # Define and use a simple function to label the plot in axes coordinates
    def label(x, color, label):
        ax = plt.gca()
        ax.text(0, .2, label, fontweight="bold", color=color,
            ha="right", va="center", transform=ax.transAxes)
    
    g.map(label, "peakShift")

    # Set the subplots to overlap
    g.fig.subplots_adjust(hspace=-.35)

    # Remove axes details that don't play well with overlap
    g.set_titles("")
    g.set(yticks=[])
    g.despine(bottom=True, left=True)
    plt.gcf().subplots_adjust(left=0.20, bottom=.15)
    plt.xlabel("Peak shift")
    plt.suptitle(title, fontsize=15)