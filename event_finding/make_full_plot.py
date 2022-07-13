import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import binned_statistic
import sys

#making a plot with 4 subplots: the cumo counts, the auto binned lc, the
#1000 s light curve and the 2000 s light curve

def make_full_plot(file,outdir,save=True,show=False,lines=None):
    if not show and not save:
        return

    trunc_name = file.split('/')[-1].split('.')[:-1]
    trunc_name = '.'.join(trunc_name)

    #read in the data from the file
    data = np.loadtxt(file, skiprows = 1)
    counts = data[::,4]

    '''
    if sum(counts) == 0:
        print('No counts')
        return
    '''

    time = data[::,2]

    #adjust the time
    time = [(i - time[0]) for i in time]

    #make the cumulative counts array
    cumo = np.cumsum(counts)

    #trim the ends off, if they have 0 counts
    i = 0
    while counts[i] == 0:
        i += 1

    j = -1
    while counts[j] == 0:
        j -= 1

    time = time[i:j]
    counts = counts[i:j]
    cumo = cumo[i:j]

    #bin the light curve data
    #bin sizes: auto, 1000 seconds, 2000 seconds
    auto_binsize = max(500,5/(sum(counts)/time[-1]))

    #start the binning:
    bin_length_auto = int(auto_binsize/3.2)
    bin_length_1 = int(1000/3.2)
    bin_length_2 = int(2000/3.2)

    bin_edges_auto = [time[i] for i in range(len(time)) if i%bin_length_auto == 0]
    bin_edges_1 = [time[i] for i in range(len(time)) if i%bin_length_1 == 0]
    bin_edges_2 = [time[i] for i in range(len(time)) if i%bin_length_2 == 0]

    bin_edges_auto.append(time[-1])
    bin_edges_1.append(time[-1])
    bin_edges_2.append(time[-1])

    binned_counts_auto = binned_statistic(time,counts,statistic='sum',bins=bin_edges_auto)[0]
    binned_counts_1 = binned_statistic(time,counts,statistic='sum',bins=bin_edges_1)[0]
    binned_counts_2 = binned_statistic(time,counts,statistic='sum',bins=bin_edges_2)[0]

    fig, axs = plt.subplots(2,2)
    cumo_plt = axs[0,0]
    auto_plt = axs[0,1]
    bin1_plt = axs[1,0]
    bin2_plt = axs[1,1]

    cumo_plt.plot(time,cumo)
    cumo_plt.set_title('Cumulative Counts')

    auto_plt.step(bin_edges_auto[:-1],binned_counts_auto)
    auto_plt.set_title(f'Bin size: {round(auto_binsize,2)}s')

    bin1_plt.step(bin_edges_1[:-1],binned_counts_1)
    bin1_plt.set_title('Bin size: 1000s')

    bin2_plt.step(bin_edges_2[:-1],binned_counts_2)
    bin2_plt.set_title('Bin size: 2000s')

    cumo_plt.set(ylabel='Counts')
    bin1_plt.set(xlabel='Time (s)',ylabel='Counts')
    bin2_plt.set(xlabel='Time (s)')

    fig.suptitle(trunc_name)

    if lines is not None:
        for ax in axs.flat:
            ax.vlines(lines,ymin=ax.get_ylim()[0],ymax=ax.get_ylim()[1],colors='r',linestyles='dotted')

    plt.subplots_adjust(hspace=.3)

    if save:
        fig.savefig(f'{outdir}/{trunc_name}_full_lc.png')
        fig.savefig(f'{outdir}/{trunc_name}_full_lc.pdf')

    if show:
        plt.show()

    plt.close()

    return

if __name__ == '__main__':
    '''
    import glob

    for file in glob.iglob('./*.txt'):
        make_full_plot(file,'.',save=False,show=True,lines = None)
    '''

    file = './lc_test.txt'
    make_full_plot(file,'.',save=False,show=True)
