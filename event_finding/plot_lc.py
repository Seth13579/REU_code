#S. Larner
#Plots a light curve, with binning

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import binned_statistic
import sys


def plot_lc(file,binsize=None,show=False,save=True,lines=None,errorbars=True):
    if not show and not save:
        return

    trunc_name = file.split('/')[-1].split('.')[:-1]
    trunc_name = ''.join(trunc_name)

    #read the file
    data = np.loadtxt(file, skiprows = 1)
    time = data[::,2]
    counts = data[::,4]
    count_rate_err = data[::,9]

    #adjust the time
    time = [(i - time[0]) for i in time]

    #trim the ends off, if they have 0 counts
    i = 0
    while counts[i] == 0:
        i += 1

    j = -1
    while counts[j] == 0:
        j -= 1

    time = time[i:j]
    counts = counts[i:j]

    if binsize == 'auto':
        c = sum(counts)
        t = time[-1]
        av_rate = c/t

        #on average, 5 counts per bin. Change as needed
        binsize = int(5/av_rate)

        binsize = max(500,binsize)

    #the length of the bins in the array
    bin_length = int(binsize/3.2)
    bin_edges = [time[i] for i in range(len(time)) if i%bin_length == 0]
    bin_edges.append(time[-1])
    binned_counts = binned_statistic(time,counts,statistic='sum',bins=bin_edges)[0]

	#currently not accurate, need to fix
    if errorbars:
        #bin the errors
        binned_errs = []
        for count,i in enumerate(range(0,len(time),bin_length)):
            next_edge = i + bin_length

            next_edge = min(len(count_rate_err),next_edge)

            err_to_bin = count_rate_err[i:next_edge]
            rate_to_bin = count_rate[i:next_edge]

            to_sum = []
            for i in range(len(err_to_bin)):
                if rate_to_bin[i] != 0:
                    to_sum.append((err_to_bin[i]/rate_to_bin[i])**2)

            binned = np.sqrt(sum(to_sum)) * binned_counts[count]

            binned_errs.append(binned)

        #fix the lower error being below zero
        binned_errs_lower = [binned_counts[i] if binned_counts[i]-binned_errs[i] < 0 else binned_errs[i] for i in range(len(binned_counts)) ]

        binned_errs_asym = np.array([binned_errs_lower,binned_errs])

    else:
        binned_errs_asym = None

    plt.errorbar(bin_edges[:-1],binned_counts,yerr=binned_errs_asym,capsize= 2)
    plt.xlabel('Time since obs start (s)')
    plt.ylabel('Counts in bin')
    plt.title(f'Light Curve {trunc_name}\nBin size: {binsize} s')

    if lines is not None:
        plt.vlines(lines,ymin=0,ymax=max(binned_counts),colors = 'r',linestyles='dotted')

    if show:
        plt.show()

    if save:
        plt.savefig(f'{trunc_name}_lc.pdf')

    plt.close()

    return

def main():
    file = sys.argv[1]
    binsize = float(sys.argv[2])
    plot_lc(file,binsize,show=True,errorbars=False)

if __name__ == '__main__':
    main()
