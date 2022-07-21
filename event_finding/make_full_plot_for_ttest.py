import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import binned_statistic
import sys
import warnings
warnings.filterwarnings("ignore")
#making a plot with 4 subplots: the cumo counts, the auto binned lc, the
#1000 s light curve and the 2000 s light curve

def bin_array(xs,ys,binsize):
    out_y =  []
    out_x = []
    y_sum = 0
    x_sum = 0

    for i in range(len(xs)):
        x = xs[i]
        y = ys[i]

        if i == 0:
            last_x = 0
        else:
            last_x = xs[i-1]

        x_sum += (x-last_x)
        y_sum += y

        if x_sum >= binsize:
            out_y.append(y_sum)
            out_x.append(x)

            x_sum = 0
            y_sum = 0
    if y_sum != 0:
        out_y.append(y_sum)
        out_x.append(x)

    return out_x,out_y


def make_full_plot(file,outdir,save=True,show=False,lines=None,trim = True):
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

    if trim:
        #trim the ends off, if they have 0 counts
        i = 0
        while counts[i] == 0:
            i += 1

        if counts[-1] == 0:
            #print('Trimming off the end')
            j = -1
            while counts[j] == 0:
                j -= 1
            j +=1
        else:
            j = None

        if j is not None:
            time = time[i:j]
            counts = counts[i:j]
            cumo = cumo[i:j]
        else:
            time = time[i:]
            counts = counts[i:]
            cumo = cumo[i:]






    #bin the light curve data
    #bin sizes: auto, 1000 seconds, 2000 seconds
    #auto_binsize = max(500,5/(sum(counts)/time[-1]))
    auto_binsize = 500

    #start the binning:
    bin_length_auto = auto_binsize
    bin_length_1 = 1000
    bin_length_2 = 2000


    bin_edges_auto,binned_counts_auto = bin_array(time,counts,bin_length_auto)
    bin_edges_1,binned_counts_1 = bin_array(time,counts,bin_length_1)
    bin_edges_2,binned_counts_2 = bin_array(time,counts,bin_length_2)

    fig, axs = plt.subplots(2,2)
    cumo_plt = axs[0,0]
    auto_plt = axs[0,1]
    bin1_plt = axs[1,0]
    bin2_plt = axs[1,1]

    cumo_plt.plot(time,cumo)
    cumo_plt.set_title('Cumulative Counts')

    auto_plt.step(bin_edges_auto,binned_counts_auto)
    auto_plt.set_title(f'Bin size: {round(auto_binsize,2)}s')

    bin1_plt.step(bin_edges_1,binned_counts_1)
    bin1_plt.set_title('Bin size: 1000s')

    bin2_plt.step(bin_edges_2,binned_counts_2)
    bin2_plt.set_title('Bin size: 2000s')

    cumo_plt.set(ylabel='Counts')
    bin1_plt.set(xlabel='Time (s)',ylabel='Counts')
    bin2_plt.set(xlabel='Time (s)')

    fig.suptitle(trunc_name)

    if lines is not None:
        print(f'Plotting lines at: {lines}')

        for ax in axs.flat:
            ax.vlines(lines,ymin=ax.get_ylim()[0],ymax=ax.get_ylim()[1],colors='r',linestyles='dotted')
            #ax.hlines()
            # TODO: Add in hline for the mean


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

    file = sys.argv[1]
    #trim = sys.argv[2]

    '''
    if trim == 'true' or trim == 'True':
        trim = True
    else:
        trim = False
    '''

    make_full_plot(file,'.',save=False,show=True,trim = False)
