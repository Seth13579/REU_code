import numpy as np
import subprocess
import glob
import matplotlib.pyplot as plt
import sys
from extract_counts import unglob

def run_BEHR(bash_file):
    subprocess.run(f'bash {bash_file}', shell = True)

#reads a behr file, returns the median HR and the error bounds
def readfile(file):
    with open(file,'r') as data:
        contents = data.read()
        line = contents.splitlines()[2].split()
        med = line[3]
        lower = line[4]
        upper = line[5]

        return med,upper,lower

def running_average(array,N): #where N is the number of points forward and back to av over
    out = []

    for i in range(len(array)):
        center_value = array[i]
        low_bound = max((0,i-N))
        high_bound = min((len(array)-1,i+N))

        cumsum = sum(array[low_bound:high_bound])

        av = cumsum/(high_bound - low_bound)

        out.append(av)

    return np.array(out).astype('float64')

#takes in the directory with the behr output files
def make_hr_plot(dir,bin_size,position,obsid,show=False,save=True,lines=None):
    print(f'Checking for files in {dir}...')

    files = glob.glob(f'{dir}/*_BEHRresults.txt')

    meds = []
    uppers = []
    lowers = []

    for file in files:
        med,upper,lower = readfile(file)


        meds.append(med)
        uppers.append(upper)
        lowers.append(lower)


    uppers = np.array(uppers).astype('float64')
    lowers = np.array(lowers).astype('float64')
    meds = np.array(meds).astype('float64')

    uppers = uppers - meds
    lowers = meds - lowers

    errors = np.stack((lowers,uppers)).astype('float64')

    x = np.array([int(bin_size) * i for i in range(len(meds))]).astype('float64')

    running_av = running_average(meds,3)

    plt.errorbar(x,meds,yerr=errors,fmt='bo',capsize=2)
    plt.plot(x,running_av,'r-')

    '''
    #now overplot the light curve
    lcfile = unglob(f'{lc_dir}/*{position}_{obsid}*')
    data = np.loadtxt(lcfile, skiprows = 1)
    time = data[::,2]
    counts = data[::,4]
    time = [(i - time[0]) for i in time]

    bin_length = int(bin_size/3.2)
    bin_edges = [time[i] for i in range(len(time)) if i%bin_length == 0]
    bin_edges.append(time[-1])
    binned_counts = binned_statistic(time,counts,statistic='sum',bins=bin_edges)[0]

    plt.step(bin_edges[:-1],binned_counts)
    '''
    plt.ylabel('(H-S)/(H+S)')
    plt.xlabel('Time (s)')
    plt.title(f'HR J{position} ObsID: {obsid} \n {bin_size}s bins')

    if lines is not None:
        ax = plt.gca()
        min,max = ax.get_ylim()

        plt.vlines(lines,ymin=min,ymax=max,colors='r',linestyles='dotted')

    if save:
        plt.savefig(f'{position}_{obsid}_HR.png')

    if show:
        plt.show()
