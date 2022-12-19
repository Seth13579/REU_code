from ciao_contrib.runtool import *
import numpy as np
import glob
import matplotlib.pyplot as plt
import os
from extract_counts import make_regions,unglob,download_obsid
from re import sub
import sys

#produces the tricolor light curve plot when given the data
def plot_lc(obsid,pos,time,hards,meds,softs,bin_size,outdir):
    fig,axs = plt.subplots(4,1)

    soft_plt = axs[0]
    med_plt = axs[1]
    hard_plt = axs[2]
    all_plt = axs[3]

    soft_plt.step(time,softs)
    soft_plt.set_title(f'Soft light curve (.3 to 1 KeV)')
    soft_plt.xaxis.set_ticklabels([])
    soft_plt.set(ylabel='Counts/bin')

    med_plt.step(time,meds)
    med_plt.set_title(f'Medium light curve (1 to 2 KeV)')
    med_plt.xaxis.set_ticklabels([])
    med_plt.set(ylabel='Counts/bin')

    hard_plt.step(time,hards)
    hard_plt.set_title(f'Hard light curve (2 to 7.5 KeV)')
    hard_plt.xaxis.set_ticklabels([])
    hard_plt.set(ylabel='Counts/bin')

    all_plt.step(time,softs,label='Soft')
    all_plt.step(time,meds,label='Medium')
    all_plt.step(time,hards,label='Hard')
    all_plt.legend()
    all_plt.set_title(f'Combined light curve')
    all_plt.set(ylabel='Counts/bin',xlabel='Time (s)')


    fig.suptitle(f'J{pos} -- ObsID: {obsid}\nBin size: {bin_size}s')

    #fig.set_size_inches(6,10,forward=True)

    plt.subplots_adjust(hspace=.75,top=.85)

    plt.savefig(f'{outdir}/{pos}_{obsid}_{bin_size}s_tricolor_lc.png',dpi=300)
    plt.show()

    plt.close()

    return


#returns soft,med,hard counts in a dictionary
#binned at the bin size specified
def tricolor_split(evt,reg,bin_size):
    soft_en = '300:1000'
    med_en = '1000:2000'
    hard_en = '2000:7500'


    energies = {'soft':'300:1000',
    'medium': '1000:2000',
    'hard':'2000:7500'
    }

    counts_dict = {}

    for key,val in energies.items():

        dmextract.punlearn()
        dmextract.infile = f'{evt}[energy={val}][sky=region({reg})][bin time=::{bin_size}]'
        dmextract.clobber = 'yes'
        dmextract.opt = 'generic'
        dmextract.outfile = 'temp.fits'
        dmextract()

        dmlist.punlearn()
        dmlist.infile = 'temp.fits[cols counts]'
        dmlist.opt = 'data, clean'

        out = np.asarray(dmlist().split('\n')[1:]).astype('int')

        counts_dict.update({key:out})

    return counts_dict

#the driver function
def make_tricolor_lc(obsid,primary,position,binsize,override=False):
    position_basic = sub('\:','',position)
    working_dir = f'{primary}/{position_basic}'

    if override:
        src_region = input('Enter path to src region file: ').strip()

    else:
        make_regions(obsid,position,working_dir)
        src_region =unglob(glob.glob(f'{working_dir}/*srcreg*'))

    evt = unglob(glob.glob(f'{primary}/*evt2*'))

    src_counts = tricolor_split(evt,src_region,binsize)

    softs = src_counts['soft']
    meds = src_counts['medium']
    hards = src_counts['hard']

    time = [i * float(binsize) for i in range(len(softs))]

    plot_lc(obsid,position_basic,time,hards,meds,softs,binsize,working_dir)



def main():
    if len(sys.argv) == 1:
        pos = input('Enter position of source: ')
        obsid = input('Enter obsid of source: ')
        bin_size = input('Enter bin size in seconds: ')

    else:
        pos = sys.argv[1]
        obsid = sys.argv[2]
        bin_size = sys.argv[3]

    dir = f'./{obsid}/primary'

    override = input('Override regions (be asked to enter your own regions instead of making them automatically)  [y/n]? ')
    if 'y' in override:
        override = True
    else:
        override = False

    download_obsid(obsid)

    make_tricolor_lc(obsid,dir,pos,bin_size,override)


if __name__ == '__main__':
    main()
