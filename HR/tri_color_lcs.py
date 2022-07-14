from ciao_contrib.runtool import *
import numpy as np
import glob
import matplotlib.pyplot as plt
import os
from extract_counts import make_regions,unglob
from re import sub
import sys

#produces the tricolor light curve plot when given the data
def plot_lc(obsid,pos,time,hards,meds,softs,bin_size,outdir):
    fig,axs = plt.subplots(4,1)
    soft_plt = axs[0,0]
    med_plt = axs[1,0]
    hard_plt = axs[2,0]
    all_plt = axs[3,0]

    soft_plt.step(time,softs)
    soft_plt.set_title(f'Soft light curve \n {bin_size}s bins')
    soft_plt.set(ylabel='Counts per bin')

    med_plt.step(time,meds)
    med_plt.set_title(f'Medium light curve \n {bin_size}s bins')
    med_plt.set(ylabel='Counts per bin')

    hard_plt.step(time,hards)
    hard_plt.set_title(f'Medium light curve \n {bin_size}s bins')
    hard_plt.set(ylabel='Counts per bin')

    all_plt.step(time,softs)
    all_plt.step(time,meds)
    all_plt.step(time,hards)
    all_plt.set_title(f'Combined light curve \n {bin_size}')
    all_plt.set(ylabel='Counts per bin',xlabel='Time (s)')

    fig.suptitle(f'J{pos} -- {obsid}')

    plt.savefig(f'{outdir}/{pos}_{obsid}_tricolor_lc.png',dpi=300)
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
def make_tricolor_lc(obsid,primary,position,binsize):
    position_basic = sub('\:','',position)
    working_dir = f'{primary}/{position_basic}'

    make_regions(obsid,position,dir)

    src_region =unglob(glob.glob(f'{working_dir}/*srcreg*'))
    #bkg_region = unglob(glob.glob(f'{working_dir}/*bkgreg.fits'))
    evt = unglob(glob.glob(f'{primary}/*evt2*'))

    src_counts = tricolor_split(evt,src_region,binsize)

    softs = src_counts['soft']
    meds = src_counts['medium']
    hards = src_counts['hard']

    plot_lc(obsid,position_basic,time,hards,meds,softs,bin_size,outdir)



def main():
    pos = sys.argv[1]
    obsid = sys.argv[2]
    dir = f'./{obsid}/primary'

    make_tricolor_lc(obsid,dir,pos,1000)







    #make_tricolor_lc(obsid,dir,pos,1000)

if __name__ == '__main__':
    main()
