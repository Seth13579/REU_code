#a script for running behr but with bins of fixed count width instead of
#fixed time width

import numpy as np
from ciao_contrib.runtool import *
from create_behr_files import region_area
import os
from run_BEHR import readfile
import glob
import sys
import matplotlib.pyplot as plt
from ciao_contrib.runtool import *
#from  extract_counts import unglob


#returns an array of soft and an array of hard counts (as ordered)
#among the +/- N at all times
def constant_count_split(evt,reg,bkgreg,divide_energy,N=16):

    try:
        dmlist.punlearn()
        dmlist.infile = f'{evt}[sky=region({reg})][cols time, energy]'
        dmlist.opt = 'data,clean'
        raw_out = np.genfromtxt(dmlist().split('\n'),dtype='str').astype('float64')

    except OSError:
        dmlist.punlearn()
        dmlist.infile = f'{evt}[sky={reg}][cols time, energy]'
        dmlist.opt = 'data,clean'
        raw_out = np.genfromtxt(dmlist().split('\n'),dtype='str').astype('float64')


    times = raw_out[::,0]
    energies = raw_out[::,1]

    dmlist.punlearn()

    try:
        dmlist.infile = f'{evt}[sky=region({bkgreg})][cols time, energy]'
        dmlist.opt = 'data,clean'

        raw_out_bkg = np.genfromtxt(dmlist().split('\n'),dtype='str').astype('float64')
    except OSError:
        dmlist.infile = f'{evt}[sky={bkgreg}][cols time, energy]'
        dmlist.opt = 'data,clean'

        raw_out_bkg = np.genfromtxt(dmlist().split('\n'),dtype='str').astype('float64')


    bkg_times = raw_out_bkg[::,0]
    bkg_energies = raw_out_bkg[::,1]

    softs = []
    hards = []
    bkg_softs = []
    bkg_hards = []

    soft_count = 0
    hard_count = 0

    soft_bkg_count = 0
    hard_bkg_count = 0

    for i in range(len(times)):
        low_bound = max((0,i-N))
        high_bound = min((len(times)-1,i+N))
        energy_slice = energies[low_bound:high_bound]

        #the first time, we have to count them all
        if i == 0:
            for count in energy_slice:
                if count <= divide_energy:
                    soft_count += 1
                else:
                    hard_count += 1

        else:
            #general case, both ends moving
            if i >= N and i+N <= len(times)-1:
                new_count = energy_slice[-1]
                old_count = energies[low_bound - 1]

                if old_count <= divide_energy:
                    soft_count -= 1
                else:
                    hard_count -= 1

                if new_count <= divide_energy:
                    soft_count += 1
                else:
                    hard_count += 1

            #left edge: adding a new count but not losing any
            elif i < N:
                new_count = energy_slice[-1]

                if new_count <= divide_energy:
                    soft_count += 1
                else:
                    hard_count += 1

            #right edge: losing an old count but not gaining any
            #elif i+N <= len(times)-1:
            else:
                old_count = energies[low_bound - 1]

                if old_count <= divide_energy:
                    soft_count -= 1
                else:
                    hard_count -= 1

        #now to deal with the background:
        #need to take only the background counts from inside the time range
        #that the source counts are taken from
        bkg_energies_bool = [1 if en < divide_energy else 0 for en in bkg_energies]

        low_time = times[low_bound]
        high_time = times[high_bound]

        bkg_low = low_bound
        bkg_high = high_bound

        for i,time in enumerate(bkg_times):
            if time < low_time:
                bkg_low = i + 1

            elif time == low_time:
                bkg_low = i

            elif time == high_time:
                bkg_high = i
                break

            elif time > high_time:
                bkg_high = i - 1
                break

        bkg_energy_slice = bkg_energies_bool[bkg_low:bkg_high]

        bkg_softs.append(max(sum(bkg_energy_slice),0))
        bkg_hards.append(max(0,len(bkg_energy_slice)-sum(bkg_energy_slice)))


        softs.append(max(soft_count,0))
        hards.append(max(hard_count,0))

    softs = np.array(softs).astype('int')
    hards = np.array(hards).astype('int')
    bkg_softs = np.array(bkg_softs).astype('int')
    bkg_hards = np.array(bkg_hards).astype('int')

    bkg_soft_nonzero = [i for i in bkg_softs if i != 0]
    bkg_hard_nonzero = [i for i in bkg_hards if i != 0]

    assert len(bkg_soft_nonzero) != 0
    assert len(bkg_hard_nonzero) != 0

    all = np.column_stack((softs,hards,bkg_softs,bkg_hards))

    np.save('debugging.npy',all)
    np.savetxt('debugging.txt',all,fmt='%s',header='#softs,hards,bkg_softs,bkg_hards',delimiter=',')

    return softs,hards,bkg_softs,bkg_hards


def make_behr(evt,srcreg,bkgreg,divide_energy,BEHR_DIR,outfile,BEHR_outdir,N=16,confidence='68.00'):
    soft_src,hard_src,soft_bkg,hard_bkg = constant_count_split(evt,srcreg,bkgreg,divide_energy,N)

    hard_area = region_area(evt,bkgreg,3000)/region_area(evt,srcreg,3000)
    soft_area = region_area(evt,bkgreg,1000)/region_area(evt,srcreg,1000)

    '''
    if os.path.exists(outfile):
        cont = ''
        while 'y' not in cont and 'n' not in cont:
            cont = input('BEHR outfile exists. Proceeding would overwrite previous work. \n Continue? (y/n): ')
            if 'n' in cont:
                raise Exception
    '''

    with open(outfile,'w') as writeto:
        writeto.write(f'cd {BEHR_DIR}')
        for i in range(len(soft_src)):
            writeto.write(f'\n echo "softsrc={soft_src[i]} hardsrc={hard_src[i]}   softbkg={soft_bkg[i]}   hardbkg={hard_bkg[i]}"')
            writeto.write(f'\n echo {i} of {len(soft_src)}')
            writeto.write(f'\n./BEHR softsrc={soft_src[i]} hardsrc={hard_src[i]}   softbkg={soft_bkg[i]}   hardbkg={hard_bkg[i]}   softarea={soft_area} hardarea={hard_area} output={BEHR_outdir}/{i}_BEHRresults level={confidence}')

def plot_BEHR_constcounts(dir,bin_size,position,obsid,evt,reg,start_time,show=False,
                        save=True,lines=None):

    dmlist.punlearn()
    dmlist.infile = f'{evt}[sky=region({reg})][cols time]'
    dmlist.opt = 'data,clean'
    try:
        x = np.array(dmlist().split('\n'))[2:]
    except OSError:
        dmlist.infile = f'{evt}[sky={reg}][cols time]'
        x = np.array(dmlist().split('\n'))[2:]

    try:
        x = x.astype('float64')
    except Exception as e:
        print(x)
        raise e


    x -= start_time

    meds = []
    uppers = []
    lowers = []

    for i in range(len(x)):
        file = f'{dir}/{i}_BEHRresults.txt'
        med,upper,lower = readfile(file)


        meds.append(med)
        uppers.append(upper)
        lowers.append(lower)

    uppers = np.array(uppers).astype('float64')
    lowers = np.array(lowers).astype('float64')
    meds = np.array(meds).astype('float64')

    if save or show:
        plt.plot(x,meds,'k-')
        plt.fill_between(x,lowers,uppers,step='mid')

        plt.ylabel('(H-S)/(H+S)')
        plt.xlabel('Time (s)')
        plt.title(f'HR J{position} ObsID: {obsid} \n Each point ± {bin_size} counts')

        if lines is not None:
            ax = plt.gca()
            min,max = ax.get_ylim()

            plt.vlines(lines,ymin=min,ymax=max,colors='r',linestyles='dotted')

        if save:
            plt.savefig(f'{position}_{obsid}_HR_constcounts.png',dpi=300)

        if show:
            plt.show()

        plt.close()

    return x,uppers,lowers,meds




if __name__ == '__main__':
    evt = 'test_evt.fits'
    reg = 'test_reg.fits'
    divide_energy = 1000

    constant_count_split(evt,reg,divide_energy)
