#A script for phase folding the hardness ratios
#python ../HR_phase_fold.py 09:56:08.06+69:41:39.59 10 10542,10543,10544,10545,16580 221260 362283895.21616

from HR_driver_const_counts import main
import numpy as np
from add_asym import add_asym
import sys
import matplotlib.pyplot as plt
import sys

class HR_graph:
    """docstring for HR_graph.
    Contains all the information from HR_driver_const_counts in a machine
    readable format
    """

    def __init__(self, text):
        #the full text from HR driver
        self.text = text

        self.times = text[::,0]

        self.uppers = text[::,1]

        self.lowers = text[::,2]

        self.meds = text[::,3]

        self.phase = None

    def make_phase(self,period,t0):
        phase = [((i-t0)%period) / period for i in self.times]

        self.phase = phase

        return


class Phase_Space_Bins:
    '''
    Represent the bins in phase space that the points are sorted into
    '''

    def __init__(self,min_phase,max_phase,med_array=None,lower_array=None,upper_array=None):
        self.min_phase = min_phase

        self.max_phase = max_phase

        self.center_phase = (min_phase + max_phase)/2

        self.med_array = med_array

        self.lower_array = lower_array

        self.upper_array = upper_array

    def collapse_bin(self):
        collapsed_array = add_asym(self.med_array,self.lower_array,self.upper_array)

        #print(self.med_array[0])

        collapsed_array = collapsed_array/len(self.med_array)

        return collapsed_array


def make_phase_bins(bin_size):
    '''returns an array of Phase_Space_Bins objects equally spaced spanning
    the phase space in even intervals of bin_size'''

    bins = []

    running_phase = 0

    while running_phase <= 1:
        min_phase = running_phase

        running_phase += bin_size

        max_phase = running_phase

        new_med = []

        new_up = []

        new_low = []

        bins.append(Phase_Space_Bins(min_phase,max_phase,new_med,new_up,new_low))

    '''
    arr = [bins[0].med_array is binn.med_array for binn in bins]
    print(arr)
    '''


    return bins


def populate_bins(HR_list,bin_size = 0.01):
    bins = make_phase_bins(bin_size)

    for object in HR_list:
        for i in range(len(object.phase)):

            med = object.meds[i]

            #convert them to their distance, rather than position
            lower = object.lowers[i] + med

            upper = object.uppers[i] - med

            phase = object.phase[i]

            matched = False
            for j,bin_ in enumerate(bins):
                if phase > bin_.min_phase and phase < bin_.max_phase:
                    bin_.med_array.append(med)
                    bin_.lower_array.append(lower)
                    bin_.upper_array.append(upper)

                    matched = True

                    #print(f'Med: {med} put into bin {j+1} with center phase {bin.center_phase}')

                    break

            assert matched


    '''
    compare_arr = [bins[2].med_array[i] == bins[3].med_array[i] for i in range(min(len(bins[0].med_array),len(bins[1].med_array)))]

    assert bins[0] is not bins[1]

    assert bins[0].med_array is not bins[1].med_array

    assert False in compare_arr
    '''

    for i,bin in enumerate(bins):
        if len(bin.med_array) == 0:
            del bins[i]

    return bins


def make_plot(bins):

    phases = []
    HRs = []
    lows = []
    highs = []

    for bin in bins:
        phase = bin.center_phase
        HR,low,high = bin.collapse_bin()

        phases.append(phase)
        HRs.append(HR)
        lows.append(low)
        highs.append(high)

    plt.errorbar(phases,HRs,yerr=[lows,highs],fmt='ko')

    plt.show()




if __name__ == '__main__':
    coords = sys.argv[1]
    N = int(sys.argv[2])
    obsids = sys.argv[3].split(',')
    period = float(sys.argv[4])
    t0 = float(sys.argv[5])

    dict = {i:None for i in obsids}

    #first run the hardness ratio calculation for each obsid
    #and read back in the files created
    for obsid in obsids:
        text = main(obsid,coords,N,None,subtract_start=False)

        HR = HR_graph(text)

        HR.make_phase(period,t0)

        dict[obsid] = HR

    #now its time to bin and plot
    bins = populate_bins(dict.values())

    make_plot(bins)
