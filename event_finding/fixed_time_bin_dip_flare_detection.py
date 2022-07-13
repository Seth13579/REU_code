
'''
S. Larner

This file contains a method for detecting dips and flares in text file light
curves.

The script converts the light curve into a cumulative count plot. It then
tests for deviations in the slope of this plot.

An increase in slope corresponds to a flare, a decrease corresponds to a dip.

The script bins the cumulative counts in constant time bins. Then for each
bin it calculates the slope a specified distance in time forward. This time is
compared to the slope from the last time the slope changed to the start of the
current bin.

The script has three controllable parameters:
    binning_size controls the size of the time bins (in seconds).

    num_forward controls the number of bins forward the slope at a given point
    is calculated from. For example, if binning_size = 100 and bins_forward = 3
    the slope at a given point is the slope over the next 300 seconds.

    threshold is the percent difference between the old slope and the new which
    indicates a change in slope.
'''

import numpy as np
import sys

'''
plot_cumo read in a Sammarth style text file to make a cumulative counts plot.

If show, the function will produce a cumulative counts plot and display it to
the screen.

If save, the function will produce a cumulative counts plot and save it to disk.

Lines is an array-like with time positions of vertical lines to put on the plot.
'''
def plot_cumo(file,show=False,save=False,lines=None):

    #read the file
    data = np.loadtxt(file, skiprows = 1)
    counts = data[::,4]

    time = data[::,2]
    bin_num = data[::,0]

    #calculate the cumulative counts
    cumo = [sum(counts[0:i]) for i in range(len(data))]

    trunc_name = file.split('/')[-1].split('.')[:-1]
    trunc_name = '.'.join(trunc_name)

    #adjust time
    time = [(i - time[0]) for i in time]

    if show or save:
        import matplotlib.pyplot as plt

        plt.plot(time,cumo)
        plt.xlabel('Time since obs start (s)')
        plt.ylabel('Cumulative Counts')
        plt.title(f'Cumulative counts for {trunc_name}')

        if lines is not None:
            plt.vlines(lines,ymin=0,ymax=max(cumo),colors = 'r',linestyles='dotted')

        if show:
            plt.show()
        if save:
            plt.savefig(f'{trunc_name}_cumo.pdf')

        plt.close()

    return

'''
slope is an object representation of the slope between two points on the
cumulative counts plot
'''
class slope:
    def __init__(self, current_time, bins_forward, value = 0, error = False, different = False):
        #current_time is the time at which the slope is caulcated from
        self.current_time = current_time

        #bins forward is the number of time bins that the slope is over
        #relative to current_time
        #if negative, the slope is calculated backwards from current_time
        self.bins_forward = bins_forward

        #value is the value of the slope
        #In units of counts/kilosecond
        self.value = value

        #error tracks if an error occured during the creation of slope
        #will be used to ensure that calculations are not perfored on bad slopes
        self.error = error

        #different tracks if this slope is different from the reference slope
        self.different = different

'''
sets the value in a slope object
'''
def set_slope_value(slope_obj, cumo_arr, time_arr, i, binning_size):
    try:
        slope_obj.value = (cumo_arr[i + slope_obj.bins_forward*binning_size] - cumo_arr[i])
        slope_obj.value /= (time_arr[i + slope_obj.bins_forward*binning_size] - time_arr[i])

    except IndexError:
        slope_obj.error = True

    return


'''
compares the slope between a test slope and a reference slope and sets
the value of different in the test slope to true or false as specified.
'''
def slope_compare(test_slope, ref_slope, thresh):
    if ref_slope.value == test_slope.value:
        test_slope.different = False

    #We should prevent the script from recognizing the first start up as a change,
    #but needs to be revisited.
    #Will probably need some way to trim the light curve to not consider the
    #very start when there is no counts yet.
    elif ref_slope.value == 0 and ref_slope.current_time == 0:
        test_slope.different = False

    elif not test_slope.error:
        percent_dif = abs(ref_slope.value - test_slope.value)/((ref_slope.value+test_slope.value)/2)
        if percent_dif >= thresh:
            test_slope.different = True
        else:
             test_slope.different = False
    else:
         test_slope.different = False

    return

'''
test_slope examines a light curve file and returns a list of all points
with a slope that deviates by more than threshold.

The function can bin the array on the fly by setting binning_size greater than 1.
binning_size must be an integer greater than zero or auto.
If 'auto' is used, the program will use bins of a size 5/average count rate.

num forward is the number of bins forward the testing slope is created from

Threshold must be a positive number.
'''
def test_slope(file, binning_size, num_forward, threshold):
    #to calulate the current slope, we need to know the last time the slope changed
    #to start, we will just use zero
    #this is an array index, not a time
    last_change = 0

    #read in from the file
    data = np.loadtxt(file, skiprows = 1)
    counts_arr = data[::,4]
    if sum(counts_arr)==0:
        return []

    time_arr = data[::,2]
    time_arr = [(i - time_arr[0]) for i in time_arr]

    if binning_size == 'auto':
        counts = sum(counts_arr)
        time = time_arr[-1]
        av_count_rate = counts/time

        #make it so on average, 5 counts per bin
        binning_size = int(5/av_count_rate)
        binning_size = max(500,binning_size)

    #call plot_cumo to create the cumo array
    cumo_arr = cumo = [sum(counts_arr[0:i]) for i in range(len(counts_arr))]

    #trim the arrays to exclude start up and endings:
    i = 0
    while counts_arr[i] == 0:
        i += 1

    j = -1
    while counts_arr[j] == 0:
        j -= 1

    time_arr = time_arr[i:j]
    counts_arr = counts_arr[i:j]
    cumo_arr = cumo_arr[i:j]

    #change_positions is the list of positions in the array where a significant
    change_positions = []

    for i in range(0,len(cumo_arr),int(binning_size/3.2)):
        if i == 0:
            continue

        #calculate and test the slope the specified number of time bins forward
        test_slope = slope(time_arr[i],num_forward)
        set_slope_value(test_slope, cumo_arr, time_arr, i, binning_size)

        #reference_slope is the slope to which we are comparing
        reference_slope = slope(time_arr[last_change], last_change - i)
        reference_slope.value = cumo_arr[last_change] - cumo_arr[i]
        reference_slope.value /= (time_arr[last_change] - time_arr[i])


        slope_compare(test_slope,reference_slope,threshold)
        if test_slope.different:
            change_positions.append(i)
            last_change = i

    #convert the list of positions in the time array to a list of times
    change_times = [time_arr[i] for i in change_positions]

    return change_times


def main():
    #a file to test on, can be changed as needed
    file = 'test_lc.txt'

    #these are the values I have found work best on this file
    change_times = test_slope(file, 250, 1, .7)

    #plot_cumo(file,show=True,save=False,lines=change_times)

    return

if __name__ == '__main__':
    main()
