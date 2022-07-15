import numpy as np
from scipy import stats
import scipy as sc
import sys
from make_full_plot_for_ttest import make_full_plot
import glob
import os
from statsmodels.stats.rates import test_poisson_2indep
import multiprocessing as mp

#represents all the observations of a single source
#and some statistics useful to the program
class Source_All:
    def __init__(self,obs):
        #obs is an array of Source objects which represent all the
        #observations of this source
        self.obs = obs

        #an array of source objects which represent all non-zero observations
        #of the source
        self.obs_nonzero = [i for i in obs if not i.is_zero()]

        #the average count rate across all observations of the source
        #units of ks^-1
        self.average_count_rate = 1000*sum([i.total_counts for i in obs])/sum([i.duration for i in obs])

        #the average count rate across all observations of the source with
        #non-zero counts
        #units of ks^-1
        self.av_count_rate_nozeros = 1000*sum([i.total_counts for i in obs_nonzero])/sum([i.duration for i in obs_nonzero])


#This object represents the abstract type of a source-obsid pair
#ie, this object represents all the information about a single obersvation
#of a single source.
class Source:
    def __init__(self, lightcurve=None,obsid=None, position=None,classification=None,
                t_interest=[]):
        #the obsid
        self.obsid = obsid

        #the position of the source is sexigesimal format
        self.position = position

        #lightcurve is a 2D array of the raw lightcurve data read in directly
        #from the light curve files
        #intended to be created with a line such as
        #np.loadtxt(file, skiprows = 1)
        self.lightcurve = lightcurve

        #arrays based on the lightcurve
        self.counts = self.lightcurve[::,4]
        self.times = self.lightcurve[::,2]

        #total counts
        self.total_counts = sum(self.counts)

        #start_time is the MJD of the start of the observation
        self.start_time  = self.times[0]

        #end_time is the MJD of the end of the observation
        self.end_time = self.times[-1]

        #duration of observation
        self.duration = self.end_time - self.start_time

        #counts per kilosecond averaged over the observation
        self.count_rate = 1000*self.total_counts/self.duration

        #some method of classifying the interest level in the source
        #TBD
        self.classification = classification

        #all times throughout the observation where something interesting is
        #found to be happening
        self.t_interest = t_interest


    #returns an unbinned (binned at frame-time) cumulative counts array for the
    #source
    def cumo_counts(self):
        return np.cumsum(self.counts)

    #returns true if the rate per ks is above "rate" false otherwise
    def count_cut(self,rate):
        return self.count_rate >= rate

    #returns true if there are zero counts, false otherwise
    def is_zero(self):
        if self.total_counts == 0:
            return True
        else:
            return False

    #returns the binned counts array binned at binsize
    def make_binned_counts(self,binsize):
        if self.is_zero():
            return None


        #trim the ends off, if they have 0 counts
        i = 0
        while self.counts[i] == 0:
            i += 1
        i -= 1

        j = -1
        while self.counts[j] == 0:
            j -= 1
        j +=1

        time = self.times[i:j]
        counts = self.counts[i:j]

        time = self.times
        counts = self.counts

        try:
            bin_length = int(binsize/3.24014)
            bin_edges = [time[i] for i in range(len(time)) if i%bin_length == 0]

            binned_counts = stats.binned_statistic(time,counts,statistic='sum',bins=bin_edges)[0]
        except Exception as d:
            raise d

            try:
                bin_edges.append(time[-1])
                binned_counts = stats.binned_statistic(time,counts,statistic='sum',bins=bin_edges)[0]

            except Exception as e:
                raise e



        return binned_counts

    #true if the low count flare condition is met
    #specifically, if there is at least one 1ks bin which is {threshold} sigma
    #above the baseline up to that bin (forward or backwards) with greater than
    #{min_counts} counts in that bin
    def find_low_flare(self,threshold=3,min_counts=7):
        #intended only to be run on low count sources. Will return false if the
        #count rate is less than 4 ks^-1
        if self.total_counts < min_counts or self.count_cut(4):
            return False

        binned_list = self.make_binned_counts(1000)

        for i in range(len(binned_list)):
            if binned_list[i]>min_counts:
                z_forward = np.abs(stats.zscore(binned_list[0:i]))
                outliers = np.where(z_forward>threshold)
                if len(outliers[0]) > 0:
                    return True

            if binned_list[-1*i] > min_counts:
                z_back = np.abs(stats.zscore(binned_list[-i:-1]))
                outliers = np.where(z_back>threshold)
                if len(outliers[0]) > 0:
                    return True

        return False

    def e_test(self,pthresh=2.867E-7,binsize=1000, debug = False):
        #default pthresh set to 5 sigma detection level

        if self.is_zero() or self.total_counts < 3:
            return None

        #first we need to chunk up the data
        #AKA dividing the data into chunks where each chunk is delimitted
        #by the light curve passing through the median of the dataset
        binned_counts = self.make_binned_counts(binsize)
        if binned_counts is None:
            return None

        #median = sc.ndimage.median(binned_counts)

        #testing unsing mean instead
        median = np.mean(binned_counts)

        chunk_edges = []

        for i, bin in enumerate(binned_counts):
            try:
                next_bin = binned_counts[i+1]
            except IndexError:
                break

            decrease = bin >= median and next_bin <= median and bin != next_bin
            increase = bin <= median and next_bin >= median and bin != next_bin

            if decrease or increase:
                chunk_edges.append(i+1)

        #turning the chunk edges into the chunks themselves
        chunks = [[binned_counts[i] for i in range(chunk_edges[j],chunk_edges[j+1])] for j in range(len(chunk_edges)) if j != len(chunk_edges)-1]

        #the above line misses the last chunk, we will add it in manually
        chunks.append([i for i in binned_counts[chunk_edges[-1]:]])


        #returns the most significantly different chunk from among all chunks
        #as determined by the e-test
        def find_interesting_chunk(self,chunks,debug):
            lowest_p = 1.0
            most_interesting = None

            #Loop through all the chunks, for each testing if the mean inside the chunk
            #is significantly different than the mean outside of it
            for i,chunk in enumerate(chunks):
                if len(chunk) == 1:
                    continue

                #build all_others, an array which contains all the chunks beside the
                #one being examined
                if i == 0:
                    all_others = chunks[1:]
                else:
                    all_others = chunks[:i]
                    seg2 = chunks[i+1:]

                    all_others.extend(seg2)

                #flatten all_others
                all_others = [bin for sub in all_others for bin in sub]

                #if there is a flare, we don't want to consider it if it has fewer than 10 counts
                if np.mean(chunk) > np.mean(all_others) and sum(chunk) < 10:
                    continue
                #if there's a dip and there are fewer than 30 total counts, we also don't want it
                elif np.mean(chunk) < np.mean(all_others) and self.total_counts < 30:
                    continue

                #TODO: Find better fix for one bin chunks than just skipping them

                #perform the e-test
                t_stat,p_val = test_poisson_2indep(sum(chunk),len(chunk),sum(all_others),len(all_others),method='etest')

                #stat, p_val = stats.mannwhitneyu(chunk,all_others)

                if p_val < pthresh:

                    if debug:
                        print(f'''Detection in chunk {i}
    This chunk: {chunk}
    all_others: {all_others}
    p_val: {p_val}
    ''')

                    if p_val < lowest_p:
                        most_interesting = i
                        lowest_p = p_val

                        if debug:
                            print('New most interesting chunk')

                else:
                    if debug:
                        print(f'''No detection in chunk {i}
    This chunk: {chunk}
    all_others: {all_others}
    p_val: {p_val}
    ''')

            return most_interesting


        interesting_chunks = []

        repeat = True
        while repeat:
            new_interesting = find_interesting_chunk(self,chunks,debug)

            if new_interesting is not None:
                interesting_chunks.append(new_interesting)

            else:
                repeat = False


        if not interesting_chunks:
            return None
        else:
            #print(interesting_chunks)

            out = [[chunk_edges[i],len(binned_counts)-1] if i == len(chunk_edges) - 1 else [chunk_edges[i],chunk_edges[i+1]]for i in interesting_chunks]
            #print(out)

            #convert out to time
            #to do so, we need a trimmed time array
            #just miimic what is done when the binned times are made
            i = 0
            while self.counts[i] == 0:
                i += 1
            i -= 1

            j = -1
            while self.counts[j] == 0:
                j -= 1
            j +=1

            trimmed_time = self.times[i:j]
            trimmed_time = [i - self.times[0] for i in trimmed_time]
            bin_length = int(binsize/3.24014)
            trimmed_times_binned = [trimmed_time[i] for i in range(len(trimmed_time)) if i%bin_length == 0]

            if debug:
                print(out)
                print(len(trimmed_times_binned))
                print(trimmed_times_binned)

            detections = [trimmed_times_binned[index] for sub in out for index in sub]

            self.t_interest.extend(detections)

        if debug:
            print(detections)
            print(out)

        return detections

def analyze_file_alt(file):
    if __name__ == '__main__':
        #verbose is an integer, controls the amount of info printed to the screen
        #0,1,2,3
        #0 prints no information
        #1 prints only the number of detections
        #2 prints every time a detection is made
        #3 prints all debug information
        #verbose = int(sys.argv[1])
        verbose = int(sys.argv[2])

        #binsize in seconds
        #I default to 500
        if len(sys.argv) < 4:
            binsize = 500
        else:
            binsize = int(sys.argv[3])

        if verbose > 2:
            debug = True
        else:
            debug = False
    else:
        debug = False
        binsize = 500
        verbose = 1


    if 'ascii' in file:
        skip = 0
    else:
        skip = 1

    source = Source(lightcurve = np.loadtxt(file, skiprows = skip))

    try:
        detections = source.e_test(binsize = binsize,debug=debug)
    except Exception as e:
        print(f'Error in file: {file}\n')
        raise e

    if detections is not None:
        if verbose > 1:
            print(f'Detection in {file}, making plot...')

        try:
            make_full_plot(file,f'./test/detections',lines=detections)
        except Exception as e:
            print(f'Error plotting file {file} with lines {detections}')
            print(f'Raising error...')
            raise e

        return 1

    else:
        return 0

def analyze_file(file):
    if __name__ == '__main__':
        #verbose is an integer, controls the amount of info printed to the screen
        #0,1,2,3
        #0 prints no information
        #1 prints only the number of detections
        #2 prints every time a detection is made
        #3 prints all debug information
        verbose = int(sys.argv[1])

        #binsize in seconds
        #I default to 500
        if len(sys.argv) < 3:
            binsize = 500
        else:
            binsize = int(sys.argv[2])

        if verbose > 2:
            debug = True
        else:
            debug = False
    else:
        debug = False
        binsize = 500
        verbose = 1


    if 'ascii' in file:
        skip = 0
    else:
        skip = 1

    try:
        source = Source(lightcurve = np.loadtxt(file, skiprows = skip))
    except StopIteration as err:
        return 0

    try:
        detections = source.e_test(binsize = binsize,debug=debug)
    except Exception as e:
        print(f'Error in file: {file}\n')
        raise e

    if detections is not None:
        if verbose > 1:
            print(f'Detection in {file}, making plot...')

        try:
            make_full_plot(file,f'{galaxy}/detections',lines=detections)
        except Exception as e:
            print(f'Error plotting file {file} with lines {detections}')
            print(f'Raising error...')
            raise e

        return 1

    else:
        return 0


"""
if __name__ == '__main__':
    #galaxies = glob.glob('./all_sk_galaxies/*')

    galaxies = ['M82','M83','M84','NGC1313',"NGC2403",'NGC247','NGC300',
                'NGC4038','NGC4631','NGC4736','NGC5907','NGC6946','NGC1052',
                'NGC2782','NGC2903']

    galaxies = [f'./all_sk_galaxies/{i}' for i in galaxies]

    #verbose is an integer, controls the amount of info printed to the screen
    #0,1,2,3
    #0 prints no information
    #1 prints only the number of detections
    #2 prints every time a detection is made
    #3 prints all debug information
    verbose = int(sys.argv[1])

    #binsize in seconds
    #I default to 500
    if len(sys.argv) < 3:
        binsize = 500
    else:
        binsize = int(sys.argv[2])

    if verbose > 2:
        debug = True
    else:
        debug = False

    for galaxy in galaxies:
        try:
            os.makedirs(f'{galaxy}/detections')
        except FileExistsError:
            pass

        if verbose > 0:
            print(f'Analyzing {galaxy}')

        files = glob.iglob(f'{galaxy}/textfiles/*.txt')

        '''
        p = mp.Pool(int(mp.cpu_count()/2))
        n_detections = p.map(analyze_file,files)
        '''

        n_detections = []
        for file in files:
            n_detections.append(analyze_file(file))

        if verbose > 0:
            print(f'Detections made in {sum(n_detections)} files')
"""

if __name__ == '__main__':
    file = sys.argv[1]

    analyze_file_alt(file)

#"""
