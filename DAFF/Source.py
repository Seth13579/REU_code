import numpy as np
import pickle
from re import sub,search
import re
import os
import matplotlib.pyplot as plt
import glob
import sys
import gc
from astropy.stats import bayesian_blocks
import warnings
from scipy import stats

#for the E-test
from statsmodels.stats.rates import test_poisson_2indep

#for the alternative E-test
#from poisson_etest import poisson_etest

#HR dependencies
from extract_counts import *
from BEHR_countbins import *
from run_BEHR import *
from create_behr_files import region_area

def poisson_twosample(count1, exposure1, count2, exposure2, ratio_null=1,method='score', alternative='two-sided'):
    '''test for ratio of two sample Poisson intensities

    If the two Poisson rates are g1 and g2, then the Null hypothesis is

    H0: g1 / g2 = ratio_null

    against one of the following alternatives

    H1_2-sided: g1 / g2 != ratio_null
    H1_larger: g1 / g2 > ratio_null
    H1_smaller: g1 / g2 < ratio_null

    Parameters
    ----------
    count1: int
        Number of events in first sample
    exposure1: float
        Total exposure (time * subjects) in first sample
    count2: int
        Number of events in first sample
    exposure2: float
        Total exposure (time * subjects) in first sample
    ratio: float
        ratio of the two Poisson rates under the Null hypothesis. Default is 1.
    method: string
        Method for the test statistic and the p-value. Defaults to `'score'`.
        Current Methods are based on Gu et. al 2008
        Implemented are 'wald', 'score' and 'sqrt' based asymptotic normal
        distribution, and the exact conditional test 'exact-cond', and its mid-point
        version 'cond-midp', see Notes
    alternative : string
        The alternative hypothesis, H1, has to be one of the following

           'two-sided': H1: ratio of rates is not equal to ratio_null (default)
           'larger' :   H1: ratio of rates is larger than ratio_null
           'smaller' :  H1: ratio of rates is smaller than ratio_null

    Returns
    -------
    stat, pvalue two-sided

    not yet
    #results : Results instance
    #    The resulting test statistics and p-values are available as attributes.


    Notes
    -----
    'wald': method W1A, wald test, variance based on separate estimates
    'score': method W2A, score test, variance based on estimate under Null
    'wald-log': W3A
    'score-log' W4A
    'sqrt': W5A, based on variance stabilizing square root transformation
    'exact-cond': exact conditional test based on binomial distribution
    'cond-midp': midpoint-pvalue of exact conditional test

    The latter two are only verified for one-sided example.

    References
    ----------
    Gu, Ng, Tang, Schucany 2008: Testing the Ratio of Two Poisson Rates,
    Biometrical Journal 50 (2008) 2, 2008

    '''

    # copied from statsmodels.stats.weightstats
    def _zstat_generic2(value, std_diff, alternative):
        '''generic (normal) z-test to save typing

        can be used as ztest based on summary statistics
        '''
        zstat = value / std_diff
        if alternative in ['two-sided', '2-sided', '2s']:
            pvalue = stats.norm.sf(np.abs(zstat))*2
        elif alternative in ['larger', 'l']:
            pvalue = stats.norm.sf(zstat)
        elif alternative in ['smaller', 's']:
            pvalue = stats.norm.cdf(zstat)
        else:
            raise ValueError('invalid alternative')
        return zstat, pvalue

    # shortcut names
    y1, n1, y2, n2 = count1, exposure1, count2, exposure2
    d = n2 / n1
    r = ratio_null
    r_d = r / d

    if method in ['score']:
        stat = (y1 - y2 * r_d) / np.sqrt((y1 + y2) * r_d)
        dist = 'normal'
    elif method in ['wald']:
        stat = (y1 - y2 * r_d) / np.sqrt(y1 + y2 * r_d**2)
        dist = 'normal'
    elif method in ['sqrt']:
        stat = 2 * (np.sqrt(y1 + 3 / 8.) - np.sqrt((y2 + 3 / 8.) * r_d))
        stat /= np.sqrt(1 + r_d)
        dist = 'normal'
    elif method in ['exact-cond', 'cond-midp']:
        from statsmodels.stats import proportion
        bp = r_d / (1 + r_d)
        y_total = y1 + y2
        stat = None
        pvalue = proportion.binom_test(y1, y_total, prop=bp, alternative=alternative)
        if method in ['cond-midp']:
            # not inplace in case we still want binom pvalue
            pvalue = pvalue - 0.5 * stats.binom.pmf(y1, y_total, bp)

        dist = 'binomial'

    if dist == 'normal':
        return _zstat_generic2(stat, 1, alternative)
    else:
        return stat, pvalue

def load_src_all(file):
    with open(file,'rb') as f:
        return pickle.load(f)

def unglob(arr,force=False):

    if len(arr) > 1 and not force:
        print('Multiple files found when 1 was expected:')
        print(arr)
        cont = None
        while cont != 'c' and cont != 'a':
            cont = input('(c)ontinue or (a)bort? ')
        if cont == 'a':
            raise Exception

    elif len(arr) == 0:
        #print('No files found')
        raise Exception

    return str(arr[0]).strip("'[]'")

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

def trimzeros(array):
    #trim the ends off, if they have 0 counts
    i = 0
    while array[i] == 0:
        i += 1

    if array[-1] == 0:
        #print('Trimming off the end')
        j = -1
        while array[j] == 0:
            j -= 1
        j +=1
    else:
        j = None

    return i,j

#represents all the observations of a single source
#and some statistics useful to the program
class Source_All:
    def __init__(self,obs,ra=None,dec=None):
        #the right ascension, in sexigesimal
        self.ra = ra

        #the declination, in sexigesimal
        self.dec = dec

        #obs is an array of Source objects which represent all the
        #observations of this source
        self.obs = obs

        _ra = sub('\:','',self.ra)
        _dec = sub('\:','',self.dec)
        self.name = f'J{_ra}{_dec}'

        #the average count rate across all observations of the source
        #units of ks^-1
        self.average_count_rate = 1000*sum([i.total_counts for i in self.obs])/sum([i.duration for i in self.obs])

    def save(self,outdir='.'):
        outfile = f'{outdir}/SOURCE_ALL_{self.name}.pkl'

        with open(outfile,'wb') as outp:
            pickle.dump(self, outp,protocol=pickle.HIGHEST_PROTOCOL)

        return

    def rename(self):
        '''re initializes the self.name field to match the current self.ra and self.dec'''

        _ra = sub('\:','',self.ra)
        _dec = sub('\:','',self.dec)
        self.name = f'J{_ra}{_dec}'

        return

    def optimal_rename(self):
        '''renames itself to the coordinates of the region within with the lowest error'''
        #cannot be called if the region field of all source objects in self.obs
        #do not have their region specified

        best_name = sorted(self.obs, key = lambda x: x.get_a())[0]

        best_pos = best_name.position

        if '+' in best_pos:
            best_ra = best_pos.split('+')[0]
            best_dec = '+' + best_pos.split('+')[1]

        elif '-' in best_pos:
            best_ra = best_pos.split('-')[0]
            best_dec = '-' + best_pos.split('-')[1]

        else:
            raise Exception

        self.ra = best_ra
        self.dec = best_dec

        self.rename()

        #also have to update the parent names in each of the child objects
        for src in self.obs:
            src.parent_name = f'{sub(":","",self.ra)}{sub(":","",self.dec)}_{src.obsid}'

        return

    def save_lcs(self,outdir='.',clobber = False):
        for source in self.obs:
            lc = source.lightcurve

            filename = f'{outdir}/{self.name}_{source.obsid}_lc.fits.txt'

            if not clobber:
                try:
                    assert not os.path.exists(filename)
                except AssertionError as e:
                    print(f'Cannot save textfile {filename}')
                    print(f'File with that name already exists')
                    raise e

            header = 'TIME_BIN TIME_MIN TIME TIME_MAX COUNTS STAT_ERR AREA EXPOSURE COUNT_RATE COUNT_RATE_ERR'
            np.savetxt(filename, lc,header=header,fmt='%s',comments='')

#This object represents the abstract type of a source-obsid pair
#ie, this object represents all the information about a single obersvation
#of a single source.
class Source:
    def __init__(self, lightcurve=None,obsid=None, position=None,classification=False,
                t_interest=None, region=None, region_object=None, parent=None):
        #the obsid
        self.obsid = obsid


        #the position of the source is sexigesimal format
        self.position = position

        name = f'{sub(":","",self.position)}_{self.obsid}'
        self.name = name

        if '+' in self.position:
            self.ra = self.position.split('+')[0]
            self.dec = '+' + self.position.split('+')[1]
        else:
            self.ra = self.position.split('-')[0]
            self.dec = '+' + self.position.split('-')[1]

        self.position_basic = sub('\:','',position)

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

        if self.total_counts != 0:

            #trim the ends off, if they have 0 counts
            i = 0
            while self.counts[i] == 0:
                i += 1

            if self.counts[-1] == 0:
                #print('Trimming off the end')
                j = -1
                while self.counts[j] == 0:
                    j -= 1
                j +=1
            else:
                j = None

            if j is not None:
                #print(f'Trimming: [{i}:{j}]')
                self.times = self.times[i:j]
                self.counts = self.counts[i:j]
            else:
                self.times = self.times[i:]
                self.counts = self.counts[i:]

        #start_time is the MJD of the start of the observation
        self.start_time  = self.times[0]

        #end_time is the MJD of the end of the observation
        self.end_time = self.times[-1]

        #duration of observation
        self.duration = self.end_time - self.start_time

        #counts per kilosecond averaged over the observation
        self.count_rate = 1000*self.total_counts/self.duration

        #True if the source is interesting, False otherwise
        self.classification = classification

        #all times throughout the observation where something interesting is
        #found to be happening
        self.t_interest = t_interest

        #optional region text describing the source
        self.region = region

        #the region object which created this source
        self.region_object = region_object

        #HR information
        self.HR = None

        #parent is the parent Source_All object to which this source belongs
        #even if unset, will auto update to the correct name when
        #the optimal_rename routine is run in the parent
        self.parent = parent

        if self.parent is not None:
            self.parent_name = f'{sub(":","",self.parent.ra)}{sub(":","",self.parent.dec)}_{self.obsid}'
        else:
            self.parent_name = None


    #returns the semimajor axis of the region which made this source
    def get_a(self):
        '''Returns the semimajor axis of the region associated with this source'''
        if self.region == None:
            return None

        else:
            line = self.region
            a = max(line[7:].strip('()').split(',')[2],line[7:].strip('()').split(',')[3]).rstrip('"')

            return a

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

    #returns the time and counts arrays binned at binsize
    def make_binned_counts(self,binsize):
        if self.is_zero():
            return None


        time = self.times - self.start_time
        counts = self.counts

        binned_times,binned_counts = bin_array(time,counts,binsize)

        #lop off the last bin to prevent the issue with thinking
        #it has drastically fewer counts just because it is much shorter
        #not the best solution, but the easiest to implement in this
        #counts-per-bin backed paradigm
        return binned_times[:-1],binned_counts[:-1]

    #makes a four panel light curve for the source
    #one cumulative plot
    #three light curves, defined by the bin sizes in binsizes
    def make_fourpanel_plot(self,outdir,binsizes=[500,1000,2000],save=True,
                            show=False,lines=None, highlight=None):
        if lines is None:
            lines = self.t_interest

        if self.parent_name is not None:
            name = f'{self.parent_name}_{self.obsid}'
        else:
            name = f'{self.name}_{self.obsid}'

        try:
            binned_times_1,binned_counts_1 = self.make_binned_counts(binsizes[0])
        except TypeError as e:
            print(self.name)
            raise e
        binned_times_2,binned_counts_2 = self.make_binned_counts(binsizes[1])
        binned_times_3,binned_counts_3 = self.make_binned_counts(binsizes[2])

        trim_low,trim_high = trimzeros(self.counts)
        trim_counts = self.counts[trim_low:trim_high]
        trim_time = self.times[trim_low:trim_high]

        cumo = np.cumsum(trim_counts)

        #adding the num and clear pars to try and prevent memory leak
        fig, axs = plt.subplots(2,2,num=1,clear=True)
        cumo_plt = axs[0,0]
        bin1_plt = axs[0,1]
        bin2_plt = axs[1,0]
        bin3_plt = axs[1,1]

        cumo_plt.plot(trim_time-trim_time[0],cumo)
        cumo_plt.set_title('Cumulative Counts')

        bin1_plt.step(binned_times_1,binned_counts_1,where='mid')
        bin1_plt.set_title(f'Bin sizes: {round(binsizes[0])}')

        bin2_plt.step(binned_times_2,binned_counts_2,where='mid')
        bin2_plt.set_title(f'Bin sizes: {round(binsizes[1])}')

        bin3_plt.step(binned_times_3,binned_counts_3,where='mid')
        bin3_plt.set_title(f'Bin sizes: {round(binsizes[2])}')

        cumo_plt.set(ylabel='Counts')
        bin2_plt.set(xlabel='Time (s)',ylabel='Counts')
        bin3_plt.set(xlabel='Time (s)')


        fig.suptitle(name)

        if lines is not None:
            for ax in axs.flat:
                ax.vlines(lines,ymin=ax.get_ylim()[0],ymax=ax.get_ylim()[1],colors='r',linestyles='dotted')

        if highlight is not None:
            for ax in axs.flat:
                lim = ax.get_ylim()
                for i in np.arange(0,len(highlight),2):
                    ax.fill_betweenx(y=lim,x1=highlight[i],x2=highlight[i+1],alpha=.3,color='green')

        plt.subplots_adjust(hspace=.3)

        if save:
            fig.savefig(f'{outdir}/{name}_full_lc.png')
            fig.savefig(f'{outdir}/{name}_full_lc.pdf')

        if show:
            plt.show()


        plt.close()

        return

    #runs the chunked e-test process
    #returns an array of time edges of interesting chunks
    #interesting is chunks with p-vals below pthresh
    def e_test(self,pthresh=2.867E-7,binsize=1000, debug = False):
        #default pthresh set to 5 sigma detection level

        if self.is_zero() or self.total_counts < 3:
            return None

        #first we need to chunk up the data
        #AKA dividing the data into chunks where each chunk is delimitted
        #by the light curve passing through the median of the dataset
        binned_time,binned_counts = self.make_binned_counts(binsize)
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

        #initialize the list of interesting chunks which pass the significance test
        interesting_chunks = []

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
            elif np.mean(chunk) < np.mean(all_others) and self.total_counts < 30:
                continue

            #perform the t-test
            #TODO: Find better fix for one bin chunks than just skipping them

            t_stat,p_val = test_poisson_2indep(sum(chunk),len(chunk),sum(all_others),len(all_others),method='etest')

            #stat, p_val = stats.mannwhitneyu(chunk,all_others)

            if p_val < pthresh:
                if not self.classification:
                    self.classification = True

                interesting_chunks.append(i)

                if debug:
                    print(f'''Detection in chunk {i}
This chunk: {chunk}
all_others: {all_others}
p_val: {p_val}
''')

            else:
                if debug:
                    print(f'''No detection in chunk {i}
This chunk: {chunk}
all_others: {all_others}
p_val: {p_val}
''')


        if not interesting_chunks:
            return None
        else:
            #print(interesting_chunks)

            out = [[chunk_edges[i],len(binned_counts)-1] if i == len(chunk_edges) - 1 else [chunk_edges[i],chunk_edges[i+1]]for i in interesting_chunks]


            try:
                detections = [binned_time[index] for sub in out for index in sub]
            except IndexError as e:
                print('\n\n************\n\n')
                print('Index error creating detections array')

                #print(out)
                #print(len(trimmed_times_binned))
                raise e

            if self.t_interest is None:
                self.t_interest = detections
            else:
                self.t_interest.extend(detections)

        if debug:
            print(f'Detections: {detections}')
            print(f'out:{out}')
            print(f'Binned counts: {binned_counts}')
            print(f'Binned times: {binned_time}')

        return detections

    #runs the chunked e-test process
    #returns an array of time edges of interesting chunks
    #interesting is chunks with p-vals below pthresh
    #this version encorporates the repeating when an interesting chunk is found
    def e_test_repeat(self,pthresh=2.867E-7,binsize=1000, debug = False):
        #default pthresh set to 5 sigma detection level

        if self.total_counts < 3:
            return None

        #first we need to chunk up the data
        #AKA dividing the data into chunks where each chunk is delimitted
        #by the light curve passing through the median of the dataset
        binned_time,binned_counts = self.make_binned_counts(binsize)

        if binned_counts is None:
            return None

        def make_chunks(binned_counts):
            median = np.mean(binned_counts)

            #chunk edges are the indeces in binned_counts where chunks are delimitted
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
            #chunks are an array of sub arrays,
            #each sub array is a list of the values of binned_counts in the chunk
            chunks = [[binned_counts[i] for i in range(chunk_edges[j],chunk_edges[j+1])] for j in range(len(chunk_edges)) if j != len(chunk_edges)-1]

            #the above line misses the last chunk, we will add it in manually
            chunks.append([i for i in binned_counts[chunk_edges[-1]:]])

            return (chunk_edges,chunks)

        chunk_edges,chunks = make_chunks(binned_counts)

        def find_best_chunk(chunks):
            #initialize the dict of interesting chunks which pass the significance test
            #key will be the index of the chunk
            #value will be the p_value
            interesting_chunks = {}

            #Loop through all the chunks, for each testing if the mean inside the chunk
            #is significantly different than the mean outside of it
            #interesting chunks get added to interesting_chunks
            for i,chunk in enumerate(chunks):
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

                if len(chunk) == 0 or len(all_others) == 0:
                    continue

                #if there is a flare, we don't want to consider it if it has fewer than 10 counts
                if np.mean(chunk) > np.mean(all_others) and sum(chunk) < 10:
                    continue
                #if there is a dip, we don't want to consider it if the observation
                #has fewer than 30 counts
                elif np.mean(chunk) < np.mean(all_others) and self.total_counts < 30:
                    continue


                try:
                    t_stat,p_val = test_poisson_2indep(sum(chunk),len(chunk),sum(all_others),len(all_others),method='etest')
                except Exception as e:
                    p_val = 1

                if p_val < pthresh:
                    if not self.classification:
                        self.classification = True

                    interesting_chunks[i] = p_val

                    if debug:
                        print(f'''Detection in chunk {i}
    This chunk: {chunk}
    all_others: {all_others}
    p_val: {p_val}
    ''')

                else:
                    if debug:
                        print(f'''No detection in chunk {i}
    This chunk: {chunk}
    all_others: {all_others}
    p_val: {p_val}
    ''')

            #if there aren't interesting chunks, return None
            if not interesting_chunks:
                return None

            #pick out the lowest p-value from among the interesting chunks
            lowest_p = 1
            for k,v in interesting_chunks.items():
                if v < lowest_p:
                    lowest_p = v
                    best_chunk = k

            if debug:
                print(f'Of {len(interesting_chunks.keys())} detections, the most interesting is chunk {k} with p-val {interesting_chunks[k]}')

            return best_chunk

        repeat = True

        interesting_chunks = []
        detections = []

        z = 0

        while repeat:
            z += 1


            new_chunk = find_best_chunk(chunks)

            if new_chunk:
                #hard break after finding 5 interesting chunks
                if z > 5:
                    repeat = False

                interesting_chunks.append(new_chunk)
                i = new_chunk

                if i == len(chunk_edges) - 1:
                    out = [chunk_edges[i],len(binned_counts)-1]
                else:
                    out = [chunk_edges[i],chunk_edges[i+1]]

                for index in out:
                    detections.append(binned_time[index])

                #remove the chunk from binned time
                try:
                    del binned_time[out[0]:out[1]+1]
                    del binned_counts[out[0]:out[1]+1]
                except IndexError as e:
                    print(f"Error deleting index [{out[0]}:{out[1]+1}]")
                    print(binned_time)
                    print(len(binned_time))
                    print(binned_counts)
                    print(len(binned_counts))

                    raise e

                chunk_edges,chunks = make_chunks(binned_counts)

            else:
                repeat = False


        if not interesting_chunks:
            return None
        else:
            if self.t_interest is None:
                self.t_interest = detections
            else:
                self.t_interest.extend(detections)

        if debug:
            print(f'Detections: {detections}')
            print(f'out:{out}')
            print(f'Binned counts: {binned_counts}')
            print(f'Binned times: {binned_time}')

        return detections

    #the same as the above but the chunks are made with Bayesian blocks
    def e_test_repeat_blocks(self,pthresh=2.867E-7,binsize=1000, debug = False):

        #default pthresh set to 5 sigma detection level

        if self.total_counts < 3:
            return None

        #first we need to chunk up the data
        #AKA dividing the data into chunks where each chunk is delimitted
        #by the light curve passing through the median of the dataset
        binned_time,binned_counts = self.make_binned_counts(binsize)

        t0 = binned_time[0]

        if binned_counts is None:
            return None

        #takes the binned counts array and turns it into chunks using the
        #bayesian blocks algorithm
        #returns the edges of the chunks, the chunks themselves
        def make_chunks(binned_counts):
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                chunk_edges = bayesian_blocks(binned_time,binned_counts,fitness='events',p0=0.001)


            #convert the chunk_edges from time space to index space
            chunk_edges = np.array(chunk_edges)/t0 - 1
            chunk_edges = chunk_edges[:-1]

            #turning the chunk edges into the chunks themselves
            #chunks are an array of sub arrays,
            #each sub array is a list of the values of binned_counts in the chunk

            #initialize the array of empty sub arrays
            chunks = [[] for i in range(len(chunk_edges))]

            for i,val in enumerate(binned_counts):
                for j,edge in enumerate(chunk_edges):

                    if edge == chunk_edges[-1]:
                        chunks[j].append(val)
                        break
                    else:
                        if i >= chunk_edges[j] and i < chunk_edges[j+1]:
                            chunks[j].append(val)
                            break

            return (chunk_edges,chunks)

        chunk_edges,chunks = make_chunks(binned_counts)
        if debug:
            print(chunk_edges)
            print('\n')
            print(chunks)
            print(len(binned_counts))
            chunk_edges_times = (chunk_edges+1)*t0
            #return [(index+1)*t0 for index in chunk_edges]

        #when given the array of chunks, applies the e_test and returns the
        #chunk with the lowest probability among all chunks which are significant
        def find_best_chunk(chunks):
            #initialize the dict of interesting chunks which pass the significance test
            #key will be the index of the chunk
            #value will be the p_value
            interesting_chunks = {}

            #Loop through all the chunks, for each testing if the mean inside the chunk
            #is significantly different than the mean outside of it
            #interesting chunks get added to interesting_chunks
            for i,chunk in enumerate(chunks):
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

                if len(chunk) == 0 or len(all_others) == 0:
                    continue

                #if there is a flare, we don't want to consider it if it has fewer than 10 counts
                if np.mean(chunk) > np.mean(all_others) and sum(chunk) < 10:
                    continue
                #if there is a dip, we don't want to consider it if the observation
                #has fewer than 30 counts
                elif np.mean(chunk) < np.mean(all_others) and self.total_counts < 30:
                    continue


                try:
                    t_stat,p_val = test_poisson_2indep(sum(chunk),len(chunk),sum(all_others),len(all_others),method='etest')
                except Exception as e:
                    p_val = 1

                if p_val < pthresh:
                    if not self.classification:
                        self.classification = True

                    interesting_chunks[i] = p_val

                    if debug:
                        print(f'''Detection in chunk {i}
    This chunk: {chunk}
    all_others: {all_others}
    p_val: {p_val}
    ''')

                else:
                    if debug:
                        print(f'''No detection in chunk {i}
    This chunk: {chunk}
    all_others: {all_others}
    p_val: {p_val}
    ''')

            #if there aren't interesting chunks, return None
            if not interesting_chunks:
                return None

            #pick out the lowest p-value from among the interesting chunks
            interesting_chunks_sorted = {k: v for k, v in sorted(interesting_chunks.items(), key=lambda item: item[1])}
            lowest_p = next(iter(interesting_chunks_sorted.values()))
            best_chunk = next(iter(interesting_chunks_sorted.keys()))

            if debug:
                print(f'Of {len(interesting_chunks.keys())} detections, the most interesting is chunk {best_chunk} with p-val {interesting_chunks[best_chunk]}')

            return best_chunk

        repeat = True


        interesting_chunks = []
        detections = []

        z = 0

        #Steps:
            #1: run find_best_chunk to find the chunk which is the most significant
            #2a: If there is no significant chunk, stop
            #2b: If there is a significant chunk, add it to the list of interesting chunks
            #and add the corresponding times to detections
            #3: delete the chunk from the binned counts array
            #4a: If this process has been repeated 5 times already, stop
            #4b: Otherwise, repeat.
        while repeat:
            z += 1


            #step 1
            #new chunk is the index of the most signicifcant chunk
            #or none if no chunks are significant
            new_chunk = find_best_chunk(chunks)


            #step 2b
            if new_chunk is not None:


                #step 4a
                if z > 5:
                    repeat = False

                interesting_chunks.append(new_chunk)
                i = new_chunk

                if i == len(chunk_edges) - 1:
                    out = [chunk_edges[i],len(binned_counts)-1]
                else:
                    out = [chunk_edges[i],chunk_edges[i+1]]

                for index in out:
                    time = (index+1)*t0
                    detections.append(time)

                if debug:
                    self.make_fourpanel_plot(outdir='.',save=False,show=True,lines=chunk_edges_times,highlight=detections)


                #step 3
                try:
                    out = [int(i+.5) for i in out]
                    del binned_time[out[0]:out[1]+1]
                    del binned_counts[out[0]:out[1]+1]
                    try:
                        del chunks[new_chunk]
                        chunk_edges = np.delete(chunk_edges,new_chunk)
                    except Exception as e:
                        print(chunks)
                        print(chunk_edges)
                        print(type(chunk_edges))
                        raise e

                except IndexError as e:
                    print(f"Error deleting index [{out[0]}:{out[1]+1}]")
                    print(binned_time)
                    print(len(binned_time))
                    print(binned_counts)
                    print(len(binned_counts))

                    raise e

                chunk_edges,chunks = make_chunks(binned_counts)
                if debug:
                    chunk_edges_times = (chunk_edges+1)*t0

            #step 2a
            else:
                repeat = False


        if not interesting_chunks:
            return None
        else:
            if self.t_interest is None:
                self.t_interest = detections
            else:
                self.t_interest.extend(detections)

        if debug:
            print(f'Detections: {detections}')
            #print(f'out:{out}')
            #print(f'Binned counts: {binned_counts}')
            #print(f'Binned times: {binned_time}')

        return detections

    #the same as the above but the chunks are made with Bayesian blocks
    #no repeating behavior
    def e_test_blocks(self,pthresh=2.867E-7,binsize=1000, debug = False):
        #default pthresh set to 5 sigma detection level

        #takes the binned counts array and turns it into chunks using the
        #bayesian blocks algorithm
        #returns the edges of the chunks, the chunks themselves
        def make_chunks(binned_counts):
            chunk_edges = bayesian_blocks(binned_time,binned_counts,fitness='events',p0=0.001)


            #convert the chunk_edges from time space to index space
            chunk_edges = np.array(chunk_edges)/t0 - 1
            chunk_edges = chunk_edges[:-1]

            #turning the chunk edges into the chunks themselves
            #chunks are an array of sub arrays,
            #each sub array is a list of the values of binned_counts in the chunk

            #initialize the array of empty sub arrays
            chunks = [[] for i in range(len(chunk_edges))]

            for i,val in enumerate(binned_counts):
                for j,edge in enumerate(chunk_edges):

                    if edge == chunk_edges[-1]:
                        chunks[j].append(val)
                        break
                    else:
                        if i >= chunk_edges[j] and i < chunk_edges[j+1]:
                            chunks[j].append(val)
                            break

            return (chunk_edges,chunks)

        #when given the array of chunks, applies the e_test and returns the
        #chunks which are below the pvalue
        def find_best_chunk(chunks):
            #initialize the dict of interesting chunks which pass the significance test
            #key will be the index of the chunk
            #value will be the p_value
            interesting_chunks = {}

            #Loop through all the chunks, for each testing if the mean inside the chunk
            #is significantly different than the mean outside of it
            #interesting chunks get added to interesting_chunks
            for i,chunk in enumerate(chunks):
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

                if len(chunk) == 0 or len(all_others) == 0 or len(chunk) == 1:
                    continue

                #if there is a flare, we don't want to consider it if it has fewer than 10 counts
                if np.mean(chunk) > np.mean(all_others) and sum(chunk) < 10:
                    continue
                #if there is a dip, we don't want to consider it if the observation
                #has fewer than 30 counts
                elif np.mean(chunk) < np.mean(all_others) and self.total_counts < 30:
                    continue


                if debug:
                    print(f'Trying E-test on chunk {i}')
                    print(f'test_poisson_2indep({sum(chunk)},{len(chunk)},{sum(all_others)},{len(all_others)},method="etest")')

                try:
                    t_stat,p_val = poisson_twosample(sum(chunk),len(chunk),sum(all_others),len(all_others),method='score')
                except Exception as e:
                    p_val = 1
                    raise e

                if p_val < pthresh:
                    if not self.classification:
                        self.classification = True

                    interesting_chunks[i] = p_val

                    if debug:
                        print(f'''Detection in chunk {i}
    This chunk: {chunk}
    all_others: {all_others}
    p_val: {p_val}
    ''')

                else:
                    if debug:
                        print(f'''No detection in chunk {i}
    This chunk: {chunk}
    all_others: {all_others}
    p_val: {p_val}
    ''')

            if not interesting_chunks:
                return None
            else:
                return [i for i in interesting_chunks.keys()]

        if self.total_counts < 3:
            return None

        #first we need to chunk up the data
        #AKA dividing the data into chunks where each chunk is delimitted
        #by the light curve passing through the median of the dataset

        if binsize is None:
            binned_time = self.times - self.start_time
            binned_counts = self.counts
        else:
            binned_time,binned_counts = self.make_binned_counts(binsize)

        t0 = binned_time[0]

        if binned_counts is None:
            return None


        if debug:
            print('Making chunks')

        chunk_edges,chunks = make_chunks(binned_counts)

        if len(chunk_edges) == 2:
            return None

        if debug:
            #print(chunk_edges)
            #print('\n')
            print('Chunks:')
            print(chunks)
            #print(len(binned_counts))
            chunk_edges_times = (chunk_edges+1)*t0
            #return [(index+1)*t0 for index in chunk_edges]


        detections = []

        if debug:
            print('finding_chunks')

        interesting_chunks = find_best_chunk(chunks)

        if interesting_chunks is None:
            return None

        for i in interesting_chunks:
            if i == len(chunk_edges) - 1:
                out = [chunk_edges[i],len(binned_counts)-1]
            else:
                out = [chunk_edges[i],chunk_edges[i+1]]

            for index in out:
                time = (index+1)*t0
                detections.append(time)


        if debug:
            self.make_fourpanel_plot(outdir='.',save=True,show=False,lines=chunk_edges_times,highlight=detections)


        if not len(detections):
            return None
        else:
            if self.t_interest is None:
                self.t_interest = detections
            else:
                self.t_interest.extend(detections)

        if debug:
            print(f'Detections: {detections}')

        return detections

    #produces arrays of information about the source which are to be put into
    #a csv.
    def make_csv_line(self):
        if self.classification:
            interesting = 'true'
            interesting_times = ';'.join([str(i) for i in self.t_interest])
        else:
            interesting = 'false'
            interesting_times = ''

        if self.parent_name is not None:
            csv_line = [self.parent_name,self.obsid,self.ra,self.dec,self.total_counts,
                            self.start_time,self.end_time,self.duration,self.count_rate,
                            interesting,interesting_times]
        else:
            csv_line = [self.name,self.obsid,self.ra,self.dec,self.total_counts,
                            self.start_time,self.end_time,self.duration,self.count_rate,
                            interesting,interesting_times]

        return np.array(csv_line)

    #fills in self.HR with HR information following the method of HR_driver_const_counts
    def make_HR(self,N,evt_dir,errorval='90',divide_energy=2000,override=False,
                subtract_start=True):
        # TODO: Write in override
        '''
        print(f'NAME: {self.name}')
        print(f'REGION: {self.region}')
        sys.exit()
        '''

        temp_region = 'TEMP.reg'

        with open(temp_region,'w') as file:
            file.write(self.region)

        working_dir = f'{evt_dir}/{self.position_basic}'
        try:
            os.makedirs(working_dir)
        except:
            pass

        BEHR_DIR = '/Users/sethlarner/BEHR_contain/BEHR'
        #BEHR_DIR = '/data/reu/slarner/BEHR_contain/BEHR'

        print('\nLooking for srcflux products...')
        try:
            bkg_region = unglob(glob.glob(f'{working_dir}/*bkgreg.fits'),True)
        except:
            print('None found, running srcflux...')

            make_regions(evt_dir,self.position,f'{working_dir}/')
            bkg_region = unglob(glob.glob(f'{working_dir}/*bkgreg.fits'),True)

        evt = unglob(glob.glob(f'{evt_dir}/*evt2*'))



        src_region = 'temp_reg.reg'

        with open('temp_reg.reg','w') as f:
            f.write(self.region)


        print('Making BEHR bash file...')

        outfile = f'{working_dir}/BEHR_bash.txt'

        BEHR_outdir = f'{BEHR_DIR}/{self.obsid}/{self.position_basic}'

        subprocess.run(f'rm -rf {BEHR_outdir}',shell=True)
        os.makedirs(BEHR_outdir)

        make_behr(evt,src_region,bkg_region,divide_energy,BEHR_DIR,outfile,BEHR_outdir,N)

        run_BEHR(outfile)


        start_time = self.start_time


        time,uppers,lowers,med = plot_BEHR_constcounts(BEHR_outdir,N,self.position_basic,self.obsid,f'./{evt}',src_region,lines=None,start_time=start_time,show=False,save=False)

        time = np.array(time).astype('float64')
        uppers = np.array(uppers).astype('float64')
        lowers = np.array(lowers).astype('float64')
        med = np.array(med).astype('float64')

        self.HR = np.column_stack((time,med,uppers,lowers))

        os.remove('temp_reg.reg')

        return

    #makes a plot of HR and LC, saves it in outdir
    #control output directory with outdir
    def plot_HR_and_lc(self,binsize,outdir):
        time = self.HR[::,0]
        med = self.HR[::,1]
        uppers = self.HR[::,2]
        lowers = self.HR[::,3]

        binned_time,binned_counts = self.make_binned_counts(binsize)

        fig,ax = plt.subplots()

        #first plot the lightcurve
        ax.step(binned_time,binned_counts)
        ax.set_xlabel('Time (s)')
        ax.set_ylabel('Counts/bin')

        ax2 = ax.twinx()
        ax2.plot(time,med,'g-')
        ax2.fill_between(time,lowers,uppers,step='mid',color='g',alpha=.2)

        if self.parent_name is not None:
            plt.savefig(f'{self.parent_name}_HR_LC.pdf')
        else:
            plt.savefig(f'{self.name}_HR_LC.pdf')

if __name__ == '__main__':
    try:
        debug = 't' in sys.argv[1] or 'T' in sys.argv[1]
    except IndexError:
        debug = False


    lc = np.loadtxt('./Survey/completed/NGC5194/textfiles/J132952.6934+471142.6143_13814_lc.fits.txt',skiprows=1)
    obsid = '13814'
    position = '132952.6934+471142.6143'

    src = Source(lightcurve=lc,obsid=obsid,position=position)

    times = src.e_test_blocks(binsize=1000,debug=debug)

    #src.make_fourpanel_plot(outdir='.',save=False,show=True,lines=times)
