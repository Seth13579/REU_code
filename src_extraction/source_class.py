import numpy as np
from scipy import stats
import pickle
from re import sub
import os


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

        #an array of source objects which represent all non-zero observations
        #of the source
        self.obs_nonzero = [i for i in obs if not i.is_zero()]

        #the average count rate across all observations of the source
        #units of ks^-1
        self.average_count_rate = 1000*sum([i.total_counts for i in self.obs])/sum([i.duration for i in self.obs])

        #the average count rate across all observations of the source with
        #non-zero counts
        #units of ks^-1
        try:
            self.av_count_rate_nozeros = 1000*sum([i.total_counts for i in self.obs_nonzero])/sum([i.duration for i in self.obs_nonzero])
        except:
            self.av_count_rate_nozeros = 0

    def save(self,outdir='.'):
        outfile = f'{outdir}/SOURCE_ALL_{ra}_{dec}.pkl'

        with open(outfile,'w') as outp:
            pickle.dump(self, outp, pickle.HIGHEST_PROTOCOL)

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

        return

    def save_lcs(self,outdir='.'):
        for source in self.obs:
            lc = source.lightcurve

            filename = f'{outdir}/{self.name}_{source.obsid}_lc.fits.txt'

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
    def __init__(self, lightcurve=None,obsid=None, position=None,classification=None,
                t_interest=None, region=None, region_object=None):
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

        #optional region text describing the source
        self.region = region

        #the region object which created this source
        self.region_object = region_object


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

    #returns the binned counts array binned at binsize
    def make_binned_counts(self,binsize):
        bin_length = int(binsize/3.24014)
        bin_edges = [self.times[i] for i in range(len(self.times)) if i%bin_length == 0]
        bin_edges.append(self.times[-1])
        binned_counts = stats.binned_statistic(self.times,self.counts,statistic='sum',bins=bin_edges)[0]

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

if __name__ == '__main__':
    file = 'test_lc.txt'

    lc = np.loadtxt(file, skiprows = 1)

    test_source = Source(lc)

    print(test_source.count_cut(4))
