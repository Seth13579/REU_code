import numpy as np
import pickle
from re import sub,search
import re
import os
import matplotlib.pyplot as plt
import glob
import sys

#for the E-test
from statsmodels.stats.rates import test_poisson_2indep

#HR dependencies
from extract_counts import *
from BEHR_countbins import *
from run_BEHR import *
from create_behr_files import region_area

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

        with open(outfile,'wb') as outp:
            pickle.dump(self, outp)

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
    def __init__(self, lightcurve=None,obsid=None, position=None,classification=False,
                t_interest=None, region=None, region_object=None):
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



        return binned_times,binned_counts

    #makes a four panel light curve for the source
    #one cumulative plot
    #three light curves, defined by the bin sizes in binsizes
    def make_fourpanel_plot(self,outdir,binsizes=[500,1000,2000],save=True,
                            show=False,lines=None):
        if lines == None:
            lines = self.t_interest

        name = f'{re.sub(":","",self.position)}_{self.obsid}'

        binned_times_1,binned_counts_1 = self.make_binned_counts(binsizes[0])
        binned_times_2,binned_counts_2 = self.make_binned_counts(binsizes[1])
        binned_times_3,binned_counts_3 = self.make_binned_counts(binsizes[2])

        trim_low,trim_high = trimzeros(self.counts)
        trim_counts = self.counts[trim_low:trim_high]
        trim_time = self.times[trim_low:trim_high]

        cumo = np.cumsum(trim_counts)

        fig, axs = plt.subplots(2,2)
        cumo_plt = axs[0,0]
        bin1_plt = axs[0,1]
        bin2_plt = axs[1,0]
        bin3_plt = axs[1,1]

        cumo_plt.plot(trim_time-trim_time[0],cumo)
        cumo_plt.set_title('Cumulative Counts')

        bin1_plt.step(binned_times_1,binned_counts_1)
        bin1_plt.set_title(f'Bin sizes: {round(binsizes[0])}')

        bin2_plt.step(binned_times_2,binned_counts_2)
        bin2_plt.set_title(f'Bin sizes: {round(binsizes[1])}')

        bin3_plt.step(binned_times_3,binned_counts_3)
        bin3_plt.set_title(f'Bin sizes: {round(binsizes[2])}')

        cumo_plt.set(ylabel='Counts')
        bin2_plt.set(xlabel='Time (s)',ylabel='Counts')
        bin3_plt.set(xlabel='Time (s)')


        fig.suptitle(name)

        if lines is not None:
            for ax in axs.flat:
                ax.vlines(lines,ymin=ax.get_ylim()[0],ymax=ax.get_ylim()[1],colors='r',linestyles='dotted')

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

    #produces arrays of information about the source which are to be put into
    #a csv.
    def make_csv_line(self):
        if self.classification:
            interesting = 'true'
            interesting_times = ';'.join([str(i) for i in self.t_interest])
        else:
            interesting = 'false'
            interesting_times = ''

        csv_line = [self.name,self.obsid,self.ra,self.dec,self.total_counts,
                        self.start_time,self.end_time,self.duration,self.count_rate,
                        interesting,interesting_times]

        return np.array(csv_line)

    #fills in self.HR with HR information following the method of HR_driver_const_counts
    def make_HR(self,N,evt_dir,errorval='90',divide_energy=2000,override=False,subtract_start=True):
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

        print('Looking for srcflux products...')
        try:
            bkg_region = unglob(glob.glob(f'{working_dir}/*bkgreg.fits'),True)
        except:

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
        print(f'binned_time[0:4]:{binned_time[0:4]}')

        fig,ax = plt.subplots()

        #first plot the lightcurve
        ax.step(binned_time,binned_counts)
        ax.set_xlabel('Time (s)')
        ax.set_ylabel('Counts/bin')

        ax2 = ax.twinx()
        ax2.plot(time,med,'g-')
        ax2.fill_between(time,lowers,uppers,step='mid',color='g',alpha=.2)

        plt.savefig(f'{self.name}_HR_LC.pdf')

if __name__ == '__main__':
    lc = np.loadtxt('NGC1313/textfiles/J031822.1213-663603.6133_2950_lc.fits.txt',skiprows=1)
    regtext='ellipse(03:18:22.2404,-66:36:3.7656,7.06016",5.4328",54.4779)'
    src = Source(lightcurve=lc,obsid='2950',position='03:18:22.2404-66:36:03.7656',region=regtext)
    src.make_HR(16,'NGC1313/2950/repro')

    src.plot_HR_and_lc(1000,'.')

    print(src.HR[0,::])
