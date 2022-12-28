#The beginning of the full process
#Dip And Flare Finder

from Source import Source, Source_All
import os
from run_wavdetect import *
from region_matching import *
import re
import pickle as pkl
import sys
import glob
import sys

import matplotlib as mpl
#mpl.use('agg')

class recursion_depth:
    '''recursion depth context manager'''
    def __init__(self, limit):
        self.limit = limit
        self.default_limit = sys.getrecursionlimit()
    def __enter__(self):
        sys.setrecursionlimit(self.limit)
    def __exit__(self, type, value, traceback):
        sys.setrecursionlimit(self.default_limit)

def restore_match(dir):
    cwd = os.getcwd()

    os.chdir(dir)

    source_alls = glob.glob('SOURCE_ALL*')

    galaxyfile = unglob(glob.glob('*galaxy_obj.pkl'))

    with open(galaxyfile,'rb') as f:
        galaxy = pkl.load(f)

    all_sources = []

    for src in source_alls:
        with open(src,'rb') as f:
            all_sources.append(pkl.load(f))

    galaxy.matches = all_sources

    os.chdir(cwd)
    return galaxy

def subreprocess(dir):
    '''function to handle computationally expensive steps'''
    obsid = dir.split('/')[-1]

    working_dir = f'./{dir}/repro'

    try:
        print("Looking for wavdetect product...")
        wavdetect_region = unglob(glob.glob(f'{working_dir}/*detect_src.reg*'))
    except:
        print('Trying wavdetect...')
        detect(working_dir)
        wavdetect_region = unglob(glob.glob(f'{working_dir}/*detect_src.reg*'))
    else:
        print('Found!')
    finally:
        all_regions_in_obsid = process_wavdetect(obsid,working_dir,wavdetect_region)

    return all_regions_in_obsid

def download_galaxy(name,ra,dec,D25):
    cwd = os.getcwd()

    try:
        os.makedirs(name)
    except:
        pass

    os.chdir(name)

    find_chandra_obsid.punlearn()

    find_chandra_obsid.arg = ra
    find_chandra_obsid.dec = dec
    find_chandra_obsid.radius = 2*D25 #must be in arcmin
    find_chandra_obsid.grating = 'none'
    find_chandra_obsid.instrument='acis'
    find_chandra_obsid.download = 'none'

    output = find_chandra_obsid().splitlines()

    if len(output) == 1: #when there are no observations found
        os.chdir(cwd)
        os.rmdir(name)
        return False

    obsids = [i.split()[0] for i in output[1:]]

    for i,ob in enumerate(obsids):
        print(f'Downloading obsid {ob}, {i+1} of {len(obsids)}')
        os.system(f'download_chandra_obsid {obsid} --exclude evt1,evt1a,vvref')

    os.chdir(cwd)

    return True

#carrys out the steps of region making and matching to create a galaxy class object
#assumes that the data is downloaded and reprocessed
#then produces textfiles in ./{galaxy}/textfiles
#then produces plots in
def process_galaxy(galaxy_name,ra,dec,D25):
    os.chdir(f'./{galaxy_name}')

    #now run wavdetect on each obsid
    dirs = [i for i in os.listdir(os.getcwd()) if re.search('\d',i)]

    ###########################
    ###########################
    print('\nRunning wavdetect...')

    all_regions_in_galaxy = []
    for i,dir in enumerate(dirs):
        print(f'\nRunning on {dir}, {i+1} of {len(dirs)}')

        try:
            new = subreprocess(dir)
        except:
            new = None

        if new is not None:
            all_regions_in_galaxy.append(new)

    galaxy = Galaxy(galaxy_name,all_regions_in_galaxy,D25=D25,ra=ra,dec=dec)

    ###########################
    ###########################
    print('\nCulling regions...')
    galaxy.eliminate_outside()

    ###########################
    ###########################
    print('\nMatching regions...')

    galaxy.region_match()

    all_sources_in_galaxy = galaxy.matches

    for source in all_sources_in_galaxy:
        source.optimal_rename()


    ###########################
    ##########################
    print('\nSaving matches...')

    try:
        os.makedirs('./matches')
    except:
        pass

    os.chdir('./matches')
    os.system('rm *.pkl')

    #Cannot save the whole galaxy object
    #Instead, I save all the source all objects by themselves and then
    #a dummy galaxy object which contains all the other data
    #This information can be used to restore the state of this program after
    #matches are made

    with recursion_depth(5000):
        for all_source in all_sources_in_galaxy:
            try:
                all_source.save()
            except RecursionError as e:
                pass

    dummy_galaxy = Galaxy(galaxy.name,galaxy.obsid_region_list)

    with open(f'{galaxy_name}_galaxy_obj.pkl','wb') as f:
        pkl.dump(dummy_galaxy,f)

    os.chdir('../')

    ###########################
    ###########################
    try:
        os.makedirs('./textfiles')
    except:
        pass

    os.chdir('./textfiles')
    os.system('rm *.txt')

    print('\nProducing textfiles...')
    for all_source in all_sources_in_galaxy:
        all_source.save_lcs()

    os.chdir('../')
    try:
        os.makedirs('./detections_lcs')
    except:
        pass
    os.chdir('./detections_lcs')

    ###########################
    ###########################
    print('\nDetecting dips and flares...')
    lines = []
    total = len([source for all_source in all_sources_in_galaxy for source in all_source.obs])

    i = 0
    for all_source in all_sources_in_galaxy:
        for source in all_source.obs:
            print(f'Testing {i} of {total}')

            try:
                #runs the repeating chunked e-test method
                #if an interesting chunk is found:
                #updates source.t_interest to the edges of the interesting times
                #updates source.classification to True
                source.e_test_repeat(binsize=1000)
            except:
                source.t_interest = None
                source.classification = False

            if source.classification:
                source.make_fourpanel_plot(outdir='.')

            csv_line = source.make_csv_line()

            lines.append(csv_line)

            i += 1

    ###########################
    ###########################
    print('\nSaving csv...')
    csv = np.row_stack(lines)
    header = 'NAME,OBSID,RA,DEC,TOTAL COUNTS,START TIME,END TIME,DURATION,COUNT RATE,INTERESTING,INTERESTING TIMES'

    os.chdir('../')

    np.savetxt('./all_csv.csv',csv,fmt='%s',delimiter=',',header=header)

    '''
    ###########################
    ###########################
    print('\nMaking HR plots...')
    try:
        os.makedirs('./HR_plots')
    except:
        pass

    os.chdir('./HR_plots')

    for all_source in all_sources_in_galaxy:
        for source in all_source.obs:

            if source.classification:
                evt_dir = f'../{source.obsid}/repro'

                try:
                    source.make_HR(16,evt_dir)
                    source.plot_HR_and_lc(1000,'.')
                except:
                    pass

    ###########################
    ###########################
    print('Saving sources...')
    os.chdir('../')
    try:
        os.makedirs('./sources')
    except:
        pass
    os.system('rm ./sources/*.pkl')
    os.chdir('./sources')

    for all_source in all_sources_in_galaxy:
        all_source.save()

    '''

    return


def lookup_galaxy(name,table_data):
    '''Looks up the given galaxy in the table, returns (ra (deg),dec (deg),D25)
    Raises exception if name not found in the table'''

    try:
        idx = np.where(table_data[::,0] == name)[0][0]
    except ValueError as e:
        print('Galaxy not found in table')

    ra_sex = ':'.join([table_data[idx,1],table_data[idx,2],table_data[idx,2]])
    dec_sex = ':'.join([''.join([table_data[idx,4],table_data[idx,5]]),table_data[idx,6],table_data[idx,7]])

    c = SkyCoord(ra_sex,dec_sex,unit=(u.hourangle, u.deg),frame='fk5')

    ra = c.ra.degree
    dec = c.dec.degree
    D25 = float(table_data[idx,8])

    return (ra,dex,D25)


if __name__ == '__main__':
    #table is a numpy file detailing all the information about the galaxies
    #it has the following columns
    #COL 0: GALAXY NAMES
    #COL 1: Hour of Right Ascension
    #COL 2: Minute of Right Ascension
    #COL 3: Second of Right Ascension
    #COL 4: Sign of the Declination
    #COL 5: Degree of Declination
    #COL 6: Arcminute of Declination
    #COL 7: Arcsecond of Declination
    #COL 8: D25 (arcmin)
    table_data = np.load(sys.argv[1])

    #second command line argument controls which galaxies to run on
    #can either list galaxies, seperated by commas
    #or use 'all' to run on all galaxies in the table
    if 'all' in sys.argv[2] or 'ALL' in sys.argv[2] or 'All' in sys.argv[2]:
        galaxies = table_data[::,0]

    else:
        galaxies = sys.argv[2].split(',')

    errors = []

    cwd = os.getcwd()

    #making the folders we will need.
    try:
        os.makedirs('./errored')
    except:
        pass

    try:
        os.makedirs('./completed')
    except:
        pass

    for i,galaxy in enumerate(galaxies):
        ra,dec,D25 = lookup_galaxy(galaxy,table_data)

        print('***********')
        print(f'DOWNLOADING {galaxy}, {i+1} OF {len(galaxies)}')
        print('***********')

        cont = download_galaxy(galaxy,ra,dec,D25)

        if cont:
            print('***********')
            print(f'PROCESSING {galaxy}, {i+1} OF {len(galaxies)}')
            print('***********')

            try:
                problematic = process_galaxy(galaxy,ra,dec,D25)
            except Exception as e:
                errors.append(galaxy)
                os.chdir(cwd)
                os.system(f'mv {galaxy} ./errored')
                print(f'{galaxy} erroed, moved to ./errored')
            else:
                os.chdir(cwd)
                os.system(f'mv {galaxy} ./completed')
                print(f'{galaxy} completed without error, moved to ./completed')
        else:
            print('NO MATCHED OBSERVATIONS')
            print('MOVING TO ERRORED')

            errors.append(galaxy)
            os.chdir(cwd)
            os.system(f'mv {galaxy} ./errored')
            print(f'{galaxy} erroed, moved to ./errored')

    with open('Error_doc.txt','w') as f:
        for gal in errors:
            f.write(gal)
