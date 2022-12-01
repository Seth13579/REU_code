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


#recursion depth context manager
class recursion_depth:
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



#function to handle computationally expensive steps
# TODO: multithreading
def subreprocess(dir):
    obsid = dir.split('/')[-1]

    working_dir = f'./{dir}/primary'

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

#carrys out the steps of region making and matching to create a galaxy class object
#assumes that the data is downloaded and reprocessed
#then produces textfiles in ./{galaxy}/textfiles
#then produces plots in
def process_galaxy(galaxy_name):
    os.chdir(f'./{galaxy_name}')

    #now run wavdetect on each obsid
    dirs = [i for i in os.listdir(os.getcwd()) if re.search('\d',i)]

    if len(dirs) == 0:
        raise Exception

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

    galaxy = Galaxy(galaxy_name,all_regions_in_galaxy)

    ###########################
    ###########################
    print('\nMatching regions...')

    galaxy.region_match()

    all_sources_in_galaxy = galaxy.matches

    for source in all_sources_in_galaxy:
        source.optimal_rename()


    try:
        os.makedirs('./textfiles')
    except:
        pass

    os.chdir('./textfiles')
    os.system('rm *.txt')

    print('\nProducing textfiles...')
    for all_source in all_sources_in_galaxy:
        all_source.save_lcs()

    return
if __name__ == '__main__':
    #first command line argument controls which galaxies to run on
    #can either list galaxies, seperated by commas
    #or use 'all' to run on all dirs in cwd

    if 'all' in sys.argv[1] or 'ALL' in sys.argv[1] or 'All' in sys.argv[1]:
        galaxies = [i for i in os.listdir(os.getcwd()) if '.' not in i and 'textfiles' not in i and 'errored' not in i and 'completed' not in i]

    else:
        galaxies = sys.argv[1].split(',')

    errors = []

    cwd = os.getcwd()

    try:
        os.makedirs('./errored')
    except:
        pass

    try:
        os.makedirs('./completed')
    except:
        pass

    for i,galaxy in enumerate(galaxies):
        print('***********')
        print(f'PROCESSING {galaxy}, {i+1} OF {len(galaxies)}')
        print('***********')

        try:
            process_galaxy(galaxy)
        except Exception as e:
            errors.append(galaxy)
            os.chdir(cwd)
            os.system(f'mv {galaxy} ./errored')
            print(f'{galaxy} erroed, moved to ./errored')
            #raise e
        else:
            os.chdir(cwd)
            os.system(f'mv {galaxy} ./completed')
            print(f'{galaxy} completed without error, moved to ./completed')

    with open('Error_doc.txt','w') as f:
        for gal in errors:
            f.write(gal)
