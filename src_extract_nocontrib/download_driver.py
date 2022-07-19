from region_matching import *
from run_wavdetect import *
from galaxy_dwnld import *
from source_class import Source, Source_All
import numpy as np
import os
import sys
import glob


#function to handle computationally expensive steps
#called in the multithreading step
#multithreading not yet functional
def subreprocess(dir):
    obsid = dir.split('/')[-1]

    working_dir = f'./{dir}/repro'

    # TODO: Multithread this
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

def process_galaxy(galaxy,raidus,rerun):
    if not rerun:
        #start by downloading and reprocessing the data
        try:
            os.makedirs(f'./{galaxy}')
        except:
            pass

        os.chdir(f'./{galaxy}')

        print('Downloading...')
        galaxy_download(galaxy,radius)

        print('Splitting...')
        split_interleaved()

        print('Reprocessing...')
        repro()

    else:
        os.chdir(f'./{galaxy}')


    ############
    #now run wavdetect on each obsid
    dirs = [i for i in os.listdir(os.getcwd()) if '.' not in i and i != 'textfiles']

    print('Running wavdetect...')

    all_regions_in_galaxy = []
    for i,dir in enumerate(dirs):
        print(f'\nRunning on {dir}, {i+1} of {len(dirs)}')

        all_regions_in_galaxy.append(subreprocess(dir))

    galaxy_object = Galaxy(galaxy,all_regions_in_galaxy)

    print('Matching regions...')
    all_sources_in_galaxy = galaxy_object.region_match()

    try:
        os.makedirs('./textfiles')
    except:
        pass

    os.chdir('./textfiles')

    print('Producing textfiles...')
    for all_source in all_sources_in_galaxy:
        all_source.save_lcs()

if __name__ == '__main__':
    #the galaxy to download and process
    galaxy = sys.argv[1]

    #the radius in arc minutes over which to search
    radius = sys.argv[2]

    rerun = 't' in sys.argv[3] or "T" in sys.argv[3]

    process_galaxy(galaxy,radius,rerun)
