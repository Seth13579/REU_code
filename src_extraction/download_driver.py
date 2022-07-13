from region_matching import *
from run_wavdetect import *
from galaxy_dwnld import *
from source_class import Source, Source_All
import numpy as np
import os
import sys
import glob


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

        print('Reprocessing...')
        repro('*')

    else:
        os.chdir(f'./{galaxy}')

    ############
    #now run wavdetect on each obsid
    all_regions_in_galaxy = []
    for dir in [i for i in os.listdir(os.getcwd()) if '.' not in i and i != 'textfiles']:
        obsid = dir.split('/')[-1]

        print(f'Running wavdetect on {obsid}')

        working_dir = f'./{dir}/repro'

        if not rerun:
            detect(working_dir)

        wavdetect_region = unglob(glob.glob(f'{working_dir}/*detect_src.reg*'))

        all_regions_in_obsid = process_wavdetect(obsid,working_dir,wavdetect_region)

        all_regions_in_galaxy.append(all_regions_in_obsid)

    galaxy_object = Galaxy(galaxy,all_regions_in_galaxy)

    print('Matching regions...')
    all_sources_in_galaxy = galaxy_object.region_match()

    if not rerun:
        os.makedirs('./textfiles')

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
