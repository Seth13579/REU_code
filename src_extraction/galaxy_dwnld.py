#functions for downloading all the obsids for a given galaxy name
import sys
import os
import re
import glob

def galaxy_download(galaxy,radius):
    os.system('punlearn find_chandra_obsid')

    search_str = f'find_chandra_obsid "{galaxy}" radius={radius} grating=none instrument=acis download=all'

    os.system(search_str)

    return

#downloads the given obsid to the current working directory
def obsid_download(obsid):
    #os.system('punlearn download_chandra_obsid')
    os.system(f'download_chandra_obsid {obsid}')

    return


#returns true if the given directory is an interleaved observation
def check_for_interleaved(dir):
    evt2 = glob.glob(f'{dir}/primary/*evt2*')
    evt1 = glob.glov(f'{dir}/secondary/*evt1*')
    if len(evt2) > 1 or len(evt1) > 1:
        return True
    else:
        return False

#runs splitobs on the given directory
def run_split(dir):
    os.system('punlearn splitobs')
    os.system(f'splitobs indir={dir} outroot={dir} clobber=yes')
    os.system(f'rm -rf {dir}')

    return

#splits all interleaved obsids in the current working directory
def split_interleaved():
    dirs = [i for i in os.listdir(os.getcwd()) if re.search('\d\d\d',i)]

    for dir in dirs:
        if check_for_interleaved(dir):
            run_split(dir)

    return

#reprocess all the obsids in the current working directory
def repro():
    #regex used to find any directory with at least three digits in a row
    #assuming that these are the obsids
    dirs = [i for i in os.listdir(os.getcwd()) if re.search('\d\d\d',i)]

    for i,dir in enumerate(dirs):
        os.system('punlearn chandra_repro')
        repro_string = f'chandra_repro {dir} outdir="" verbose=0 clobber=yes'

        print(f'Reprocessing {dir}, {i+1} of {len(dirs)}')

        os.system(repro_string)

    return

if __name__ == '__main__':
    galaxy = 'NGC 1313'

    galaxy_download(galaxy,20)
    repro('*')
