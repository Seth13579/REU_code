#functions for downloading all the obsids for a given galaxy name

from ciao_contrib.runtool import chandra_repro
import sys
import os
import re

def galaxy_download(galaxy,radius):
    search_str = f'find_chandra_obsid "{galaxy}" radius={radius} grating=none instrument=acis download=all'

    os.system(search_str)

    return



#reprocess all the obsids in the current working directory
def repro():
    dirs = [i for i in os.listdir(os.getcwd()) if not re.search('[a-zA-Z]',i)]

    for i,dir in enumerate(dirs):
        repro_string = f'chandra_repro {dir} outdir="" verbose=0 clobber=yes'

        print(f'Reprocessing {dir}, {i+1} of {len(dirs)}')

        os.system(repro_string)

if __name__ == '__main__':
    galaxy = 'NGC 1313'

    galaxy_download(galaxy,20)
    repro('*')
