#functions for downloading all the obsids for a given galaxy name

from ciao_contrib.runtool import chandra_repro
import sys
import os

def galaxy_download(galaxy,radius):
    search_str = f'find_chandra_obsid "{galaxy}" radius={radius} grating=none instrument=acis download=all'

    os.system(search_str)

    return

def repro(dir = '.'):

    os.chdir(dir)
    for _dir in os.listdir('.'):
        chandra_repro.punlearn()

if __name__ == '__main__':
    galaxy = 'NGC 1313'

    galaxy_download(galaxy,20)
    repro('*')
