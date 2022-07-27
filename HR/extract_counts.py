import os
from ciao_contrib.runtool import *
import glob
import numpy as np


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

def download_obsid(obsid):
    if not os.path.isdir(obsid):
        os.system(f'download_chandra_obsid {obsid} --exclude vvref,evt1,evt1a')
    return

def make_regions(obsid,pos,outroot):
    evt2 = unglob(glob.glob(f'./{obsid}/primary/*evt2*'))

    try:
        #no need to run it if the products already exist
        src_region =unglob(glob.glob(f'{outroot}/*srcreg.fits'),force=True)
        bkg_region = unglob(glob.glob(f'{outroot}/*bkgreg.fits'),force=True)
    except:

        srcflux.punlearn()
        srcflux.infile = evt2
        srcflux.pos = pos
        srcflux.outroot = outroot
        srcflux.clobber = 'yes'
        #srcflux.psfmethod = 'arfcorr'
        srcflux()

#returns an array of counts in the selected band
def split_events(evt,reg,bin_size,band):
    if band == 'soft':
        energy = '500:2000'
    elif band == 'hard':
        energy = '2000:7500'
    else:
        raise Exception

    dmextract.punlearn()
    dmextract.infile = f'{evt}[energy={energy}][sky=region({reg})][bin time=::{bin_size}]'
    dmextract.clobber = 'yes'
    dmextract.opt = 'generic'
    dmextract.outfile = 'temp.fits'
    dmextract()

    dmlist.punlearn()
    dmlist.infile = 'temp.fits[cols counts]'
    dmlist.opt = 'data, clean'

    try:
        out = np.asarray(dmlist().split('\n')[1:]).astype('int')
    except Exception as e:
        print(dmlist().split('\n')[1:])
        raise e
    return out
