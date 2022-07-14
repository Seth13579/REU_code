#functions and classes for running wavdetect on an obsid

from ciao_contrib.runtool import *
import sys
import os
import glob
import numpy as np
from source_class import Source
from re import sub
import pandas as pd
from astropy.io import fits
from astropy.table import Table

#a class representing all the regions in a given observation
class Obsid_all_regions:

    def __init__(self,obsid,all_regions):

        #the obsid this object represents
        self.obsid = obsid

        #a list of all the regions
        self.all_regions = all_regions

#a class representing a single region
class Region:

    def __init__(self, path, regtext, ra, dec, a, b, obsid, evt):
        #the path to the region file
        self.path = path

        #the actual string of the region, so that it can be recreated as needed
        self.regtext = regtext

        #the ra of the center of the region
        #in sexigesimal units
        self.ra = ra

        #the dec of the center of the region
        #in sexigesimal units
        self.dec = dec

        #the semimajor axis of the region
        self.a = a

        #the semiminor axis of the region
        self.b = b

        #the obsid of the region
        self.obsid = obsid

        #the path to the events file of the region
        self.evt = evt

        ra_mod = sub('\:','',self.ra)
        dec_mod = sub('\:','',self.ra)

        self.name = f'J{ra_mod}{dec_mod}'

        '''
        #used by other code, not meant to be interacted with by the user
        #tracks if the region has been matched to a source
        self.matched = False
        '''

    def make_lc(self,outdir='.'):

        dmextract.punlearn()

        dmextract.infile = f'{self.evt}[sky=region({self.path})][bin time=::3.24104]'
        dmextract.outfile = f'{outdir}/{self.name}_{self.obsid}_lc.fits'
        dmextract.opt = 'ltc1'
        dmextract.clobber = 'yes'

        dmextract()

        # list of columns
        cols = "TIME_BIN TIME_MIN TIME TIME_MAX COUNTS STAT_ERR AREA EXPOSURE COUNT_RATE COUNT_RATE_ERR"
        cols = cols.split(" ")

        # accessing fits data
        hdu_list = fits.open(f'{outdir}/{self.name}_{self.obsid}_lc.fits', memmap=True)
        evt_data = Table(hdu_list[1].data)

        # initialising DataFrame
        df = pd.DataFrame()

        # writing to dataframe
        for col in cols:
            df[col] = list(evt_data[col])

         # writing to file
        df.to_csv(f"{outdir}/TEMP_{self.name}_{self.obsid}_lc.fits.txt", index=False, sep=" ")

        lc = np.loadtxt(f"{outdir}/TEMP_{self.name}_{self.obsid}_lc.fits.txt", skiprows = 1)

        #cleanup
        os.remove(f"{outdir}/TEMP_{self.name}_{self.obsid}_lc.fits.txt")
        os.remove(f'{outdir}/{self.name}_{self.obsid}_lc.fits')

        hdu_list.close()

        return lc

    #method to make a source object from the region object
    def make_source(self):
        return Source(lightcurve=self.make_lc(),obsid=self.obsid, position=f'{self.ra}{self.dec}')

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
        print('No files found')
        raise Exception
    return str(arr[0]).strip("'[]'")

def detect(dir):#detect is used to run fluximage and wavdetect in sequence on an obsid
    #needs a level 2 directory with a chandra evt2 file
    fluximage.punlearn()
    evt = glob.glob(f'{dir}/*evt2*')
    evt = unglob(evt,True)
    fluximage.infile= evt
    fluximage.outroot = f'{dir}/detect'
    fluximage.bands = '0.3:7.5:2.3'
    fluximage.psfecf=0.9 #one sigma of 2D gaussian (see 'running wavdetect')
    fluximage.clobber = 'yes'
    fluximage()

    wavdetect.punlearn()
    img = glob.glob(f'{dir}/*thresh.img')
    img = unglob(img)
    wavdetect.infile = img
    psf = glob.glob(f'{dir}/*thresh.psfmap')
    psf = unglob(psf)
    wavdetect.psffile = psf
    wavdetect.outfile = f'{dir}/detect_src.fits'
    wavdetect.imagefile = f'{dir}/detect_imgfile.fits'
    wavdetect.regfile = f'{dir}/detect_src.reg'
    wavdetect.defnbkgfile = f'{dir}/detect_nbgd.fits'
    wavdetect.scellfile = f'{dir}/detect_scell.fits'
    wavdetect.clobber = 'yes'
    wavdetect()

    return

#when given as argument a region from wavdetect, returns an Obsid_all_regions
#object
def process_wavdetect(obsid,region_dir,region_file):
    evt = glob.glob(f'{region_dir}/*evt2*')
    evt = unglob(evt)

    regphystocel.punlearn()
    regphystocel.infile = region_file
    regphystocel.outfile = f'{region_dir}/wavdetect_wcs_reg.fits'
    regphystocel.wcsfile = evt
    regphystocel.clobber = 'yes'
    regphystocel()

    regiontxt = np.loadtxt(f'{region_dir}/wavdetect_wcs_reg.fits',dtype='str',skiprows=3)

    #incredibly ugly way of taking the data in the array and mapping it into
    #multiple region objects
    regions = [Region(region_file,line,line[7:].strip('()').split(',')[0], line[7:].strip('()').split(',')[1],max(line[7:].strip('()').split(',')[2],line[7:].strip('()').split(',')[3]).rstrip('"'), min(line[7:].strip('()').split(',')[2],line[7:].strip('()').split(',')[3]).rstrip('"'),obsid,evt) for line in regiontxt]

    out = Obsid_all_regions(obsid,regions)

    return out

if __name__ == '__main__':
    obsid = sys.argv[1]
    dir = sys.argv[2]

    detect(dir)

    regions = process_wavdetect(obsid,dir,'detect_src.fits').all_regions

    #print(regions[0].regtext)
    print(regions[0].ra)
    print(regions[0].dec)
    #print(regions[0].a)
    #print(regions[0].b)
    #print(regions[0].obsid)
