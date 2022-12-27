#functions and classes for running wavdetect on an obsid

from ciao_contrib.runtool import *
import sys
import os
import glob
import numpy as np
from Source import Source,Source_All
from re import sub
import pandas as pd
from astropy.io import fits
from astropy.table import Table

def countcheck(evt,reg):
    '''When given a evt file and a region, returns the number of counts between .3 and 7.5 KeV'''

    dmlist.punlearn()
    dmlist.infile = f'{evt}[sky=region({reg}),energy=300:7500]'
    dmlist.opt = 'counts'
    return int(dmlist())


#a class representing all the regions in a given observation
class Obsid_all_regions:

    def __init__(self,obsid,all_regions):

        #the obsid this object represents
        self.obsid = obsid

        #a list of all the regions
        self.all_regions = all_regions

#a class representing a single region
class Region:

    def __init__(self, path, regtext, ra, dec, a, b, obsid, evt, all_reg_path, bkg_reg = None):
        #the path to the region file
        self.path = path

        #path to the wavdetect made all region path
        self.all_reg_path = all_reg_path

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
        dec_mod = sub('\:','',self.dec)

        self.name = f'J{ra_mod}{dec_mod}'

        # a list of match objects corresponding to the matches the region recieved
        self.matches = None

        #a match object corresponding to the best match the region object recieved
        self.best_match = None

        #The text of the background region
        #can be filled in with self.create_bkg_region
        self.bkg_reg = bkg_reg

    def create_bkg_region(self,ncounts = 50,outdir='.'):
        '''
        Creates a background region to go with the source region represented
        by this object.

        Background region is created by making an annulus around the source region.
        The outer radius is then increased in small steps until the region encircles
        ncounts. The region excludes all point sources detected by wavdetect.

        This region will be used to subtract the background from the source in
        order to (hopefully) reduce the impact of systematic effects.

        The resulting background region will be save as text in self.bkg_reg
        '''

        def make_bkg(src_reg,all_reg,size,outdir):
            '''
            Makes a background annulus centered on "src_reg" with outer radius
            equal to "size" arcmins excluding all the point source regions in
            "all_reg".

            Resulting region saved to "{outdir}/bkg.reg"
            Returns path to that file.
            '''

            outfile = f'{outdir}/bkg.reg'

            region = str(np.loadtxt(src_reg,dtype='str'))[8:-1]

            all_region = str(np.loadtxt(all_reg,dtype='str',skiprows=3)).rstrip('#').strip('[]').strip("''")
            all_region = all_region.split('\n')

            region_list = region.split(",")
            region_axis_1 = float(region_list[2][:-1])
            region_axis_2 = float(region_list[3][:-1])
            #have to convert the axis values from sec to min
            region_axis_1 *= 0.01667 #converts to arcminutes
            region_axis_2 *= 0.01667
            if region_axis_1 > region_axis_2:
                bkg_radius = region_axis_1
            else:
                bkg_radius = region_axis_2

            with open(outfile,'w') as out:
                bkg_reg = f"annulus({region_list[0]},{region_list[1]},{bkg_radius*1.5}',{(bkg_radius*1.5)+size}')"
                out.write(bkg_reg)
                out.write('\n')

                for line in all_region:
                    line = line.strip()
                    line = line.strip("''")
                    out.write(f'\n -{line}')

            return outfile

        stop = False
        i = 0

        evt = self.evt
        sec_reg = self.path
        all_reg = self.all_reg_path
        #List of radii to try, 1 arcsec to 1.5 arc minutes in 1 arcsec steps
        sizes = np.arange(.01667,1.5,.01667)

        while not stop:
            #repeatedly make backgrounds increasing in size until the counts
            #threshold is reached

            try:
                bkg = make_bkg(src_reg,all_reg,sizes[i],outdir)
                i+=1

                check = countcheck(evt,bkg)
                if check > limit:
                    stop=True

            except IndexError:
                stop = True

        with open(bkg,'r') as f:
            bkg_reg_text = f.read()

        os.remove(bkg)

        self.bkg_reg = bkg_reg_text

        return


    def make_lc(self,subtract=True,outdir='.'):

        self_region = f"{outdir}/TEMP_{self.name}_{self.obsid}_reg.reg"
        with open(self_region,'w') as reg_file:
            reg_file.write(self.regtext)

        dmextract.punlearn()

        dmextract.infile = f'{self.evt}[sky=region({self_region})][bin time=::3.24104]'

        #If the background is to be subtracted, need to make a region and then
        #Update the dmextract parameters
        if subtract:
            create_bkg_region()

            bkg_region = f"{outdir}/TEMP_{self.name}_{self.obsid}_BKG_reg.reg"

            with open(bkg_region,'w') as reg_file:
                reg_file.write(self.bkg_reg)

            dmextract.bkg = f'{self.evt}[sky=region({bkg_region})][bin time=::3.24104]'

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
        os.remove(self_region)

        if subtract:
            os.remove(bkg_region)

        hdu_list.close()

        return lc

    #method to make a source object from the region object
    def make_source(self,parent=None):
        #if parent is known, can specify it here
        return Source(lightcurve=self.make_lc(),
                        obsid=self.obsid,
                        position=f'{self.ra}{self.dec}',
                        region = self.regtext,
                        region_object = self,
                        parent = parent)

    def make_new_source(self,all_dict):
        '''called when no suitable match is found and a new source must be made'''
        new_source = self.make_source()

        new_source_all = Source_All([new_source],self.ra,self.dec)

        new_source.parent = new_source_all

        all_dict[new_source_all] = [self]

        return


    def apply_match(self, all_dict):
        '''Called when a match between regions is found to handle the matching process'''

        #the first time this is called, the closest region will be the match
        if self.best_match is None:
            self.matches.sort(key = lambda x: x.sep)
            closest_match = self.matches[0]

            self.best_match = closest_match

        else:
            current_match_index = self.matches.index(self.best_match)

            try:
                closest_match = self.matches[current_match_index + 1]
            except IndexError:
                for arr in all_dict.values():
                    try:
                        assert self not in arr
                    except AssertionError:
                        arr.remove(self)

                print('No suitable match, making new source.')
                self.make_new_source(all_dict)

                self.best_match = None

                return

            self.best_match = closest_match

        obsids_in_closest_match = [x.obsid for x in closest_match.source.obs]

        if self.obsid in obsids_in_closest_match:
            #test the seperation of both matches:
            #alt_region is the region which was previously matched to
            #the source which shares an obsid with self

            alt_source = [x for x in closest_match.source.obs if x.obsid == self.obsid][0]

            alt_region = alt_source.region_object

            if self.best_match.sep < alt_region.best_match.sep:
                closest_match.source.obs.remove(alt_source)

                assert alt_source not in closest_match.source.obs

                closest_match.source.obs.append(self.make_source(parent = closest_match.source))

                all_dict[self.best_match.source].append(self)

                #have to re match the alt source
                alt_region.apply_match(all_dict)

            else:
                self.apply_match(all_dict)

        else:
            #make the region into a source and add it to the source_all object
            closest_match.source.obs.append(self.make_source(parent = closest_match.source))
            all_dict[self.best_match.source].append(self)

        return

def unglob(arr,force=False):
    #force is used to force the program to take the first option of multiple files
    #and to not ask the user whether or not to continue
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
    evt = unglob(evt,True)

    regphystocel.punlearn()
    regphystocel.infile = region_file
    regphystocel.outfile = f'{region_dir}/wavdetect_wcs_reg.fits'
    regphystocel.wcsfile = evt
    regphystocel.clobber = 'yes'
    regphystocel()

    regiontxt = np.loadtxt(f'{region_dir}/wavdetect_wcs_reg.fits',dtype='str',skiprows=3,ndmin=1)

    if len(regiontxt) == 0:
        return None

    def make_dec(line):
        dec = line[7:].strip('()').split(',')[1]

        #the default for declination is positive,
        #in which case we have to add a plus sign for pyasl
        if '-' not in dec and not '+' in dec:
            dec = f'+{dec}'

        dec_pieces = dec.split(':')

        piece_arr = []
        for i,p in enumerate(dec_pieces):
            if i == 0:
                if len(p) != 3:
                    p = p[0] + '0' + p[1:]
            else:
                if len(p) != 2 and len(p) != 7:
                    p = '0' + p

            piece_arr.append(p)

        return ':'.join(piece_arr)

    def make_ra(line):
        ra = line[7:].strip('()').split(',')[0]

        ra_pieces = ra.split(':')

        piece_arr = []
        for p in ra_pieces:
            if len(p) != 2 and len(p) != 7:
                p = '0' + p

            piece_arr.append(p)

        return ':'.join(piece_arr)


    regions = [Region(region_file, #path
                line, #reg text
                make_ra(line), #ra
                make_dec(line), #dec
                max(line[7:].strip('()').split(',')[2],line[7:].strip('()').split(',')[3]).rstrip('"'), #a
                min(line[7:].strip('()').split(',')[2], line[7:].strip('()').split(',')[3]).rstrip('"'), #b
                obsid, #obsid
                evt, #evt
                f'{region_dir}/wavdetect_wcs_reg.fits') #all_reg_path
                for line in regiontxt]

    out = Obsid_all_regions(obsid,regions)

    return out

if __name__ == '__main__':
    region_dir = sys.argv[1]

    obsid = sys.argv[2]

    regions = process_wavdetect(obsid,region_dir,f'{region_dir}/detect_src.fits').all_regions


    for reg in regions:
        print(reg.ra)
        print(reg.dec)
        print(reg.regtext)
        print(reg.name)
        print('\n')
