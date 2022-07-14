import numpy as np
import sys
import os
from extract_counts import *
from create_behr_files import *
from run_BEHR import *
from re import sub


def create_hr_plot(obsid,position,bin_size,lines):

    #change to match where BEHR is contained on your machine
    BEHR_DIR = '/data/reu/slarner/BEHR_contain/BEHR'
    #BEHR_DIR = '/Users/sethlarner/BEHR_contain/BEHR'
    download_obsid(obsid)

    primary = f'./{obsid}/primary'

    #remove the colons from the position
    position_basic = sub('\:','',position)
    working_dir = f'{primary}/{position_basic}'

    print('Running srcflux to make regions...')
    make_regions(obsid,position,f'{working_dir}/')

    #split the source and background counts by time and energy
    #should take care that the old region files are being overwirtten or delted between runs
    src_region =unglob(glob.glob(f'{working_dir}/*srcreg.fits'))
    bkg_region = unglob(glob.glob(f'{working_dir}/*bkgreg.fits'))
    evt = unglob(glob.glob(f'{primary}/*evt2*'))

    print('Binning and splitting events...')
    soft_src = split_events(evt,src_region,bin_size,'soft')
    hard_src = split_events(evt,src_region,bin_size,'hard')
    soft_bkg = split_events(evt,bkg_region,bin_size,'soft')
    hard_bkg = split_events(evt,bkg_region,bin_size,'hard')

    print('Calculating areas...')
    hard_area_ratio = region_area(evt,bkg_region,3000)/region_area(evt,src_region,3000)
    soft_area_ratio = region_area(evt,bkg_region,1000)/region_area(evt,src_region,1000)

    outfile = f'{working_dir}/BEHR_bash.txt'

    BEHR_outdir = f'{BEHR_DIR}/{obsid}/{position_basic}'

    subprocess.run(f'rm -rf {BEHR_outdir}',shell=True)
    os.makedirs(BEHR_outdir)

    make_behr_file(BEHR_DIR,soft_src,hard_src,soft_bkg,hard_bkg,soft_area_ratio,hard_area_ratio,outfile,BEHR_outdir)

    print('Running BEHR...')

    run_BEHR(outfile)

    make_hr_plot(BEHR_outdir,lc_dir,bin_size,position_basic,obsid,show = True,lines=lines)

if __name__ == '__main__':

    obsid = sys.argv[1]
    position = sys.argv[2]
    bin_size = sys.argv[3]
    lc_dir = sys.argv[4]

    try:
        lines = np.array(sys.argv[5].split(',')).astype('float64')
    except (IndexError,TypeError,ValueError):
        lines = None

    create_hr_plot(obsid,position,bin_size,lines)
