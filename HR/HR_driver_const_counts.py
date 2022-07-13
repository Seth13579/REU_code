import numpy as np
import sys
import os
from extract_counts import *
from BEHR_countbins import *
from run_BEHR import *
from re import sub
from create_behr_files import region_area

def main(obsid,position,N,lines,divide_energy=2000):
    #BEHR_DIR = '/data/reu/slarner/BEHR_contain/BEHR'
    #or use:
    BEHR_DIR = '/Users/sethlarner/BEHR_contain/BEHR'

    download_obsid(obsid)

    primary = f'./{obsid}/primary'

    #remove the colons from the position
    position_basic = sub('\:','',position)
    working_dir = f'{primary}/{position_basic}'

    print('Running srcflux to make regions...')
    make_regions(obsid,position,f'{working_dir}/')

    src_region = unglob(glob.glob(f'{working_dir}/*srcreg.fits'))
    bkg_region = unglob(glob.glob(f'{working_dir}/*bkgreg.fits'))
    evt = unglob(glob.glob(f'{primary}/*evt2*'))

    #check to see if the background is wrong
    #when the srcflux source region extends over the edge of a chip, the
    #corresponding background region will be far far too large
    #I am going to detect this by comparing the area of each

    area = region_area(evt,bkg_region,1000)/region_area(evt,src_region,1000)

    try:
        assert area < 100
    except AssertionError as e:
        print('The background region is wrong, please draw a suitible background by hand and run the program with ')


    print('Making BEHR bash file...')

    outfile = f'{working_dir}/BEHR_bash.txt'

    BEHR_outdir = f'{BEHR_DIR}/{obsid}/{position_basic}'

    subprocess.run(f'rm -rf {BEHR_outdir}',shell=True)
    os.makedirs(BEHR_outdir)

    make_behr(evt,src_region,bkg_region,divide_energy,BEHR_DIR,outfile,BEHR_outdir,N)

    print('Running BEHR...')

    run_BEHR(outfile)

    print('Making plot...')
    time,uppers,lowers,meds = plot_BEHR_constcounts(BEHR_outdir,N,position_basic,obsid,evt,src_region,show = True,lines=lines)


    print('Saving...')
    to_save = np.column_stack((time,meds,uppers,lowers))

    np.savetxt(f'./HR_saved_{position}_{obsid}.txt',to_save,delimiter=',',header='Time,Median,Upper Error,Lower Error')



if __name__ == '__main__':
    obsid = sys.argv[1]
    position = sys.argv[2]
    N = int(sys.argv[3])

    if len(sys.argv) > 4:
        bkg_override = sys.argv[4]
        if bkg_override == '!bkg_override':
            override = True
        else:
            override = False 

    try:
        lines = np.array(sys.argv[4].split(',')).astype('float64')
    except (IndexError,TypeError,ValueError):
        lines = None

    main(obsid,position,N,lines)
