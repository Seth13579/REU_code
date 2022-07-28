import numpy as np
import sys
import os
from extract_counts import *
from BEHR_countbins import *
from run_BEHR import *
from re import sub
from create_behr_files import region_area

def main(obsid,position,N,lines,divide_energy=2000,override=False,subtract_start=True):
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

    src_region = unglob(glob.glob(f'{working_dir}/*srcreg.fits'),True)

    if override:
        print('To avoid background issue, you must draw by hand a background region')
        print('Please do that now and then input the absolute path to that region below')
        print('Make sure to save the region as a fits')
        bkg_region = input('Insert path to hand-drawn region: ')
    else:
        bkg_region = unglob(glob.glob(f'{working_dir}/*bkgreg.fits'),True)

    evt = unglob(glob.glob(f'{primary}/*evt2*'))

    #check to see if the background is wrong
    #when the srcflux source region extends over the edge of a chip, the
    #corresponding background region will be far far too large
    #I am going to detect this by comparing the area of each

    area = region_area(evt,bkg_region,1000)/region_area(evt,src_region,1000)

    try:
        assert area < 50
    except AssertionError as e:
        print('The background region is wrong, please run the program with "bkgoverride" at the end and follow the instruction')
        raise e


    print('Making BEHR bash file for 68%...')

    outfile = f'{working_dir}/BEHR_bash.txt'

    BEHR_outdir = f'{BEHR_DIR}/{obsid}/{position_basic}'

    subprocess.run(f'rm -rf {BEHR_outdir}',shell=True)
    os.makedirs(BEHR_outdir)

    make_behr(evt,src_region,bkg_region,divide_energy,BEHR_DIR,outfile,BEHR_outdir,N)

    print('Running BEHR for 68%...')

    run_BEHR(outfile)

    if subtract_start:
        dmextract.punlearn()
        dmextract.infile = f'{evt}[bin time=::3.24014]'
        dmextract.clobber = 'yes'
        dmextract.opt = 'generic'
        dmextract.outfile = 'temp.fits'
        dmextract()

        dmlist.punlearn()
        dmlist.infile = 'temp.fits[cols time]'
        dmlist.opt = 'data, clean'

        start_time = float(dmlist().splitlines()[1])
    else:
        start_time = 0

    time_68,uppers_68,lowers_68,meds_68 = plot_BEHR_constcounts(BEHR_outdir,N,position_basic,obsid,evt,src_region,lines=lines,start_time=start_time,show=False,save=False)

    print('Making BEHR bash file for 90%...')

    subprocess.run(f'rm -rf {BEHR_outdir}',shell=True)
    os.makedirs(BEHR_outdir)

    make_behr(evt,src_region,bkg_region,divide_energy,BEHR_DIR,outfile,BEHR_outdir,N,confidence='90.00')

    print('Running BEHR for 90%...')

    run_BEHR(outfile)

    time_90,uppers_90,lowers_90,meds_90 = plot_BEHR_constcounts(BEHR_outdir,N,position_basic,obsid,evt,src_region,lines=lines,start_time=start_time,show=False,save=False)

    print('Plotting...')

    plt.plot(time_68,meds_68,'k-')
    plt.fill_between(time_90,lowers_90,uppers_90,step='mid',color='g')
    plt.fill_between(time_68,lowers_68,uppers_68,step='mid',color='b')    

    plt.ylabel('(H-S)/(H+S)')
    plt.xlabel('Time (s)')
    plt.title(f'HR J{position} ObsID: {obsid} \n Each point ± {N} counts')

    plt.savefig(f'{position}_{obsid}_HR_constcounts_68and90.png',dpi=300)

    plt.show()

    print('Saving...')
    to_save_68 = np.column_stack((time_68,meds_68,uppers_68,lowers_68))
    np.savetxt(f'./HR_saved_{position}_{obsid}_68.txt',to_save_68,delimiter=',',header='Time,Median,Upper Error,Lower Error')

    to_save_90 = np.column_stack((time_90,meds_90,uppers_90,lowers_90))
    np.savetxt(f'./HR_saved_{position}_{obsid}_90.txt',to_save_90,delimiter=',',header='Time,Median,Upper Error,Lower Error')


    return (to_save_68,to_save_90)

if __name__ == '__main__':
    obsid = sys.argv[1]
    position = sys.argv[2]
    N = int(sys.argv[3])

    if len(sys.argv) >= 5:
        try:
            divide = float(sys.argv[4])
        except:
            divide = 2000
    else:
        divide = 2000


    '''
    try:
        lines = np.array(sys.argv[4].split(',')).astype('float64')
        if 'bkgoverride' in lines:
            lines = None
    except (IndexError,TypeError,ValueError):
        lines = None
    '''

    lines = None

    if 'bkgoverride' in sys.argv:
        main(obsid,position,N,lines,override=True,divide_energy=divide)
    else:
        main(obsid,position,N,lines,override=False,divide_energy=divide)
