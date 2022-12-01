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

    if override:
        src_region = input('Enter path to src region file: ').strip()
        bkg_region = input('Enter path to bkg region file: ').strip()

    else:
        print('Running srcflux to make regions...')
        #If not overriden, we make the regions with src flux
        make_regions(obsid,position,f'{working_dir}/')

        src_region = unglob(glob.glob(f'{working_dir}/*srcreg.fits'),True)
        bkg_region = unglob(glob.glob(f'{working_dir}/*bkgreg.fits'),True)

    evt = unglob(glob.glob(f'{primary}/*evt2*'))

    #check to see if the background is wrong
    #when the srcflux source region extends over the edge of a chip, the
    #corresponding background region will be far far too large
    #I am going to detect this by comparing the area of each

    try:

        area = region_area(evt,bkg_region,1000)/region_area(evt,src_region,1000)
    except Exception as e:
        print('If OSError because the virtual file for a region could not be opened')
        print('Remember to save the bkg region in ciao format instead of ds9')
        print('Talk me if you need help')

        raise e

    try:
        assert area < 50
    except AssertionError as e:
        print('The background region is wrong, please run the program with "bkgoverride" at the end and follow the instruction')
        raise e


    print('Making BEHR bash file for 68%...')

    outfile = f'{working_dir}/BEHR_bash.txt'

    BEHR_outdir = f'{BEHR_DIR}/{obsid}/{position_basic}'

    try:
        subprocess.run(f'rm -rf {BEHR_outdir}',shell=True)
    except:
        pass
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
    print('Remember to change the BEHR_DIR paramter located in this file to where your BEHR instillation is located.')

    obsid = sys.argv[1]
    position = sys.argv[2]
    N = int(sys.argv[3])

    #can manually enter an energy to divide hard and soft as the last command
    #line paramter
    if len(sys.argv) >= 5:
        try:
            divide = float(sys.argv[4])
        except:
            divide = 2000
    else:
        divide = 2000

    override = input('Override regions (be asked to enter your own regions instead of making them automatically)  [y/n]? ')
    if 'y' in override:
        override = True
    else:
        override = False

    lines = None

    main(obsid,position,N,lines,override=override,divide_energy=divide)
