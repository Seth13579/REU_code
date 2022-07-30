#The beginning of the full process
#Dip And Flare Finder

from Source import Source, Source_All
import os
from run_wavdetect import *
from region_matching import *
import re


#function to handle computationally expensive steps
def subreprocess(dir):
    obsid = dir.split('/')[-1]

    working_dir = f'./{dir}/repro'

    # TODO: Multithread this
    try:
        print("Looking for wavdetect product...")
        wavdetect_region = unglob(glob.glob(f'{working_dir}/*detect_src.reg*'))
    except:
        print('Trying wavdetect...')
        detect(working_dir)
        wavdetect_region = unglob(glob.glob(f'{working_dir}/*detect_src.reg*'))
    else:
        print('Found!')
    finally:
        all_regions_in_obsid = process_wavdetect(obsid,working_dir,wavdetect_region)

    return all_regions_in_obsid




#carrys out the steps of region making and matching to create a galaxy class object
#assumes that the data is downloaded and reprocessed
#then produces textfiles in ./{galaxy}/textfiles
#then produces plots in
def process_galaxy(galaxy_name):
    os.chdir(f'./{galaxy_name}')

    #now run wavdetect on each obsid
    dirs = [i for i in os.listdir(os.getcwd()) if re.search('\d',i)]

    ###########################
    ###########################
    print('\nRunning wavdetect...')

    all_regions_in_galaxy = []
    for i,dir in enumerate(dirs):
        print(f'\nRunning on {dir}, {i+1} of {len(dirs)}')

        new = subreprocess(dir)

        if new is not None:
            all_regions_in_galaxy.append(new)

    galaxy = Galaxy(galaxy_name,all_regions_in_galaxy)

    ###########################
    ###########################
    print('\nMatching regions...')
    all_sources_in_galaxy = galaxy.region_match()

    for source in all_sources_in_galaxy:
        source.optimal_rename()

    try:
        os.makedirs('./textfiles')
    except:
        pass

    os.chdir('./textfiles')
    os.system('rm *.txt')

    ###########################
    ###########################
    print('\nProducing textfiles...')
    for all_source in all_sources_in_galaxy:
        all_source.save_lcs()

    os.chdir('../')
    try:
        os.makedirs('./detections_lcs')
    except:
        pass
    os.chdir('./detections_lcs')

    ###########################
    ###########################
    print('\nDetecting dips and flares...')
    lines = []
    for all_source in all_sources_in_galaxy:
        for source in all_source.obs:
            source.t_interest = source.e_test(binsize=500)

            if source.classification:
                source.make_fourpanel_plot(outdir='.')

            csv_line = source.make_csv_line()

            lines.append(csv_line)

    ###########################
    ###########################
    print('\nSaving csv...')
    csv = np.row_stack(lines)
    header = 'NAME,OBSID,RA,DEC,TOTAL COUNTS,START TIME,END TIME,DURATION,COUNT RATE,INTERESTING,INTERESTING TIMES'

    os.chdir('../')

    np.savetxt('./all_csv.csv',csv,fmt='%s',delimiter=',',header=header)

    ###########################
    ###########################
    print('\nMaking HR plots...')
    try:
        os.makedirs('./HR_plots')
    except:
        pass

    os.chdir('./HR_plots')

    for all_source in all_sources_in_galaxy:
        for source in all_source.obs:

            if source.classification:
                evt_dir = f'../{source.obsid}/repro'

                source.make_HR(16,evt_dir)
                source.plot_HR_and_lc(1000,'.')

    return

if __name__ == '__main__':
    process_galaxy('NGC1313')
