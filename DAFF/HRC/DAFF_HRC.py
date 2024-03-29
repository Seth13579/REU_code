#The beginning of the full process
#Dip And Flare Finder

from Source import Source, Source_All
import os
from region_matching import *
import re
from run_wavdetect_HRC import *


#function to handle computationally expensive steps
def subreprocess(dir):
    obsid = dir.split('/')[-1]

    working_dir = f'./{dir}/repro'

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

        try:
            new = subreprocess(dir)
        except:
            new = None

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
    total = len([source for all_source in all_sources_in_galaxy for source in all_source.obs])
    i = 0
    for all_source in all_sources_in_galaxy:
        for source in all_source.obs:

            try:
                source.t_interest = source.e_test_repeat(binsize=1000)
            except:
                source.t_interest = None
                source.classification = False

            if source.classification:
                source.make_fourpanel_plot(outdir='.')

            csv_line = source.make_csv_line()

            lines.append(csv_line)

            i += 1

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

                try:
                    source.make_HR(16,evt_dir)
                    source.plot_HR_and_lc(1000,'.')
                except:
                    pass

    ###########################
    ###########################
    print('Saving sources...')
    os.chdir('../')
    try:
        os.makedirs('./sources')
    except:
        pass
    os.system('rm ./sources/*.pkl')
    os.chdir('./sources')

    for all_source in all_sources_in_galaxy:
        all_source.save()

    return

if __name__ == '__main__':
    #first command line argument controls which galaxies to run on
    #can either list galaxies, seperated by commas
    #or use 'all' to run on all dirs in cwd

    if 'all' in sys.argv[1] or 'ALL' in sys.argv[1] or 'All' in sys.argv[1]:
        galaxies = [i for i in os.listdir(os.getcwd()) if '.txt' not in i and i != 'NGC5194']

    else:
        galaxies = sys.argv[1].split(',')

    errors = []

    cwd = os.getcwd()

    try:
        os.makedirs('./errored')
    except:
        pass

    try:
        os.makedirs('./completed')
    except:
        pass

    for i,galaxy in enumerate(galaxies):
        print('***********')
        print(f'PROCESSING {galaxy}, {i+1} OF {len(galaxies)}')
        print('***********')

        try:
            process_galaxy(galaxy)
        except Exception as e:
            errors.append(galaxy)
            os.chdir(cwd)
            os.system(f'mv {galaxy} ./errored')
            raise e
            print(f'{galaxy} erroed, moved to ./errored')
        else:
            os.chdir(cwd)
            os.system(f'mv {galaxy} ./completed')
            print(f'{galaxy} completed without error, moved to ./completed')

    with open('Error_doc.txt','w') as f:
        for gal in errors:
            f.write(gal)
