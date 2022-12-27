from Source import *
import gc

#example textfile name: J051952.2412-692553.9081_14030_lc.fits

def filename_coords(filename):
    filename_list = filename.split('/')[-1].split('_')

    coords = filename_list[0].strip("J")

    return coords

def filename_obsid(filename):
    filename_list = filename.split('/')[-1].split('_')

    obsid = filename_list[1]

    return obsid

def find_same_source(file_list):
    out_dict = {}

    for i,file in enumerate(file_list):
        coords = filename_coords(file)

        if coords not in out_dict.keys():
            out_dict[coords] = [i]
        else:
            out_dict[coords].append(i)

    return out_dict

def main():
    galaxy_dirs = sys.argv[1].split(',')

    if 'all' in galaxy_dirs:
        all = True
        galaxy_dirs = [i for i in os.listdir() if '.' not in i and i != 'Detected']
    else:
        all = False

    cwd = os.getcwd()

    for j,galaxy_dir in enumerate(galaxy_dirs):
        print(f'''***********
ANALYZING {galaxy_dir}, {j+1} OF {len(galaxy_dirs)}
***********''')

        os.chdir(galaxy_dir)

        if all:
            if os.path.exists('./detections_lcs'):
                os.chdir(cwd)
                continue

        try:
            os.makedirs('detections_lcs')
        except:
            pass
        os.chdir('textfiles')

        textfiles = glob.glob('*_lc.fits.txt')

        for i,file in enumerate(textfiles):
            print(f'Testing {file}, {i+1} of {len(textfiles)}')

            coords = filename_coords(file)
            obsid = filename_obsid(file)

            lc = np.loadtxt(file,skiprows=1)

            src = Source(lightcurve = lc, obsid = obsid, position = coords)

            src.e_test_blocks(binsize = 500)

            if src.classification:
                print('Match found!')
                src.make_fourpanel_plot(outdir='../detections_lcs')

            del src

            #gc.collect()

        os.chdir(cwd)

if __name__ == '__main__':
    main()
