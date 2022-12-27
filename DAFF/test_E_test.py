from DAFF_detect_only import *

def main(file):
    coords = filename_coords(file)
    obsid = filename_obsid(file)

    lc = np.loadtxt(file,skiprows=1)

    src = Source(lightcurve = lc, obsid = obsid, position = coords)

    src.e_test_blocks(binsize = 500,debug=True)

    if src.classification:
        print('Match found!')
        src.make_fourpanel_plot(outdir='.',show=True)

    del src


if __name__ == '__main__':
    file = sys.argv[1]

    main(file)
