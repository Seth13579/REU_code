import sys
import os
from e_test_testing import *

def main():

    file = '/Volumes/External/REU_main/t_test/all_sk_galaxies/NGC300/textfiles/J005501.00-373440.43_12238_lc.fits.txt'

    if 'ascii' in file:
        skip = 0
    else:
        skip = 1

    source = Source(lightcurve = np.loadtxt(file, skiprows = skip))


    binned_time,binned_counts = source.make_binned_counts(2000)



    analyze_file(file,True,500,3,'test')

    os.system('open /Volumes/External/REU_main/t_test/test/detections/005501.00-373440.43_12238_full_lc.pdf')



if __name__ == '__main__':


    main()
