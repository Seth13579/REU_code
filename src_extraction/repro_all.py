#for reprocessing a bunch of galaxies
import pickle as pkl
import os
from galaxy_dwnld import *

if __name__ == '__main__':
    with open('observed_galaxies_dict','rb') as f:
        dict = pkl.load(f)

    cwd = os.getcwd()


    i = 0
    for galaxy, obsids in dict.items():
        i += 1
        print(f'Processing {galaxy}, {i} of {len(dict.keys())}')
        os.makedirs(f'./{galaxy}')
        os.chdir(f'./{galaxy}')



        for obsid in obsids:
            obsid_download(obsid)

        split_interleaved()

        repro()

        os.chdir(cwd)
