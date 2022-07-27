import numpy as np
import os
from ciao_contrib.runtool import *

#returns the area of the given region in the given events file at the given energy
def region_area(events,region,energy):
    dmlist.punlearn()
    dmlist.infile = f'{events}[energy={energy}][sky=region({region})]'
    dmlist.opt = 'subspace,clean'
    dmlist.cells = '1:2'

    dmlist_out = dmlist()

    for line in dmlist_out.splitlines():
        if 'Region area' in line:
            area_line = line

    area = float(area_line[area_line.find('Region area'):].strip('Region area = ').rstrip())
    return area

#takes as arguments arrays of soft/hard source and background counts, the area ratios, and some file keeping
def make_behr_file(BEHR_DIR,soft_src,hard_src,soft_bkg,hard_bkg,soft_area,
                    hard_area,outfile,BEHR_outfile):
    '''
    if os.path.exists(outfile):
        cont = ''
        while 'y' not in cont and 'n' not in cont:
            cont = input('BEHR outfile exists. Proceeding would overwrite previous work. \nContinue? (y/n): ')
            if 'n' in cont:
                raise Exception
    '''

    with open(outfile,'w') as writeto:
        writeto.write(f'cd {BEHR_DIR}')
        for i in range(len(soft_src)):
            writeto.write(f'\n./BEHR softsrc={soft_src[i]} hardsrc={hard_src[i]}   softbkg={soft_bkg[i]}   hardbkg={hard_bkg[i]}   softarea={soft_area} hardarea={hard_area} output={BEHR_outfile}/{i}_BEHRresults')
