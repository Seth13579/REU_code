from run_wavdetect import *
from source_class import Source_All, Source
import numpy as np
import os
import sys
import astropy.units as u
from astropy.coordinates import SkyCoord
from PyAstronomy import pyasl

#when given two regions, returns their seperation in arcsec
def seperation(reg1,reg2):
    #first convert sexigesimal coordinates to degrees
    coord_str_1 = " ".join([reg1.ra,reg1.dec])
    coord_str_2 = " ".join([reg2.ra,reg2.dec])

    try:
        ra1,dec1 = pyasl.coordsSexaToDeg(coord_str_1)
    except Exception as e:
        print(reg1.ra)
        print(reg1.dec)
        print(reg1.path)

        raise e
    try:
        ra2,dec2 = pyasl.coordsSexaToDeg(coord_str_2)
    except Exception as e:
        print(reg2.ra)
        print(reg2.dec)
        print(reg2.path)

        raise e




    c1 = SkyCoord(ra1,dec1,unit=u.deg,frame='fk5')
    c2 = SkyCoord(ra2,dec2,unit=u.deg,frame='fk5')

    return c1.separation(c2).arcsecond

'''
When given two regions, tests if they represent the same source
The criteria of which is as follows:
    If the seperation betweem the regions is less than 2 * either of the
    semimajor axes of the regions
'''
def match_test(reg1,reg2):
    if reg1.obsid == reg2.obsid:
        return False

    sep = seperation(reg1,reg2)

    try:
        if sep <= 2*float(reg1.a) or sep <= 2*float(reg2.a):
            return True
        else:
            return False
    except TypeError as e:
        print(sep)
        print(reg1.a)
        print(reg2.a)

        raise e


#class representing all the observations of one galaxy
class Galaxy:

    def __init__(self,name,obsid_region_list):
        #the name of the galaxy
        self.name = name

        #the list of "Obsid_all_regions" objects
        #represents the output of wavdetect for all the obsids in the galaxy
        self.obsid_region_list = obsid_region_list

        #a list of the obsids in the galaxy
        self.obsids = [i.obsid for i in obsid_region_list]


    #perfroms region matching and returns a list of Source_All objects which
    #represents all the distinct sources in the galaxy
    def region_match(self):
        #list of type Obsid_all_regions
        obsid_list = self.obsid_region_list

        mast_region_list = [reg for ob in obsid_list for reg in ob.all_regions]

        '''
        A dictionary which represents the sources and the regions which have
        been matched to that source.

        Each key is a Source_All object

        Each value is a list of region objects corresponding to the regions
        which were matched to the key
        '''
        all_dict = {}

        for i,region1 in enumerate(mast_region_list):

            print(f'\nLooking to match {region1.path}, {i+1} of {len(mast_region_list)}:')
            matched = False

            for source in all_dict.keys():
                test_against = all_dict[source]

                for region2 in test_against:

                    if match_test(region1,region2):
                        #add the region to the dictionary so we can later check against it
                        all_dict[source].append(region1)

                        #make the region into a source and add it to the source_all object
                        source.obs.append(region1.make_source())

                        matched = True

                        print('Match!')

                        break

                if matched:
                    break

            #after we check all the sources, if no match, then we have a new source
            if not matched:

                print('No match, making new source...')

                new_source = region1.make_source()

                new_source_all = Source_All([new_source],region1.ra,region1.dec)

                all_dict[new_source_all] = [region1]

        return all_dict.keys()
