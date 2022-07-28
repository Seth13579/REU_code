from run_wavdetect import *
from source_class import Source_All, Source
import numpy as np
import os
import sys
import astropy.units as u
from astropy.coordinates import SkyCoord

#when given two regions, returns their seperation in arcsec
def seperation(reg1,reg2):

    ra1 = reg1.ra
    dec1 = reg1.dec
    ra2 = reg2.ra
    dec2 = reg2.dec


    c1 = SkyCoord(ra1,dec1,unit=(u.hourangle, u.deg),frame='fk5')
    c2 = SkyCoord(ra2,dec2,unit=(u.hourangle, u.deg),frame='fk5')

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


class Match:
    '''Describing the match between a region and a source'''
    def __init__(self,region,matched_region,source):
        #the region which is being match to a source
        self.region = region

        #the region to which the region matched to
        self.matched_region = matched_region

        #the source to which the region is matched
        self.source = source

        self.sep = seperation(region,matched_region)


class Galaxy:
    '''class representing all the observations of one galaxy'''

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

        #error_rate = 0
        #clean_rate = 0
        for i,region1 in enumerate(mast_region_list):

            #an array of Match objects which represent all the matches
            #the source got
            match_array = []

            print(f'\nLooking to match a region in ObsID {region1.obsid}, {i+1} of {len(mast_region_list)}:')
            matched = False

            #a debugging tool to learn how many regions are being compared against
            z = 0
            for source in all_dict.keys():
                #he list of region objects which correspond to the regions
                #which make up source
                test_against = all_dict[source]

                for region2 in test_against:
                    z += 1

                    if match_test(region1,region2):
                        match_array.append(Match(region1,region2,source))

                        matched = True


            #after we check all the sources, if no match, then we have a new source
            if not matched:

                print('No match, making new source...')

                region1.make_new_source(all_dict)

            #if a match is made, append to the source which is the closest match
            else:
                region1.matches = match_array

                print(f"{len(region1.matches)} matches found, finding closest and appending...")

                region1.apply_match(all_dict)



            print(f'Region was tested against {z} other regions in total')



        return all_dict.keys()
