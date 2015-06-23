#!/usr/bin/env python

#Author: Maurice Wilson
'''  
      Purpose
Search for satisfactory comparison stars.

Input resultant stars into final "target list" file.


      Command Line Call
./query_stars.py eta_merged_list.txt
OR
>>> import sys
>>> sys.argv = ['', 'eta_merged_list.txt']
>>> execfile('query_stars.py')

Note: make sure that any initial target list of stars, such as the 'eta_merged_list.txt' file, maintains that tabulated (i.e., delimiter= \t) table format.
'''


import sys

from astropy.table import Table, Column
import numpy as np

from astroquery.simbad import Simbad

import json, collections, datetime

infile = sys.argv[1]

def search_stars(infile):
    with open(infile, 'r') as p_star_file:
        file_strings = p_star_file.readlines()
    #p_star_file.close() # this is unnecessary

    # find which line in the file starts the table
    finisher = 1
    list_index = 0
    while finisher == 1:    
        check = file_strings[list_index]
        if check[0:4] == '0\tHD':
            finisher = 0
            table_start_index = list_index
        list_index += 1

    # make each line in file one string array
    p_star_table_strings = [file_strings[ind].replace('\t\t', '\t').rstrip('\n') for ind in xrange(table_start_index, len(file_strings))]

    # in each string array make an array of strings divided by the columns
    p_star_table_arrays = [p_star_table_strings[ind].split('\t') for ind in xrange(len(p_star_table_strings))]
    
    # separate information based on the columns
    p_star_indice = [int(p_star_table_arrays[ind][0]) for ind in xrange(len(p_star_table_arrays))]
    
    p_star_id = [p_star_table_arrays[ind][1] for ind in xrange(len(p_star_table_arrays))]

    p_star_spec = [p_star_table_arrays[ind][2] for ind in xrange(len(p_star_table_arrays))]

    p_star_mass = [float(p_star_table_arrays[ind][3]) for ind in xrange(len(p_star_table_arrays))]

    p_star_RA = [float(p_star_table_arrays[ind][4]) for ind in xrange(len(p_star_table_arrays))]

    p_star_Dec = [float(p_star_table_arrays[ind][5]) for ind in xrange(len(p_star_table_arrays))]

    p_star_Vmag = [float(p_star_table_arrays[ind][6]) for ind in xrange(len(p_star_table_arrays))]

    # make a table of the primary target stars' information
    p_star_table = Table([p_star_indice, p_star_id, p_star_spec, p_star_mass, p_star_RA, p_star_Dec, p_star_Vmag], names=('Indice', 'ID', 'Spec Type', 'Mass (solar)', 'RA (deg)', 'Dec (deg)', 'V mag'))

    #print(p_star_table)



    # customize the output of the query results
    customS = Simbad()
    customS.remove_votable_fields('coordinates')
    customS.add_votable_fields('coo(d)')
    customS.add_votable_fields('flux(B)')
    customS.add_votable_fields('flux(V)')
    customS.add_votable_fields('otype(V)')
    customS.ROW_LIMIT = 4  # search for only 4 results/comparison stars
    customS.TIMEOUT = 100

    mag_limit = 1.5 # the range at which a comparison star's magnitude must be within

    # obtain the magnitude of the primary target stars
    mags_V = [p_star_table['V mag'][p_star] for p_star in xrange(len(p_star_table['ID']))] # p_star is a number

    # find upper and lower limit for each star's magnitude
    mags_V_p = np.array(mags_V)+mag_limit
    mags_V_m = np.array(mags_V)-mag_limit

    # determine B - V color of primary star
    # add/subtract to find upper and lower limits of color
    #PROBABLY UNNECESSARY
                  
    # criteria-query for comparison stars based on the parameters: V magnitude, vicinity, and if it is not variable, (and maybe color)
    def find_comparison(primary_S_indice):
        return customS.query_criteria('region(circle, '+p_star_table['ID'][primary_S_indice]+', 5d)', 'maintypes=*', 'maintypes!=V*', 'Vmag < '+str(mags_V_p[primary_S_indice]), 'Vmag > '+str(mags_V_m[primary_S_indice]))

    decision = raw_input('Do you want to search for a comparison star for all '+str(len(p_star_table['ID']))+' primary target stars?  (Type 1)\nOr do you want to only search around 1 primary target star for a comparison star (Type 2)\n---> ')

                   
    if decision == '1':
        comparison_S_tables =  [customS.query_criteria('region(circle, '+p_star_table['ID'][star]+', 5d)', 'maintypes=*', 'maintypes!=V*', 'Vmag < '+str(mags_V_p[star]), 'Vmag > '+str(mags_V_m[star])) for star in xrange(len(p_star_table['ID']))]
    elif decision == '2':
        primary_S_indice  = input('What is the indice of the primary star within the '+sys.argv[1]+' file?\n---> ')
        comparison_table = find_comparison(primary_S_indice)
        print(comparison_table)
    else:
        print('You have typed an invalid number.  Crash program...')
        return

    # input the primary and comparison stars into file formatted at a readable Target List


    telnum = 1

    tonight = datetime.date.today() + datetime.timedelta(days=1)
    # specify the calibrations
    filename = datetime.datetime.strftime(tonight,'n%Y%m%d') + '.T' + str(telnum) + '.txt'
    with open(filename, "w") as outfile:
        outfile.write('{"nbias": 11, "ndark": 11, "nflat": 9, "darkexptime": [300], "flatFilters": ["ip","R"]\, "WaitForMorning": true}\n')
        outfile.write('{"nbiasEnd": 11, "ndarkEnd": 11, "nflatEnd": 9}\n')


    if telnum == 1:
        name = p_star_table['ID'][0]
        rahrs = p_star_table['RA (deg)'][0]
        decdegs = p_star_table['Dec (deg)'][0]
    elif telnum == 2:
        name = comparison_S_tables[0]['MAIN_ID'][0]
        rahrs = comparison_S_tables[0]['RA_d'][0]
        decdegs = comparison_S_tables[0]['DEC_d'][0]

    
    # Add a target to the file:
    target = collections.OrderedDict()
    target['name'] = name
    target['ra'] = rahrs
    target['dec'] = decdegs
    target['starttime'] = datetime.datetime.strftime(starttime.datetime(),'%Y-%m-%d %H:%M:%S')
    target['endtime'] = datetime.datetime.strftime(endtime.datetime(),'%Y-%m-%d %H:%M:%S')
    target['filter'] = ["ip"]
    target['exptime'] = [300]
    target['num'] = [1]
    target['defocus'] = 0.0
    target['selfguide'] = True
    target['guide'] = False
    target['cycleFilter'] = True
    target['positionAngle'] = 0.0
    with open(filename, "a") as outfile:
        json.dump(target,outfile)
        outfile.write('\n')





#search_stars(infile)
        
'''    #Scraps

# customize the output of the query results
customS = Simbad()
customS.remove_votable_fields('coordinates')
customS.add_votable_fields('flux(B)')
customS.add_votable_fields('flux(V)')
customS.add_votable_fields('coo(s)')
customS.ROW_LIMIT = 2

primaries = ["rho Cnc", "sirius"] # save the list of primary target stars to variable
mag_limit = 3
  # loop for each star
  # query for primary star  
primary_S_tables = [customS.query_object(primaries[star]) for star in xrange(len(primaries))] # star is a number
#print(primary_S_tables)

  # add/subtract to find upper and lower limits to the magnitudes B and V
mags_B = [primary_S_tables[p_star]['FLUX_B'][0] for p_star in xrange(len(primary_S_tables))] # p_star is a number
mags_B_p = np.array(mags_B)+mag_limit
mags_B_m = np.array(mags_B)-mag_limit

mags_V = [primary_S_tables[p_star]['FLUX_V'][0] for p_star in xrange(len(primary_S_tables))] # p_star is a number
mags_B_p = np.array(mags_B)+mag_limit
mags_B_m = np.array(mags_B)-mag_limit

  # determine B - V color of primary star
    # add/subtract to find upper and lower limits of color
                  #PROBABLY UNNECESSARY
                  
  # criteria-query for comparison stars based on the parameters: magnitudes, vicinity, (and color)
    # choose first 3 or 4 results
comparison_S_tables =  [customS.query_criteria('region(circle, '+primaries[star]+', 1d)', 'maintypes!=V*', 'Bmag < '+str(mags_B_p[star]), 'Bmag > '+str(mags_B_m[star])) for star in xrange(len(primaries))]

#comparison_S_tables =  [customS.query_criteria('region(circle, '+primaries[star]+', 5d)') for star in xrange(len(primaries))]

print(comparison_S_tables)



'''
