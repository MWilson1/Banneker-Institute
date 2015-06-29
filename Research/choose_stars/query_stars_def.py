#!/usr/bin/env python

#Author: Maurice Wilson
'''  
      Purpose
Search for satisfactory comparison stars.

Input resultant stars into final "target list" file.


      Command Line Call

./query_stars.py eta_merged_list.txt

OR, in python

import sys
sys.argv = ['', 'eta_merged_list.txt']
execfile('query_stars_def.py')

Note: make sure that any initial target list of stars, such as the 'eta_merged_list.txt' file, maintains that tabulated (i.e., delimiter= \t) table format.
'''


import sys

from astropy.table import Table, Column
import numpy as np

from astroquery.simbad import Simbad

import json, collections, datetime

import ephem

infile = str(sys.argv[1])

with open(infile, 'r') as p_star_file:
    file_strings = p_star_file.readlines()
#p_star_file.close() # this is unnecessary

# find which line in the file starts the table
finisher = 1
list_index = 0
while finisher == 1:    
    check = file_strings[list_index]
    if check[0:2] == '0\t':
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



# Get information on sunset, sunrise and astronomical twilight.  Account for location of telescope.
for day_counter in xrange(14):
    
    #tonight = datetime.date.today() + datetime.timedelta(days=1)

    obs = ephem.Observer()
    obs.lat = '31.680407'
    obs.lon = '-110.878977'
    obs.horizon = '-12'    # Nautical Twilight

    '''
    today = datetime.datetime.utcnow()   # UTC TIME !!!!!        Arizona time is UTC - 7 hours             # UTC TIME !   EST is UTC - 6 hours (daylight savings time)!!!
    startNightTime = datetime.datetime(today.year, today.month, today.day, 17) - datetime.timedelta(days=1)

    tonight = datetime.date.today() + datetime.timedelta(days=1)
    '''
    
    today = datetime.datetime.utcnow() + datetime.timedelta(days= day_counter)

    startNightTime = datetime.datetime(today.year, today.month, today.day, 17) - datetime.timedelta(days=1)

    tonight = datetime.date.today() + datetime.timedelta(days= day_counter+1)  # this 'today()' is different that the variable "today"



    sunrise = obs.next_rising(ephem.Sun(), start=startNightTime,use_center=True).datetime()
    sunset = obs.next_setting(ephem.Sun(), start=startNightTime,use_center=True).datetime()

    subsequent_sunrise = obs.next_rising(ephem.Sun(), start=startNightTime + datetime.timedelta(days=1),use_center=True).datetime()
    subsequent_sunset = obs.next_setting(ephem.Sun(), start=startNightTime + datetime.timedelta(days=1),use_center=True).datetime()

    if subsequent_sunrise > subsequent_sunset:
        not_vis_latter = obs.next_rising(ephem.Sun(), start=startNightTime + datetime.timedelta(days=2),use_center=True).datetime() # latest sunrise
        not_vis_former = subsequent_sunrise
    else:
        not_vis_latter = obs.next_setting(ephem.Sun(), start=startNightTime + datetime.timedelta(days=1),use_center=True).datetime() # latest sunset
        not_vis_former = subsequent_sunset




    # Go through each primary target star and find its rise and set time.
    def choose_primary(p_star_table, p_star_indice=0):        
        name = p_star_table['ID'][p_star_indice]
        rahrs = p_star_table['RA (deg)'][p_star_indice]
        decdegs = p_star_table['Dec (deg)'][p_star_indice]
    
        # compute the rise/set times of the target
    
        obs.horizon = '30.0' 
    
        body = ephem.FixedBody()
        #body._ra = str(target['ra'])
        #body._dec = str(target['dec'])
        body._ra = str(rahrs)
        body._dec = str(decdegs)
        body._epoch = '2000.0'
        body.compute()


        observable = 'Not sure'
        time_l_visible = None
        start_time = None
        end_time = None
        midnight_t = None
    
        # the "NeverUpError" can occur so I need this try statement
        try:
            risetime = obs.next_rising(body,start=startNightTime).datetime() # start=startNightTime + datetime.timedelta(days=1)    
            settime = obs.next_setting(body,start=startNightTime).datetime()
            transit_t = obs.next_transit(body, start=startNightTime).datetime()
        
        except ephem.NeverUpError:
            print('The NeverUpError occurred for p_star_indice = '+str(p_star_indice))
            observable = False
            return observable, time_l_visible
    
        except ephem.AlwaysUpError:
            print('The AlwaysUpError occurred for p_star_indice = '+str(p_star_indice))
            observable = False
            return observable, time_l_visible

    
        def longer_than_24hrs_and_Diff(diff_of_times): # must be a datetime type
            return np.abs(diff_of_times.days*24.0 + diff_of_times.seconds/3600.00) > 24.0 , diff_of_times.days*24.0 + diff_of_times.seconds/3600.00

        settime_counter = 0
        while settime < today:
            settime = obs.next_setting(body,start=startNightTime + datetime.timedelta(days=settime_counter)).datetime()   #if so, then find the very next settime
            settime_counter += 1
            #print(str(settime_counter)+'  \n')
        if settime < risetime: settime = obs.next_setting(body,start=startNightTime + datetime.timedelta(days=1)).datetime()   #if so, then find the very next settime
    
        if longer_than_24hrs_and_Diff(settime - risetime)[1] < 8760.0 and longer_than_24hrs_and_Diff(settime - risetime)[0] == True:
            risetime = obs.next_rising(body,start=startNightTime + datetime.timedelta(days=1)).datetime()   # if difference is less than a year but more than a day, then find the next risetime.  This might be useful because of how many times (settime_counter) I might have changed the settime.  However, it is more than likely that settime_counter will , at most, equal 2. Therefore, I will have only need to bump up risetime once, which is exactly what I am doing.
            transit_t = obs.next_transit(body, start=startNightTime + datetime.timedelta(days=1)).datetime()
        
    
        # at this point settime is after RIGHT NOW-time and after the risetime.  AND sunrise can be after sunset or before sunset
     
    
        if longer_than_24hrs_and_Diff(settime - risetime)[0] == True: observable = False # if the set and rise times are more than one day apart, then this star won't be visible until a year later

        if risetime < sunset and risetime > sunrise and settime < sunset and settime > sunrise: observable = False

        if risetime < subsequent_sunset and risetime > sunrise and settime < subsequent_sunset and settime > sunrise: observable = False
        
        if risetime > not_vis_former and risetime < not_vis_latter and settime > not_vis_former and settime < not_vis_latter: observable = False

    
        # Determine if star is visible between sunset and sunrise.  Also, find the length of time it is visible.
        if observable != False:
            # at this point settime is after RIGHT NOW-time, after the risetime, and after sunset time.   
            observable = True
            if longer_than_24hrs_and_Diff(risetime - sunset)[1] > 0:     # if this Diff is positive, the risetime is after sunset (but BEFORE subsequent_sunset)
                find_cases = np.array([np.abs(longer_than_24hrs_and_Diff(settime - risetime)[1]), np.abs(longer_than_24hrs_and_Diff(subsequent_sunrise - risetime)[1]), np.abs(longer_than_24hrs_and_Diff(settime - subsequent_sunset)[1]), np.abs(longer_than_24hrs_and_Diff(sunrise - risetime)[1]) ])
            
                case = np.where(find_cases == np.min(find_cases))[0][0]

                time_l_visible = find_cases[case]
                if case == 0:
                    start_time = risetime
                    end_time = settime
                    
                elif case == 1:                
                    start_time = risetime
                    end_time = subsequent_sunrise
                
                elif case == 2:
                    start_time = subsequent_sunset
                    end_time = settime
                
                elif case == 3:
                    start_time = risetime
                    end_time = sunrise
                    
            else:                                                        # the risetime is before sunset
                find_least = np.array([np.abs(longer_than_24hrs_and_Diff(settime - sunset)[1]), np.abs(longer_than_24hrs_and_Diff(sunrise - sunset)[1]) ])

                the_least = np.where(find_least == np.min(find_least))[0][0]
                
                time_l_visible = find_least[the_least]
                if the_least == 0:
                    start_time = sunset
                    end_time = settime
                elif the_least == 1:
                    start_time = sunset
                    end_time = sunrise
                    
    
        #trans_meridian_time = datetime.timedelta(0, longer_than_24hrs_and_Diff(settime - risetime)[1] /2.0 *3600.00) + risetime
    
        #following_sunrise = obs.next_rising(ephem.Sun(), start=startNightTime + datetime.timedelta(days=2), use_center=True).datetime()

        if risetime > sunrise and sunrise > sunset:
            midnight_t = datetime.timedelta(0, np.abs(longer_than_24hrs_and_Diff(subsequent_sunrise - subsequent_sunset)[1]) /2.0 *3600.00) + subsequent_sunset
        elif sunrise > sunset and risetime < sunrise:
            midnight_t = datetime.timedelta(0, np.abs(longer_than_24hrs_and_Diff(sunrise - sunset)[1]) /2.0 *3600.00) + sunset
        
        diff_transit_meridian_and_midnight = longer_than_24hrs_and_Diff(transit_t - midnight_t)[1]

    
        #if p_star_indice==116: print('Sunrise : '+str(sunrise)+'\n\nSunset : '+str(sunset)+'\n\nMidnight : '+str(midnight_time)+'\n\n')
        return observable, p_star_indice, time_l_visible, risetime, settime, diff_transit_meridian_and_midnight, transit_t, start_time, end_time

    #all_observables = [choose_primary(p_star_table, p_star_indice=iterate) for iterate in xrange(len(p_star_table)) if choose_primary(p_star_table, p_star_indice=iterate)[0]  == True]

    #best_observables = [all_observables[iterate] for iterate in xrange(len(all_observables)) if all_observables[iterate][2] > 5.0 and np.abs(all_observables[iterate][5]) < 3.0]

    #best_observables = [all_observables[iterate] for iterate in xrange(len(all_observables)) if all_observables[iterate][2] > 4.0 ]
    
    best_observables =  [choose_primary(p_star_table, p_star_indice=iterate) for iterate in [9]]

    #sorted_best_observables = l

    #'''


    def query_stars(p_star_table, best_observables):
        # customize the output of the query results
        customS = Simbad()
        customS.remove_votable_fields('coordinates')
        customS.add_votable_fields('coo(d)')
        customS.add_votable_fields('flux(B)')
        customS.add_votable_fields('flux(V)')
        customS.add_votable_fields('otype(V)')
        customS.ROW_LIMIT = 3  # search for only 4 results/comparison stars
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
            radius = 5
            Comp_table = customS.query_criteria('region(circle, '+p_star_table['ID'][primary_S_indice]+', '+str(radius)+'d)', 'maintypes=*', 'maintypes!=V*', 'Vmag < '+str(mags_V_p[primary_S_indice]), 'Vmag > '+str(mags_V_m[primary_S_indice]))
            while len(Comp_table) < 3:
                radius += 5
                Comp_table = customS.query_criteria('region(circle, '+p_star_table['ID'][primary_S_indice]+', '+str(radius)+'d)', 'maintypes=*', 'maintypes!=V*', 'Vmag < '+str(mags_V_p[primary_S_indice]), 'Vmag > '+str(mags_V_m[primary_S_indice])) # if so, then just search within a 10 degree radius rather than a 5 degree radius
                
            return Comp_table


        comparison_S_tables = [find_comparison(best_observables[p_S_ind][1]) for p_S_ind in xrange(len(best_observables))]
        return comparison_S_tables

    if day_counter == 0: comparison_S_tables = query_stars(p_star_table, best_observables)




    # input the primary and comparison stars into file formatted at a readable Target List

    for telnum in np.array([1, 2, 3, 4]):
    
        # specify the calibrations
        filename = datetime.datetime.strftime(tonight,'n%Y%m%d') + '.T' + str(telnum) + '.txt'
        with open(filename, "w") as outfile:
            outfile.write('{"nbias": 11, "ndark": 11, "nflat": 9, "darkexptime": [300], "flatFilters": ["ip","R"]\, "WaitForMorning": true}\n')
            outfile.write('{"nbiasEnd": 11, "ndarkEnd": 11, "nflatEnd": 9}\n')
        
        for ind in xrange(len(comparison_S_tables)):


            starttime = best_observables[ind][-2]
            #endtime = starttime + datetime.timedelta(0, 1800.0, 000000)  # I just added by 30 minutes
            endtime = best_observables[ind][-1]

        
            if telnum == 1:
                name = p_star_table['ID'][best_observables[ind][1]]
                rahrs = p_star_table['RA (deg)'][best_observables[ind][1]]
                decdegs = p_star_table['Dec (deg)'][best_observables[ind][1]]
            else:
                name = comparison_S_tables[ind]['MAIN_ID'][telnum - 2]
                rahrs = comparison_S_tables[ind]['RA_d'][telnum - 2]
                decdegs = comparison_S_tables[ind]['DEC_d'][telnum - 2]
            
            
            # Add a target to the file:
            target = collections.OrderedDict()
            target['name'] = name
            target['ra'] = rahrs
            target['dec'] = decdegs
            target['starttime'] = datetime.datetime.strftime(starttime,'%Y-%m-%d %H:%M:%S')
            target['endtime'] = datetime.datetime.strftime(endtime,'%Y-%m-%d %H:%M:%S')
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

#'''
        
'''    
#Scraps



        year_s_e = [best_observables[ind][3].year, best_observables[ind][4].year]
        year_s_e = [2015, 2015]
        month_s_e = [06,06]
        day_s_e = [26,27]
        hour_s_e = [18,06]
        minute_s_e = [15,15]
        seconds_s_e = [00,00]
                    
        starttime = datetime.datetime(year_s_e[0],month_s_e[0],day_s_e[0],hour_s_e[0],minute_s_e[0])  # this is just risetime
        endtime = datetime.datetime(year_s_e[1],month_s_e[1],day_s_e[1],hour_s_e[1],minute_s_e[1])   # this is just endtime
        




    decision = raw_input('Do you want to search for a comparison star for all '+str(len(p_star_table['ID']))+' primary target stars?  (Type 1)\n\nOr do you want to only search around 1 primary target star for a comparison star (Type 2)\n---> ')

                   
    if decision == '1':
        comparison_S_tables =  [customS.query_criteria('region(circle, '+p_star_table['ID'][star]+', 5d)', 'maintypes=*', 'maintypes!=V*', 'Vmag < '+str(mags_V_p[star]), 'Vmag > '+str(mags_V_m[star])) for star in xrange(len(p_star_table['ID']))]
    elif decision == '2':
        primary_S_indice  = input('What is the indice of the primary star within the '+sys.argv[1]+' file?\n---> ')
        comparison_table = find_comparison(primary_S_indice)
        print(comparison_table)
    else:
        print('You have typed an invalid number.  Crashing program...')
        return









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
