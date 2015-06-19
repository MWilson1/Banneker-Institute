#!/usr/bin/env python

'''  
Purpose
Query for comparison stars

Input resultant stars into final "target list" file
'''

from astroquery.simbad import Simbad
import numpy as np

# customize the output of the query results
customS = Simbad()
customS.remove_votable_fields('coordinates')
customS.add_votable_fields('flux(B)')
customS.add_votable_fields('flux(V)')
customS.ROW_LIMIT = 7

primaries = ["rho Cnc", "sirius"] # save the list of primary target stars to variable
mag_limit = 1.5
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
comparison_S_tables =  [customS.query_criteria('region(circle, '+primaries[star]+', 5d)', 'maintypes!=V*', 'Bmag < '+str(mags_B_p[star]), 'Bmag > '+str(mags_B_m[star])) for star in xrange(len(primaries))]

print(comparison_S_tables)
