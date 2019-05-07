#!/usr/bin/env python

'''
By David French (frenchd@astro.wisc.edu)

$Id:  make_table2.py, v 1.5 10/07/16

Make table 2 for the pilot paper - list of targets with associated galaxies and line info

v1.1: updates for LG_correlation_combined5_8_edit2.csv, where l_min = 0.001
        - (02/17/2016)

v1.2: updates for LG_correlation_combined5_9_edit2.csv, using v_hel instead of vcorr
        - (05/17/2016)
        
v1.3: updated for LG_correlation_combined5_11_25cut_edit.csv (7/21/16)
        
v1.4: updated again for LG_correlation_combined5_11_25cut_edit4.csv (9/27/16)

v1.5: reorder the columns according to the notes in pilot_paper9_bart.pdf (10/07/16)
        

'''

import sys
import os
import csv

from pylab import *
# import atpy
from math import *
from utilities import *
import getpass
import pickle

# from astropy.io.votable import parse,tree

# from vo.table import parse
# import vo.tree


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter

from matplotlib import rc
# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
# ## for Palatino and other serif fonts use:
# #rc('font',**{'family':'serif','serif':['Palatino']})
# rc('text', usetex=True)


def round_to_n(x,n):
    # round to the nth sig fig
    print 'trying: ',x
    return round(x, -int(floor(log10(x))) + (n - 1))
    

###########################################################################

    
def main():
    # assuming 'theFile' contains one name per line, read the file
    
    
    if getpass.getuser() == 'David':
        resultsFilename = '/Users/David/Research_Documents/inclination/git_inclination/LG_correlation_combined5_11_25cut_edit4.csv'
        saveDirectory = '/Users/David/Research_Documents/inclination/git_inclination/pilot_paper_code/'

    elif getpass.getuser() == 'frenchd':
        resultsFilename = '/usr/users/frenchd/inclination/git_inclination/LG_correlation_combined5_11_25cut_edit4.csv'
        saveDirectory = '/usr/users/frenchd/inclination/git_inclination/pilot_paper_code/'

    else:
        print 'Could not determine username. Exiting.'
        sys.exit()

    
    results = open(resultsFilename,'rU')
    reader = csv.DictReader(results)
    
    outFilename = 'table2__5_11_25cut_edit4_2.txt'
    outFile = open(saveDirectory + outFilename,'wt')
    
    virInclude = False
    cusInclude = False
    finalInclude = True
    
    # if match, then the includes in the file have to MATCH the includes above. e.g., if 
    # virInclude = False, cusInclude = True, finalInclude = False, then only systems
    # matching those three would be included. Otherwise, all cusInclude = True would be included
    # regardless of the others
    match = False
    
    # loop through the table now
    for l in reader:
    
        # decide which rows to include
        include_vir = eval(l['include_vir'])
        include_cus = eval(l['include_custom'])
        include = eval(l['include'])
        
        go = False
        if match:
            if virInclude == include_vir and cusInclude == include_cus and finalInclude == include:
                go = True
            else:
                go = False
                
        else:
            if virInclude and include_vir:
                go = True
                
            elif cusInclude and include_cus:
                go = True
                
            elif finalInclude and include:
                go = True
            
            else:
                go = False
        
        
        if go:
            agnName = str(l['AGNname']).replace('_','\_')
            galaxyName = str(l['galaxyName']).replace('_','\_')
            env = str(l['environment'])
            AGNra_dec = str(l['degreesJ2000RA_DecAGN'])
            galaxyRA_Dec = str(l['degreesJ2000RA_DecGalaxy'])
            likelihood = str(l['likelihood'])
            likelihoodm15 = str(l['likelihood_1.5'])
            virialRadius = str(l['virialRadius'])
            m15 = str(l['d^1.5'])
            impact = str(l['impactParameter (kpc)'])
            vcorr = str(l['vcorrGalaxy (km/s)'])
            v_hel = str(l['radialVelocity (km/s)'])
            vel_diff = str(l['vel_diff'])
            galaxyDist = str(l['distGalaxy (Mpc)'])
            maj = str(l['majorAxis (kpc)'])
            min = str(l['minorAxis (kpc)'])
            inc = str(l['inclination (deg)'])
            pa = str(l['positionAngle (deg)'])
            az = str(l['azimuth (deg)'])
            RC3pa = str(l['RC3pa (deg)'])
            morph = str(l['morphology']).replace('_','\_')
            agnRedshift = str(l['AGNredshift'])

            lyaV = str(l['Lya_v'])
            lyaW = str(l['Lya_W']).partition('pm')[0]
            lyaW_err = str(l['Lya_W']).partition('pm')[2]
            b = str(l['b']).partition('pm')[0]
            b_err = str(l['b']).partition('pm')[2]
            na = str(l['Na']).partition(' pm ')[0]
            na_err = str(l['Na']).partition(' pm ')[2]
            lyaWplusErr = str(lyaW)+'$\pm$'+str(lyaW_err)
            
            # the spacer between values
            s = '  &  '
            
            if float(likelihood) >0:
                likelihood = str(round_to_n(float(likelihood),2))
                
            if float(likelihoodm15) >0:
                likelihoodm15 = str(round_to_n(float(likelihoodm15),2))
                
            if float(virialRadius) >0:
                virialRadius = str(round(float(virialRadius),0)).replace('.0','')
                
            impact = str(round(float(impact),0)).replace('.0','')
            vcorr = str(round(float(vcorr),0)).replace('.0','')
            vel_diff = vel_diff
            
            if isNumber(inc):
                inc = str(round(float(inc),0)).replace('.0','')
            
            if isNumber(az):
                az = str(round(float(az),0)).replace('.0','')
                
            if float(likelihood) >= float(likelihoodm15):
                final_likelihood = likelihood
            else:
                final_likelihood = likelihoodm15 + '*'
                
            # the line to be written to file
            line = agnName + s + \
            galaxyName + s + \
            virialRadius + s + \
            v_hel + s + \
            inc + s + \
            az + s + \
            impact + s + \
            lyaV + s + \
            lyaWplusErr + s + \
            vel_diff + s + \
            final_likelihood + '  \\\\'
            
            outFile.write(line + '\n')
            
    outFile.close()
    results.close()

##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################


if __name__=="__main__":
    # do the work
    main()
    