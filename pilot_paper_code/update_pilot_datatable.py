#!/usr/bin/env python
"""
By David French (frenchd@astro.wisc.edu)

$Id: update_pilot_table.py, v 1.0 09/21/2016

Some entries have been hand-updated in the pilot paper data table: 
    LG_correlation_combined5_11_25cut_edit3.csv
    
This program recalculates things like likelihood, inclination, azimuth, based on updated
diameter & pa values and remakes the table as:
    LG_correlation_combined5_11_25cut_edit4.csv

"""

# from __future__ import division
import optparse
import sys
import os
# import tempfile
import csv
import string
from math import *
import numpy
import getpass
import correlateSingle7 as correlateSingle
from utilities import *

    
################################################################

def main():
    # This function reformats Bart's file
    
    if getpass.getuser() == 'David':
        fileName = '/Users/David/Research_Documents/inclination/git_inclination/LG_correlation_combined5_11_25cut_edit3.csv'
        outName = '/Users/David/Research_Documents/inclination/git_inclination/LG_correlation_combined5_11_25cut_edit4.csv'

    elif getpass.getuser() == 'frenchd':
        fileName = '/usr/users/frenchd/inclination/git_inclination/LG_correlation_combined5_11_25cut_edit3.csv'
        outName = '/usr/users/frenchd/inclination/git_inclination/LG_correlation_combined5_11_25cut_edit4.csv'
        
    else:
        print 'Could not determine username. Exiting.'
        sys.exit()
    
    file = open(fileName,'rU')
    reader = csv.DictReader(file)
    
    fieldnames = ('AGNname',\
    'center',\
    'galaxyName',\
    'environment',\
    'degreesJ2000RA_DecAGN',\
    'degreesJ2000RA_DecGalaxy',\
    'likelihood',\
    'likelihood_1.5',\
    'virialRadius',\
    'd^1.5',\
    'impactParameter (kpc)',\
    'redshiftDistances',\
    'vcorrGalaxy (km/s)',\
    'radialVelocity (km/s)',\
    'vel_diff',\
    'distGalaxy (Mpc)',\
    'AGN S/N',\
    'majorAxis (kpc)',\
    'minorAxis (kpc)',\
    'acos[b/a] (deg)',\
    'inclination (deg)',\
    'positionAngle (deg)',\
    'azimuth (deg)',\
    'RC3flag',\
    'RC3type',\
    'RC3inc (deg)',\
    'RC3pa (deg)',\
    'morphology',\
    'final_morphology',\
    'galaxyRedshift',
    'AGNredshift',\
    'spectrumStatus',\
    'include',\
    'include_vir',\
    'include_custom',\
    'Lya_v',\
    'vlimits',\
    'Lya_W',\
    'Na',\
    'b',\
    'identified',\
    'comment')
    
    writerOutFile = open(outName,'wt')
    writer = csv.DictWriter(writerOutFile, fieldnames=fieldnames)
    headers = dict((n,n) for n in fieldnames)
    writer.writerow(headers)

    for l in reader:
        q0=0.2
    
        major = float(l['majorAxis (kpc)'])
        minor = float(l['minorAxis (kpc)'])
        AGNra,AGNdec = eval(l['degreesJ2000RA_DecAGN'])
        galRA,galDec = eval(l['degreesJ2000RA_DecGalaxy'])
        dist = float(l['distGalaxy (Mpc)'])
        pa = float(l['positionAngle (deg)'])
        impact = float(l['impactParameter (kpc)'])
        az = l['azimuth (deg)']
        vel_dif = float(l['vel_diff'])
        
        aCos = calculateInclination(major,minor)
        inc = calculateFancyInclination(major,minor,q0)
        az = calculateAzimuth(galRA,galDec,AGNra,AGNdec,dist,pa)
        
        rVir = calculateVirialRadius(major)

        # try this "sphere of influence" value instead
        m15 = major**1.5

        # first for the virial radius
        likelihood = math.exp(-(impact/rVir)**2) * math.exp(-(vel_dif/200.)**2)
        
        if rVir>= impact:
            likelihood = likelihood*2
            
        # then for the second 'virial like' m15 radius
        likelihoodm15 = math.exp(-(impact/m15)**2) * math.exp(-(vel_dif/200.)**2)
        
        if m15>= impact:
            likelihoodm15 = likelihoodm15*2
                
            # should be like 33% at 1R_v, and linear down after that.
            # 66% at 0.5R_v, something quadratic, or logarithmic
            # look up the "sphere of influence" for where the probability drops
            # down to 10%
            
            # use M/L, surpress environment by M/L ratios
            # also include velocity difference to make the virial radius 3D - 
            # really a 3D impact parameter
            
        
        objectInfoList = [l['AGNname'],\
        l['center'],\
        l['galaxyName'],\
        l['environment'],\
        l['degreesJ2000RA_DecAGN'],\
        l['degreesJ2000RA_DecGalaxy'],\
        likelihood,\
        likelihoodm15,\
        rVir,\
        m15,\
        impact,\
        l['redshiftDistances'],\
        l['vcorrGalaxy (km/s)'],\
        l['radialVelocity (km/s)'],\
        l['vel_diff'],\
        dist,\
        l['AGN S/N'],\
        major,\
        minor,\
        aCos,\
        inc,\
        pa,\
        az,\
        l['RC3flag'],\
        l['RC3type'],\
        l['RC3inc (deg)'],\
        l['RC3pa (deg)'],\
        l['morphology'],\
        l['morphology'],\
        l['galaxyRedshift'],\
        l['AGNredshift'],\
        l['spectrumStatus'],\
        l['include'],\
        l['include_vir'],\
        l['include_custom'],\
        l['Lya_v'],\
        l['vlimits'],\
        l['Lya_W'],\
        l['Na'],\
        l['b'],\
        l['identified'],\
        l['comment']]

        row = dict((f,o) for f,o in zip(fieldnames,objectInfoList))
        writer.writerow(row)
                
    file.close()
    writerOutFile.close()

if __name__=="__main__":
    main()
