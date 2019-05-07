#!/usr/bin/env python
"""
By David French (frenchd@astro.wisc.edu)

$Id: createTargetTable_tex.py, v 1.0 07/20/2016

Create table 1 for the pilot paper -> list of targets with ra and dec, z, program ID, 
grating, obs ID, obs date, texp and s/n.

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
        filename = '/Users/David/Research_Documents/correlation/TARGETLIST_6_23_16_reduced.csv'
        resultsName = '/Users/David/Research_Documents/inclination/git_inclination/LG_correlation_combined5_11_25cut_edit.csv'
        outName = '/Users/David/Research_Documents/inclination/git_inclination/table1.txt'

    elif getpass.getuser() == 'frenchd':
        filename = '/usr/users/frenchd/correlation/TARGETLIST_6_23_16_reduced.csv'
        resultsName = '/usr/users/frenchd/inclination/git_inclination/LG_correlation_combined5_11_25cut_edit.csv'
        outName = '/usr/users/frenchd/inclination/git_inclination/table1.txt'
        
    else:
        print 'Could not determine username. Exiting.'
        sys.exit()
    
    
    file = open(filename,'rU')
    reader = csv.DictReader(file)
    
    results = open(resultsName,'rU')
    resultsReader = csv.DictReader(results)
    
    output = open(outName,'wt')
    
    # run through the results file and collect the names of included targets
    targets = {}
    for l in resultsReader:
        include = eval(l['include'])
        
        if include:
            targetName = l['AGNname']
            if not targets.has_key(targetName):
                targets[targetName]=1
                
    targetList = targets.keys()
    
    for t in targetList:
        print t
    
    print len(targetList)
    
    summaryList = []
    summaryList.append(('Target',\
    'R.A.',\
    'Dec.',\
    'z',\
    'Program',\
    'Grating',\
    'Obs ID',\
    'Obs Date',\
    'T_exp [ks]',\
    'S/N [1238]'))
    
    nameList = ['Target']
    raList = ['R.A.']
    decList = ['Dec.']
    zList = ['z']
    programList = ['Program']
    gratingList = ['Grating']
    obsIDList = ['Obs ID']
    obsDateList = ['Obs Date']
    texpList = ['T_exp [ks]']
    snList = ['S/N [1238]']
    
    for l in reader:
        target = l['targetName']
        
        for t in targetList:
            
            if t == target:
                ra = l['ra']
                dec = l['dec']
                z = l['z']
#                 sn = l['SNratio']
                programID = l['programID']
#                 programPI = l['programPI']
#                 instrument = l['instrument']
                grating = l['grating']
                texp = l['exposureTime']
                sn = l['SNexpected']
                
                nameList.append(target)
                raList.append(ra)
                decList.append(dec)
                zList.append(z)
                programList.append(programID)
                gratingList.append(grating)
                obsIDList.append('obsID')
                obsDateList.append('Obs Date')
                texpList.append(texp)
                snList.append(sn)
                
                summaryList.append((target,\
                ra.replace('(','').replace(')','').replace(',',' '),\
                dec.replace('(','').replace(')','').replace(',',' '),\
                z,\
                programID,\
                grating,\
                'Obs ID',\
                'Obs Date',\
                texp,\
                sn))
                

    padding = 4
    widths =[\
    max(len(str(d)) for d in nameList) + padding,\
    max(len(str(d)) for d in raList) + padding,\
    max(len(str(d)) for d in decList) + padding,\
    max(len(str(d)) for d in zList) + padding,\
    max(len(str(d)) for d in programList) + padding,\
    max(len(str(d)) for d in gratingList) + padding,\
    max(len(str(d)) for d in obsIDList) + padding,\
    max(len(str(d)) for d in obsDateList) + padding,\
    max(len(str(d)) for d in texpList) + padding,\
    max(len(str(d)) for d in snList) + padding]

    for row in summaryList:
        output.write("".join(str(i + '  &').ljust(width) for i,width in zip(row,widths))+'\\\\\n')


    file.close()
    results.close()
    output.close()

if __name__=="__main__":
    main()