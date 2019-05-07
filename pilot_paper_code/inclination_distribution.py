#!/usr/bin/env python

'''
By David French (frenchd@astro.wisc.edu)

$Id:  inclination_distribution.py, v 1.0 02/24/2015

See how targets are distributed on the sky compared to the galaxies (are we missing
a region of inclination space?)

'''

import sys
# sys.path.append("/usr/users/frenchd/local/lib/python2.6/site-packages/")
import os
import csv
# import string
# import warnings
# import urllib
# import numpy
from pylab import *
import math
from scipy import stats

# from astropy.io.votable import parse,tree

# from vo.table import parse
# import vo.tree

def isNumber(s):
    try:
        float(s)
        if math.isnan(float(s)):
            return False
        else:
            return True
            
    except Exception:
        return False
        
def bfind(k,s):
    # better find (returns true/false)
    if str(k).find(s) == -1:
        return False
    else:
        return True

    
def parseGalaxyNames(nameList):
    # format galaxy names for the url
    
    newNameList = []
    for name in nameList:
        nname = name.strip()
        nname = urllib.quote_plus(nname)
        nname = nname.replace('\n','')
        newNameList.append(nname)
    return newNameList
    

def createCSVTable(outFile,fieldnames):
    # creates and returns a DictReader object populated with header names
    
    writer = csv.DictWriter(outFile, fieldnames=fieldnames)
    headers = dict((n,n) for n in fieldnames)
    writer.writerow(headers)
    return writer


def calculateAngularSeparation(ra1,dec1,ra2,dec2,dist):
    angSep = acos(cos(dec1)*cos(dec2) + sin(dec1)*sin(dec2)*cos(ra1-ra2))
    distSep = float(angSep)*float(dist)*1000
    return distSep
    
    
def fancyInclination(maj,min):
    # calculate more advanced version of inclination
    
    q0 = 0.2
    if float(min) > float(maj):
        maj,min = min, maj
        
    print min, maj
    q = float(min)/float(maj)
    print (q**2 - q0**2)/(1-q0**2)
    inc = math.acos(math.sqrt((q**2 - q0**2)/(1-q0**2)))
    return inc


###########################################################################

    
def main():
    # assuming 'theFile' contains one name per line, read the file
    filename = '/Users/David/Research_Documents/correlationUpdate2.csv'
    theFile = open(filename,'rU')
    reader = csv.DictReader(theFile)

    # instantiate the arrays of data to be collected
#     incArray = []
#     impactArray = []

    incArray = array([])
    aIncArray = array([])
    cosIncArray = array([])
    impactArray = array([])
    
    # loop through the data table and append data to arrays
    for line in reader:
        AGNname = line['AGNname']
        galaxyName = line['galaxyName']
        impact = line['impactParameter (kpc)']
        vcorr = line['vcorrGalaxy (km/s)']
        if isNumber(vcorr):
            vcorr = float(vcorr)
        maj = line['majorAxis (kpc)']
        min = line['minorAxis (kpc)']
        inc = line['inclination (deg)']
        
        
        # ignore 'x' values in the table
        if isNumber(inc) and isNumber(impact):
            if AGNname != galaxyName and float(impact) >1.0:
#                 incArray.append(float(inc))
#                 impactArray.append(float(impact))
                
                if isNumber(maj) and isNumber(min):
                    aInc = fancyInclination(maj,min)
                    cosInc= math.cos(float(aInc)*math.pi/180.)

                    incArray = append(incArray,float(inc))
                    aIncArray = append(aIncArray,float(aInc))
                    impactArray = append(impactArray,float(impact))
                    cosIncArray = append(cosIncArray,float(cosInc))

    print 'mode(inclination): ',stats.mode(incArray)
    print
    print
    print 'mode(impact): ',stats.mode(impactArray)
    
    
    counts, bins = histogram(cosIncArray,bins=9)
    step = abs(bins[0]-bins[1])
    
########################################################################################
    
    # plot the results
    fig = figure()
    ax = gca()
#     ax = fig.add_subplot(111)

#     ax.hist(cosIncArray,bins = 9,histtype='step')
#     show()

    plot1 = ax.scatter(impactArray,aIncArray)
    xlabel('Impact Parameter (kpc)')
    ylabel('New Inclination (deg)')
    xlim(0,500)
    ylim(-1, 91)

    show()
            
            
    # close galaxy table   
    theFile.close()
    
###############################################################################

if __name__=="__main__":
    # do the work
    main()