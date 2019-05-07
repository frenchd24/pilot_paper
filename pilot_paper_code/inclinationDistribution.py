 #!/usr/bin/env python

'''
By David French (frenchd@astro.wisc.edu)

$Id:  inclinationDistribution.py, v 1.0 04/01/2015

Plot the distribution of inclination angles in the newest correlation tables:
    NGT3-TG6_500Correlation_full_500cutoff.csv
    NGT3-TG6_2000Correlation_full_500cutoff.csv
    
Compare to the values in the full table, and to a random sampling

'''

import sys
import os
import csv
# import string
# import warnings
# import urllib
# import numpy
from pylab import *
# import atpy
import math
import getpass
import random

# import scipy.optimize as optimization
# import pickle
# import itertools


# import numpy as np
# import matplotlib.pyplot as plt
# from matplotlib.ticker import NullFormatter


def isNumber(s):
    # determine if something is a number or not. returns a boolean

    try:
        float(s)
        if math.isnan(float(s)):
            return False
        else:
            return True
    except Exception,e:
        return False


def bfind(obj, st):
    # "better find". Just returns true or false instead of an index
    
    if str(obj).find(st) !=-1:
        return True
    else:
        return False



def fancyInclination(maj,min):
    # calculate more advanced version of inclination
    
    q0 = 0.2
    if float(min) > float(maj):
        maj,min = min, maj
        
#     print min, maj
    q = float(min)/float(maj)
#     print (q**2 - q0**2)/(1-q0**2)
    inc = math.acos(math.sqrt((q**2 - q0**2)/(1-q0**2)))
    return inc



##########################################################################################
##########################################################################################

def main():
    
    
    # the correlation table filename and path   
    if getpass.getuser() == 'David':
        correlationFilename = '/Users/David/Research_Documents/gt/NGT3-TG6_500Correlation_full_500cutoff.csv'
        galaxyFilename = '/Users/David/Research_Documents/gt/NewGalaxyTable4.csv'
        
    elif getpass.getuser() == 'frenchd':
        correlationFilename = '/usr/users/frenchd/gt/NGT3-TG6_500Correlation_full_500cutoff.csv'
        galaxyFilename = '/usr/users/frenchd/gt/NewGalaxyTable4.csv'
    else:
        print 'Could not determine username. Exiting.'
        sys.exit()
    
    # open the correlation file
    theFile = open(correlationFilename,'rU')
    galaxyFile = open(galaxyFilename,'rU')
    reader = csv.DictReader(theFile)
    galaxyReader = csv.DictReader(galaxyFile)

    # loop through the correlation table
    incArray = array([])
    fIncArray = array([])
    for l in reader:
        inc = l['inclination (deg)']
        maj = l['majorAxis (kpc)']
        min = l['minorAxis (kpc)']
        if isNumber(maj) and isNumber(min):
            try:
                fInc = fancyInclination(maj,min)
                fIncArray = append(fIncArray,float(fInc))
            except Exception, e:
#                 sys.stdout.write('{0}; tried for {1}, {2}.'.format(e,maj,min))
                pass
            
        if isNumber(inc):
            incArray = append(incArray,float(inc))


    count = 0
    tots = 0
    for n in incArray:
        if n >=70:
            count +=1
            
        tots+=1
    
    print 'total: ',tots
    print 'counts: ',count
    print 'fraction with inc > 70: ',float(count)/tots

    # the full table
    fullIncArray = array([])
    fullfIncArray = array([])
    for l in galaxyReader:
        inc = l['inclination (deg)']
        maj,min = eval(l['angDiameters (arcsec)'])
        if isNumber(maj) and isNumber(min):
            try:
                fInc = fancyInclination(maj,min)
                fullfIncArray = append(fullfIncArray,float(fInc))
            except Exception, e:
#                 sys.stdout.write('{0}; tried for {1}, {2}.'.format(e,maj,min))
                pass
        
        if isNumber(inc):
            fullIncArray = append(fullIncArray, float(inc))
            
            
    # create a pseudo-random sample
    randIncArray = array([])
    for i in range(len(fullIncArray)):
        num = random.uniform(0,90.0)
        
        randIncArray = append(randIncArray,num)
            

    # bin up all the counts
#     counts, bins = histogram(currentData,bins=50)
#     step = abs(bins[0]-bins[1])


    fig = figure()
    ax = gca()
    ax = fig.add_subplot(311)
    
    ax.hist(incArray,bins = 45,histtype='step')
    ylabel('Number')
#     xlabel('inclination angle')
    title('inclination distribution for correlation'.format(correlationFilename))
#     legend(loc='lower center')

    ax = fig.add_subplot(312)
    
    ax.hist(fullIncArray,bins = 45,histtype='step')
    ylabel('Number')
#     xlabel('inclination angle')
    title('full inclination distribution'.format(correlationFilename))
#     legend(loc='lower center')

    ax = fig.add_subplot(313)
    
    ax.hist(randIncArray,bins = 45,histtype='step')
    ylabel('Number')
    xlabel('inclination angle')
    title('random inclination distribution')
#     legend(loc='lower center')

    show()
    
    
    # close the correlation table
    theFile.close()
    galaxyFile.close()


if __name__=="__main__":
    main()
    
    