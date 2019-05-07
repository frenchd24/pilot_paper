 #!/usr/bin/env python

'''
By David French (frenchd@astro.wisc.edu)

$Id:  azimuthDistribution.py, v 1.0 03/18/2015

Plot the distribution of azimuth angles in the newest correlation tables:
    NGT3-TG6_500Correlation_full_500cutoff.csv
    NGT3-TG6_2000Correlation_full_500cutoff.csv

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


##########################################################################################
##########################################################################################

def main():
    
    
    # the correlation table filename and path   
    if getpass.getuser() == 'David':
        correlationFilename = '/Users/David/Research Documents/gt/NGT3-TG6_500Correlation_full_500cutoff.csv'
        
    elif getpass.getuser() == 'frenchd':
        correlationFilename = '/usr/users/frenchd/gt/NGT3-TG6_500Correlation_full_500cutoff.csv'
    else:
        print 'Could not determine username. Exiting.'
        sys.exit()
    
    # open the correlation file
    theFile = open(correlationFilename,'rU')
    reader = csv.DictReader(theFile)

    # loop through the tables
    azArray = array([])
    for l in reader:
        az = l['azimuth (deg)']
        
        if isNumber(az):
            azArray = append(azArray,float(az))
    


    # bin up all the counts
#     counts, bins = histogram(currentData,bins=50)
#     step = abs(bins[0]-bins[1])


    fig = figure()
    ax = gca()

#     plot1 = ax.plot(bins[:-1],counts,lw=0,marker='.',ms=5)

    ax.hist(azArray,bins = 45,histtype='step')
    
#     plot2 = ax.plot(m,s,c='Red')
#     plot3 = ax.plot(m,sFit,c='Green',label='phif={0}, alphaf={1}, mstarf={2}'.format(int(phiFit),alphaFit,round(mstarFit,3)))
#     plot3 = ax.plot(m,sGuess,c='red',label='phi={0}, alpha={1}, mstar={2}'.format(phiGuess,alphaGuess,mstarGuess))
#     ax.invert_xaxis()

    ylabel('Number')
    xlabel('Azimuth angle')
    ax.set_yscale('log')
    title('Azimuth angle distribution for {0}'.format(correlationFilename))
    legend(loc='lower center')
    ax.invert_xaxis()

    show()
    
    
    
    # close the correlation table
    theFile.close()


if __name__=="__main__":
    main()
    
    