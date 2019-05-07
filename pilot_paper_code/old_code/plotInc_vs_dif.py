#!/usr/bin/env python

'''
By David French (frenchd@astro.wisc.edu)

$Id:  plotInc_vs_dif.py, v 4.0 05/13/2015

This is the plotInc_vs_dif bit from histograms3.py. Now is separated, and loads in a pickle
file of the relevant data, as created by "buildDataLists.py"


Previous (from histograms3.py):
    Plot some stuff for the 100largest initial results

    Make plots for AAS winter 2014 poster. Uses LG_correlation_combined2.csv file

    Updated for the pilot paper (05/06/15)


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


    
def returnLinDiameters(major,minor,distance):
    # input major and minor in arcsec, distance in Mpc
    # outputs major and minor in kpc
    newMajor = math.tan(math.radians(float(major)))*(distance*1000)
    newMinor = math.tan(math.radians(float(minor)))*(distance*1000)
    return (newMajor,newMinor)
    
    

def returnAngDiameters(major,minor,distance):
    # input distances in mpc, major and minor is in kpc
    # outputs angular diameters in arcsec
    newMajor = math.atan((float(major)/1000)/float(distance))*(1/3600)
    newMinor = math.atan((float(minor)/1000)/float(distance))*(1/3600)
    return (newMajor,newMinor)
    
    

###########################################################################

    
def main():
    # assuming 'theFile' contains one name per line, read the file
    
    
    if getpass.getuser() == 'David':
        pickleFilename = '/Users/David/Research_Documents/inclination/pilotData.p'
        saveDirectory = '/Users/David/Research_Documents/inclination/pilot_paper/figures'

    elif getpass.getuser() == 'frenchd':
        pickleFilename = '/usr/users/inclination/pilotData.p'
        saveDirectory = '/usr/users/inclination/pilot_paper/figures'

    else:
        print 'Could not determine username. Exiting.'
        sys.exit()
    
    pickleFile = open(pickleFilename,'rU')
    fullDict = pickle.load(pickleFile)
    
    pickleFile.close()
    
    
    # save each plot?
    save = False
    
    
    # overall structure: fullDict is a dictionary with all the lines and their data in it
    # separated into 'associated' and 'ambiguous' as the two keys. Associated contains
    # all the lists of data for lines associated with a galaxy. Ambiguous contains all
    # the lists of data for lines not unambiguously associated (could be many galaxies
    # or none)
    
    
    # all the lists to be used for associated lines
    lyaVList = fullDict['lyaVList']
    lyaWList = fullDict['lyaWList']
    lyaErrorList = fullDict['lyaErrorList']
    naList = fullDict['naList']
    bList = fullDict['bList']
    impactList = fullDict['impactList']
    azList = fullDict['azList']
    newAzList = fullDict['newAzList']
    incList = fullDict['incList']
    fancyIncList = fullDict['fancyIncList']
    cosIncList = fullDict['cosIncList']
    fancyCosIncList = fullDict['fancyCosIncList']
    paList = fullDict['paList']
    vcorrList = fullDict['vcorrList']
    majList = fullDict['majList']
    difList = fullDict['difList']
    envList = fullDict['envList']
    morphList = fullDict['morphList']
    galaxyNameList = fullDict['galaxyNameList']
    
    # all the lists to be used for ambiguous lines
    lyaVListAmb = fullDict['lyaVListAmb']
    lyaWListAmb = fullDict['lyaWListAmb']
    lyaErrorListAmb = fullDict['lyaErrorListAmb']
    naListAmb = fullDict['naListAmb']
    bListAmb = fullDict['bListAmb']
    impactListAmb = fullDict['impactListAmb']
    azListAmb = fullDict['azListAmb']
    newAzListAmb = fullDict['newAzListAmb']
    incListAmb = fullDict['incListAmb']
    fancyIncListAmb = fullDict['fancyIncListAmb']
    cosIncListAmb = fullDict['cosIncListAmb']
    fancyCosIncListAmb = fullDict['fancyCosIncListAmb']
    paListAmb = fullDict['paListAmb']
    vcorrListAmb = fullDict['vcorrListAmb']
    majListAmb = fullDict['majListAmb']
    difListAmb = fullDict['difListAmb']
    envListAmb = fullDict['envListAmb']
    morphListAmb = fullDict['morphListAmb']
    galaxyNameListAmb = fullDict['galaxyNameListAmb']
    
    # lists for the full galaxy dataset
    allPA = fullDict['allPA']
    allInclinations = fullDict['allInclinations']
    allCosInclinations = fullDict['allCosInclinations']
    allFancyInclinations = fullDict['allFancyInclinations']
    allCosFancyInclinations = fullDict['allCosFancyInclinations']
    
    total = 0
    totalNo = 0
    totalYes = 0
    totalIsolated = 0
    totalGroup = 0
    

########################################################################################
########################################################################################

    # I have no idea what this is or was supposed to do.
    plotInc_vs_dif = True
    
    if plotInc_vs_dif:
        fig = figure()
#         subplots_adjust(hspace=0.200)
        ax = fig.add_subplot(211)
#         bins = [0,10,20,30,40,50,60,70,80,90]
#         bins = [5,15,25,35,45,55,65,75,85,95]

    #     bins = [5,15,25,35,45,55,65,75,85]
#         bins = [0,15,30,45,60,75,90]
#         bins = [0,30,60,90]
        bins = [0,30,60,90]
        incBlue = []
        incRed = []
        for i,d,l in zip(incList,difList,lyaWList):
            if d>0:
                incBlue.append(i)
            else:
                incRed.append(i)
                
        n, bins, patches = hist(incBlue, bins)
        setp(patches, 'facecolor', 'blue', 'alpha', 0.5)               
#         plot1 = hist(incBlue,bins=bins,histtype='bar',c='blue',alpha=0.5)
        title('Blue Shifted Absorbers')
        xlabel('Inclination (deg) / W (AA)')
        ylabel('Number')
#         xlim(0,90)
        
        ax = fig.add_subplot(212)
#         bins = [5,15,25,35,45,55,65,75,85,95]
    #     bins = [5,15,25,35,45,55,65,75,85]
#         bins = [0,15,30,45,60,75,90]
    
        n, bins, patches = hist(incRed, bins)
        setp(patches, 'facecolor', 'red', 'alpha', 0.5)
#         plot2 = hist(incRed,bins=bins,histtype='bar',c='red',alpha=0.5)
        title('Red Shifted Absorbers')
        xlabel('Inclination (deg) / W (AA)')
        ylabel('Number')
#         ylim(0,1)
#         xlim(0,90)
        tight_layout()


        # give me the stats:
        incList2 = []
        azList2 = []
        for i,a in zip(incList,azList):
            if i>=50:
                incList2.append(i)
            if a>50:
                azList2.append(a)
        print '{0} of {1} are inclined >=50%'.format(len(incList2),len(incList))
        print '{0} of {1} have azimuth >=50%'.format(len(azList2),len(azList))
        print 'a: ',a
        print 'len: ',len(incList)
        
        if save:
            savefig('{0}/inc_vs_dif.pdf'.format(saveDirectory),format='pdf')
        else:
            show()
    
    
###############################################################################
###############################################################################
###############################################################################
###############################################################################


if __name__=="__main__":
    # do the work
    main()
    