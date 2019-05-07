#!/usr/bin/env python

'''
By David French (frenchd@astro.wisc.edu)

$Id:  plotFancyCosIncDifHist_full.py, v 4.0 05/18/2015

This is the plotFancyCosIncDifHist_full bit from histograms3.py. Now is separated, and loads in a pickle
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
    save = True
    
    
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

    # cos(inclination) histograms for redshifted vs blueshifted distributions of absorbers
    plotCosIncDifHist_full = True
    
    if plotCosIncDifHist_full:
    
#         fig = figure(figsize=(12,10))
#         fig = figure()
#         subplots_adjust(hspace=0.200)

#         subplots_adjust(hspace=0.7)
#         ax = fig.add_subplot(311)
#         subplot(311)
        bins = [0,.10,.20,.30,.40,.50,.60,.70,.80,.90]
        blue = []
        red = []
        
        blueLya = []
        blueLyaErr = []

        redLya = []
        redLyaErr = []

        
        for d,i,l,e in zip(difList,fancyCosIncList,lyaWList,lyaErrorList):
            if d >0:
                # blue shifted absorber, but galaxy is REDSHIFTED
                print 'd: ',d
                blue.append(i)
                blueLya.append(l)
                blueLyaErr.append(e)
            else:
                # red shifted absorber, but galaxy is BLUESHIFTED
                red.append(i)
                redLya.append(l)
                redLyaErr.append(e)
        
#         plot1 = hist([red,blue],bins=bins,histtype='bar',color=['Red','blue'],alpha=0.5)
 #        title('Red shifted aborption: Galaxies')
#         xlabel('Inclination (deg)')
#         ylim(0,6)
#         ylabel('Number')
        fig = figure(figsize=(10,2))
        ax = fig.add_subplot(111)
        subplots_adjust(top = 0.85,bottom = 0.25)
        
        hist(red,bins=bins,histtype='bar',color='red',alpha = 0.9)
        print 'average red: ',average(redLya)
        print 'median red: ', median(redLya)
        print 'average(error): ',average(redLyaErr)
        print 'median(error): ',median(redLyaErr)
        print
        title('Red shifted absorption: Galaxies')
        xlabel('Cos(inclination) = b/a')
#         ylim(0,6)
        ylabel('Number')
        
        if save:
            savefig('{0}/hist(cos(fancy_inclination))_dif_red.pdf'.format(saveDirectory),format='pdf')
        else:
            show()


#         fig.add_subplot(212)
#         subplot(212)
        fig = figure(figsize=(10,2))
        ax = fig.add_subplot(111)
        subplots_adjust(top = 0.85,bottom=0.25)

        hist(blue,bins=bins,histtype='bar',color='Blue',alpha = 0.9)
        print 'average blue: ',average(blueLya)
        print 'median blue: ',median(blueLya)
        print 'avereage(error): ',average(blueLyaErr)
        print 'median(error): ',median(blueLyaErr)
        
        title('Blue shifted absorption: Galaxies')
        xlabel('Cos(inclination) = b/a')
#         ylim(0,5)
        ylabel('Number')
        
        if save:
            savefig('{0}/hist(cos(fancy_inclination))_dif_blue.pdf'.format(saveDirectory),format='pdf')
        else:
            show()

#         ax = fig.add_subplot(312)
#         subplot(312)    
        fig = figure(figsize=(10,2))
        ax = fig.add_subplot(111)
        subplots_adjust(top = 0.85,bottom=0.25)

        hist(allCosFancyInclinations,bins=bins,histtype='bar',color = 'green',alpha=0.9)
        title('Full galaxy sample inclinations')
        xlabel('Cos(inclination) = b/a')
        ylabel('Number')
#         tight_layout()

        if save:
            savefig('{0}/hist(cos(fancy_inclination))_fulldataset.pdf'.format(saveDirectory),format='pdf')
        else:
            show()
    
    
###############################################################################
###############################################################################
###############################################################################
###############################################################################


if __name__=="__main__":
    # do the work
    main()
    