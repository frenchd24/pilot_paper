#!/usr/bin/env python

'''
By David French (frenchd@astro.wisc.edu)

$Id:  plotImpactHist_Diam.py, v 4.0 05/13/2015

This is the plotImpactHist_Diam bit from histograms3.py. Now is separated, and loads in a pickle
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
    

#########################################################################################
#########################################################################################
    
    # make a histogram of the distribution of impact parameters for associated galaxies 
    # normalized by major diameter
    plotImpactHist_Diam = True
    
    if plotImpactHist_Diam:
        fig = figure(figsize=(10,2))
#         subplots_adjust(hspace=0.200)
        ax = fig.add_subplot(111)
        
#         bins = [0,.10,.20,.30,.40,.50,.60,.70,.80,.90]
    #     bins = [5,15,25,35,45,55,65,75,85]
    #     bins = [0,15,30,45,60,75,90]
        bins = [0,10,20,30,40,50,60,70,80,90]
        
        impactArray = np.array([])
        majArray = np.array([])
        
        # make sure all the values are okay
        for i,m in zip(impactList,majList):
            if isNumber(i) and isNumber(m):
                if i !=-99 and m !=-99:
                    impactArray = append(impactArray,float(i))
                    majArray = append(majArray,float(m))
        
        normalizedImpactArray = impactArray/majArray
        
        plot1 = hist(normalizedImpactArray,bins=bins,histtype='bar')
        
        title('Distribution of impact parameters')
        xlabel('Impact Parameter (kpc)')
        ylabel('Number')
        ax.tick_params(axis='x', labelsize=0)
        ax.tick_params(axis='y',labelsize=8)
#         ylim(0,5)

        if save:
            savefig('{0}/hist(Impact_Diameter).pdf'.format(saveDirectory),format='pdf')
        else:
            show()


###############################################################################
###############################################################################
###############################################################################
###############################################################################


if __name__=="__main__":
    # do the work
    main()
    