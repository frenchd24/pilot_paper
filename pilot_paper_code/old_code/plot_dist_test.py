#!/usr/bin/env python

'''
By David French (frenchd@astro.wisc.edu)

$Id:  plot_dist_test.py, v 1.0 08/25/2015

Compare the distributions of inclination angles between the datasets: full galaxy table,
blueshifted, and redshifted absorption


'''

import sys
import os
import csv

from pylab import *
from scipy import stats
# import atpy
from math import *
from utilities import *
import getpass
import pickle

# from astropy.io.votable import parse,tree
# import numpy as np
# import matplotlib.pyplot as plt
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
    
##########################################################################################
##########################################################################################
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
    
    
    lyaV_blue = []
    lyaV_red = []
    lyaW_blue = []
    lyaW_red = []
    lyaErr_blue = []
    lyaErr_red = []
    na_blue = []
    na_red = []
    b_blue = []
    b_red = []
    impact_blue = []
    impact_red = []
    az_blue = []
    az_red = []
    newAz_blue = []
    newAz_red = []
    inc_blue = []
    inc_red = []
    fancyInc_blue = []
    fancyInc_red = []
    cosInc_blue = []
    cosInc_red = []
    fancyCosInc_blue = []
    fancyCosInc_red = []
    pa_blue = []
    pa_red = []
    vcorr_blue = []
    vcorr_red = []
    maj_blue = []
    maj_red = []
    dif_blue = []
    dif_red = []
    env_blue = []
    env_red = []
    morph_blue = []
    morph_red = []
    
    c = -1
    for d in difList:
        c +=1
        if d > 0:
            # blueshifted absorption
            lyaV_blue.append(lyaVList[c])
            lyaW_blue.append(lyaWList[c])
            lyaErr_blue.append(lyaErrorList[c])
            na_blue.append(naList[c])
            b_blue.append(bList[c])
            impact_blue.append(impactList[c])
            az_blue.append(azList[c])
            newAz_blue.append(newAzList[c])
            inc_blue.append(incList[c])
            fancyInc_blue.append(fancyIncList[c])
            cosInc_blue.append(cosIncList[c])
            fancyCosInc_blue.append(fancyCosIncList[c])
            pa_blue.append(paList[c])
            vcorr_blue.append(vcorrList[c])
            maj_blue.append(majList[c])
            dif_blue.append(difList[c])
            env_blue.append(envList[c])
            morph_blue.append([c])
            
        else:
            # redshifted absorption
            lyaV_red.append(lyaVList[c])
            lyaW_red.append(lyaWList[c])
            lyaErr_red.append(lyaErrorList[c])
            na_red.append(naList[c])
            b_red.append(bList[c])
            impact_red.append(impactList[c])
            az_red.append(azList[c])
            newAz_red.append(newAzList[c])
            inc_red.append(incList[c])
            fancyInc_red.append(fancyIncList[c])
            cosInc_red.append(cosIncList[c])
            fancyCosInc_red.append(fancyCosIncList[c])
            pa_red.append(paList[c])
            vcorr_red.append(vcorrList[c])
            maj_red.append(majList[c])
            dif_red.append(difList[c])
            env_red.append(envList[c])
            morph_red.append(morphList[c])
    
        
        
##########################################################################################
##########################################################################################
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
    
    lyaV_blueAmb = []
    lyaV_redAmb = []
    lyaW_blueAmb = []
    lyaW_redAmb = []
    lyaErr_blueAmb = []
    lyaErr_redAmb = []
    na_blueAmb = []
    na_redAmb = []
    b_blueAmb = []
    b_redAmb = []
    impact_blueAmb = []
    impact_redAmb = []
    az_blueAmb = []
    az_redAmb = []
    newAz_blueAmb = []
    newAz_redAmb = []
    inc_blueAmb = []
    inc_redAmb = []
    fancyInc_blueAmb = []
    fancyInc_redAmb = []
    cosInc_blueAmb = []
    cosInc_redAmb = []
    fancyCosInc_blueAmb = []
    fancyCosInc_redAmb = []
    pa_blueAmb = []
    pa_redAmb = []
    vcorr_blueAmb = []
    vcorr_redAmb = []
    maj_blueAmb = []
    maj_redAmb = []
    dif_blueAmb = []
    dif_redAmb = []
    env_blueAmb = []
    env_redAmb = []
    morph_blueAmb = []
    morph_redAmb = []
    
    
    c = -1
    for d in difListAmb:
        c +=1
        if d > 0:
            # blueshifted absorption
            lyaV_blueAmb.append(lyaVListAmb[c])
            lyaW_blueAmb.append(lyaWListAmb[c])
            lyaErr_blueAmb.append(lyaErrorListAmb[c])
            na_blueAmb.append(naListAmb[c])
            b_blueAmb.append(bListAmb[c])
            impact_blueAmb.append(impactListAmb[c])
            az_blueAmb.append(azListAmb[c])
            newAz_blueAmb.append(newAzListAmb[c])
            inc_blueAmb.append(incListAmb[c])
            fancyInc_blueAmb.append(fancyIncListAmb[c])
            cosInc_blueAmb.append(cosIncListAmb[c])
            fancyCosInc_blueAmb.append(fancyCosIncListAmb[c])
            pa_blueAmb.append(paListAmb[c])
            vcorr_blueAmb.append(vcorrListAmb[c])
            maj_blueAmb.append(majListAmb[c])
            dif_blueAmb.append(difListAmb[c])
            env_blueAmb.append(envListAmb[c])
            morph_blueAmb.append(morphListAmb[c])
            
        else:
            # redshifted absorption
            lyaV_redAmb.append(lyaVListAmb[c])
            lyaW_redAmb.append(lyaWListAmb[c])
            lyaErr_redAmb.append(lyaErrorListAmb[c])
            na_redAmb.append(naListAmb[c])
            b_redAmb.append(bListAmb[c])
            impact_redAmb.append(impactListAmb[c])
            az_redAmb.append(azListAmb[c])
            newAz_redAmb.append(newAzListAmb[c])
            inc_redAmb.append(incListAmb[c])
            fancyInc_redAmb.append(fancyIncListAmb[c])
            cosInc_redAmb.append(cosIncListAmb[c])
            fancyCosInc_redAmb.append(fancyCosIncListAmb[c])
            pa_redAmb.append(paListAmb[c])
            vcorr_redAmb.append(vcorrListAmb[c])
            maj_redAmb.append(majListAmb[c])
            dif_redAmb.append(difListAmb[c])
            env_redAmb.append(envListAmb[c])
            morph_redAmb.append(morphListAmb[c])
    
    
##########################################################################################
##########################################################################################
    
    # lists for the full galaxy dataset
    
    allPA = fullDict['allPA']
    allInclinations = fullDict['allInclinations']
    allCosInclinations = fullDict['allCosInclinations']
    allFancyInclinations = fullDict['allFancyInclinations']
    allFancyCosInclinations = fullDict['allCosFancyInclinations']
    
    total = 0
    totalNo = 0
    totalYes = 0
    totalIsolated = 0
    totalGroup = 0
    

########################################################################################
########################################################################################

    # plot histograms of the cos(inclinations) for both associated galaxies and the 
    # full galaxy data set, combining both redshifted and blueshifted
    plot_dist_cosInc = False
    
    if plot_dist_cosInc:
    
        '''
        Here's an example:
        
        n1 = 200
        n2 = 300
        
        rvs1 = stats.norm.rvs(size=n1, loc=0., scale=1)
        rvs2 = stats.norm.rvs(size=n2, loc=0.5, scale=1.5)
        
        
        ans = stats.ks_2samp(rvs1, rvs2)
        print 'ans: ',ans
        print
        
        '''
    
        # define the datasets
        rvs1all = cosIncList
        rvs1 = []
        rvs2 = allCosInclinations

        # remove -99 'no-data' values
        for i in rvs1all:
            if float(i) >=0:
                rvs1.append(i)
                
        # do the K-S test and print the results
        ans = stats.ks_2samp(rvs1, rvs2)
        print 'real ans: ',ans 
        
        
        # plot the distributions
        fig = figure()
        
        bins = 15
        
        ax1 = fig.add_subplot(211)
        plot1 = hist(rvs1,bins=bins)
        xlim(0,1)
        
        ax2 = fig.add_subplot(212)
        plot1 = hist(rvs2,bins=bins)
        xlim(0,1)
        
        show()

        

########################################################################################
########################################################################################

    # plot histograms of the inclinations for both associated galaxies and the 
    # full galaxy data set, combining both redshifted and blueshifted
    plot_dist_inc = False
    
    if plot_dist_inc:
    
        # define the datasets
        rvs1all = incList
        rvs1 = []
        rvs2 = allInclinations
        
        # remove -99 'no-data' values
        for i in rvs1all:
            if float(i) >=0:
                rvs1.append(i)
                
                
        # perform the K-S test
        ans = stats.ks_2samp(rvs1, rvs2)
        print 'real ans: ',ans 
        
        
        # plot the distributions 
        fig = figure()
        
        bins = 15
        
        ax1 = fig.add_subplot(211)
        plot1 = hist(rvs1,bins=bins)
        
        ax2 = fig.add_subplot(212)
        plot1 = hist(rvs2,bins=bins)
        
        show()
        


########################################################################################
########################################################################################

    # plot histograms of the inclinations for both associated galaxies and the 
    # full galaxy data set, combining both redshifted and blueshifted
    plot_dist_fancyInc = False
    
    if plot_dist_fancyInc:
    
        # define the datasets
        rvs1all = fancyIncList
        rvs1 = []
        rvs2 = allFancyInclinations
        
        # remove -99 'no-data' values
        for i in rvs1all:
            if float(i) >=0:
                rvs1.append(i)
                
                
        # perform the K-S test
        ans = stats.ks_2samp(rvs1, rvs2)
        print 'real ans: ',ans 
        
        
        # plot the distributions 
        fig = figure()
        
        bins = 15
        
        ax1 = fig.add_subplot(211)
        plot1 = hist(rvs1,bins=bins)
        
        ax2 = fig.add_subplot(212)
        plot1 = hist(rvs2,bins=bins)
        
        show()
        

########################################################################################
########################################################################################

    # plot histograms of the inclinations for both associated galaxies and the 
    # full galaxy data set, combining both redshifted and blueshifted
    plot_dist_fancyCosInc = False
    
    if plot_dist_fancyCosInc:
    
        # define the datasets
        rvs1all = fancyCosIncList
        rvs1 = []
        rvs2 = allFancyCosInclinations
        
        # remove -99 'no-data' values
        for i in rvs1all:
            if float(i) >=0:
                rvs1.append(i)
                
                
        # perform the K-S test
        ans = stats.ks_2samp(rvs1, rvs2)
        print 'real ans: ',ans 
        
        
        # plot the distributions 
        fig = figure()
        
        bins = 15
        
        ax1 = fig.add_subplot(211)
        plot1 = hist(rvs1,bins=bins)
        
        ax2 = fig.add_subplot(212)
        plot1 = hist(rvs2,bins=bins)
        
        show()


########################################################################################
########################################################################################

    # plot histograms of the inclinations for both associated galaxies and the 
    # full galaxy data set, combining both redshifted and blueshifted
    plot_dist_fancyCosInc_red_blue = False
    
    if plot_dist_fancyCosInc_red_blue:
    
        # define the datasets
        rvs1all = fancyCosInc_blue
        rvs1 = []
        rvs2all = fancyCosInc_red
        rvs2 = []
        rvs3 = allFancyCosInclinations
        
        # remove -99 'no-data' values
        for i in rvs1all:
            if float(i) >=0:
                rvs1.append(i)

        for k in rvs2all:
            if float(k) >=0:
                rvs2.append(k)
                
                
        # perform the K-S test
        ans1 = stats.ks_2samp(rvs1, rvs2)
        print 'blue vs red: ',ans1
        
        ans2 = stats.ks_2samp(rvs1, rvs3)
        print 'blue vs all: ',ans2
        
        ans3 = stats.ks_2samp(rvs2, rvs3)
        print 'red vs all: ',ans3
        
        # plot the distributions 
        fig = figure()
        
        bins = 15
        
        ax1 = fig.add_subplot(311)
        plot1 = hist(rvs1,bins=bins)
        title('blueshifted Cos(fancy_inc)')
        
        ax2 = fig.add_subplot(312)
        plot2 = hist(rvs2,bins=bins)
        title('redshifted Cos(fancy_inc)')
        
        ax3 = fig.add_subplot(313)
        plot3 = hist(rvs3,bins=bins)
        title('Full galaxy table Cos(fancy_inc)')
        
        show()



########################################################################################
########################################################################################

    # plot histograms of the inclinations for both associated galaxies and the 
    # full galaxy data set, combining both redshifted and blueshifted
    plot_dist_cosInc_red_blue = True
    
    if plot_dist_cosInc_red_blue:
    
        # define the datasets
        rvs1all = cosInc_blue
        rvs1 = []
        rvs2all = cosInc_red
        rvs2 = []
        rvs3 = allCosInclinations
        
        # remove -99 'no-data' values
        for i in rvs1all:
            if float(i) >=0:
                rvs1.append(i)

        for k in rvs2all:
            if float(k) >=0:
                rvs2.append(k)
                
                
        # perform the K-S test
        ans1 = stats.ks_2samp(rvs1, rvs2)
        ans1a = stats.anderson_ksamp([rvs1,rvs2])
        print 'blue vs red: ',ans1, ans1a
        
        ans2 = stats.ks_2samp(rvs1, rvs3)
        print 'blue vs all: ',ans2
        
        ans3 = stats.ks_2samp(rvs2, rvs3)
        print 'red vs all: ',ans3
        
        # plot the distributions 
        fig = figure(figsize=(8,8))
        subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.4)

        bins = 15
        
        ax1 = fig.add_subplot(311)
        plot1 = hist(rvs1,bins=bins)
        title('blueshifted Cos(inc)')
        
        ax2 = fig.add_subplot(312)
        plot2 = hist(rvs2,bins=bins)
        title('redshifted Cos(inc)')
        
        ax3 = fig.add_subplot(313)
        plot3 = hist(rvs3,bins=bins)
        title('Full galaxy table Cos(inc)')
        
        show()



        
#         fig = figure()
# #         subplots_adjust(hspace=0.200)
#         ax = fig.add_subplot(211)
#         bins = [0,.10,.20,.30,.40,.50,.60,.70,.80,.90]
#         subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.4)
# 
#         plot1 = hist(cosIncList,bins=bins,histtype='bar')
#         title('Absorber-associated galaxies cos(inclination)')
#         xlabel('Inclination (deg)')
#         ylabel('Number')
# 
#         ax = fig.add_subplot(212)
#         bins = [0,.10,.20,.30,.40,.50,.60,.70,.80,.90]
#         plot1 = hist(allCosInclinations,bins=bins,histtype='bar')
#         title('Full galaxy sample cos(inclination)')
#         xlabel('Inclination (deg)')
#         ylabel('Number')
# #         tight_layout()
# 
#         if save:
#             savefig('{0}/inc_dist.pdf'.format(saveDirectory),format='pdf')
#         else:
#             show()
    
    
###############################################################################
###############################################################################
###############################################################################
###############################################################################


if __name__=="__main__":
    # do the work
    main()
    