#!/usr/bin/env python

'''
By David French (frenchd@astro.wisc.edu)

$Id:  plot_dist_test2.py, v 2.0 12/28/2015

Compare the distributions of inclination angles between the datasets: full galaxy table,
blueshifted, and redshifted absorption

Also for inclination, fancy_inclination, cos(inclination), cos(fancy_inclination)


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


###########################################################################

    
def main():

    
    if getpass.getuser() == 'David':
        pickleFilename = '/Users/David/Research_Documents/inclination/git_inclination/pilot_paper_code/pilotData2.p'
        resultsFilename = '/Users/David/Research_Documents/inclination/git_inclination/LG_correlation_combined5_3.csv'
        saveDirectory = '/Users/David/Research_Documents/inclination/git_inclination/pilot_paper_code/plots/'

    elif getpass.getuser() == 'frenchd':
        pickleFilename = '/usr/users/frenchd/inclination/git_inclination/pilot_paper_code/pilotData2.p'
        resultsFilename = '/usr/users/frenchd/inclination/git_inclination/LG_correlation_combined5_3.csv'
        saveDirectory = '/usr/users/frenchd/inclination/git_inclination/pilot_paper_code/plots/'

    else:
        print 'Could not determine username. Exiting.'
        sys.exit()
    

    pickleFile = open(pickleFilename,'rU')
    fullDict = pickle.load(pickleFile)
    
    pickleFile.close()
    
    # save each plot?
    save = False
    
    results = open(resultsFilename,'rU')
    reader = csv.DictReader(results)
    
    virInclude = False
    cusInclude = False
    finalInclude = True
    
    # if match, then the includes in the file have to MATCH the includes above. e.g., if 
    # virInclude = False, cusInclude = True, finalInclude = False, then only systems
    # matching those three would be included. Otherwise, all cusInclude = True would be included
    # regardless of the others
    match = False
    
    # all the lists to be used for associated lines
    lyaVList = []
    lyaWList = []
    naList = []
    bList = []
    impactList = []
    azList = []
    incList = []
    fancyIncList = []
    cosIncList = []
    cosFancyIncList = []
    paList = []
    vcorrList = []
    majList = []
    difList = []
    envList = []
    morphList = []
    m15List = []
    virList = []
    likeList = []
    likem15List = []
    
    
    for l in reader:
        include_vir = eval(l['include_vir'])
        include_cus = eval(l['include_custom'])
        include = eval(l['include'])
        
        go = False
        if match:
            if virInclude == include_vir and cusInclude == include_cus:
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
            AGNra_dec = eval(l['degreesJ2000RA_DecAGN'])
            galaxyRA_Dec = eval(l['degreesJ2000RA_DecGalaxy'])
            lyaV = l['Lya_v']
            lyaW = l['Lya_W'].partition('pm')[0]
            lyaW_err = l['Lya_W'].partition('pm')[2]
            env = l['environment']
            galaxyName = l['galaxyName']
            impact = l['impactParameter (kpc)']
            galaxyDist = l['distGalaxy (Mpc)']
            pa = l['positionAngle (deg)']
            RC3pa = l['RC3pa (deg)']
            morph = l['morphology']
            vcorr = l['vcorrGalaxy (km/s)']
            maj = l['majorAxis (kpc)']
            min = l['minorAxis (kpc)']
            inc = l['inclination (deg)']
            az = l['azimuth (deg)']
            b = l['b'].partition('pm')[0]
            b_err = l['b'].partition('pm')[2]
            na = eval(l['Na'].partition(' pm ')[0])
            print "l['Na'].partition(' pm ')[2] : ",l['Na'].partition(' pm ')
            na_err = eval(l['Na'].partition(' pm ')[2])
            likelihood = l['likelihood']
            likelihoodm15 = l['likelihood_1.5']
            virialRadius = l['virialRadius']
            m15 = l['d^1.5']
            vel_diff = l['vel_diff']
            
            if isNumber(RC3pa) and not isNumber(pa):
                pa = RC3pa
            
            if isNumber(inc):
                cosInc = cos(float(inc) * pi/180.)
                
                if isNumber(maj) and isNumber(min):
                    q0 = 0.2
                    fancyInc = calculateFancyInclination(maj,min,q0)
                    cosFancyInc = cos(fancyInc * pi/180)
                else:
                    fancyInc = -99
                    cosFancyInc = -99
            else:
                cosInc = -99
                inc = -99
                fancyInc = -99
                cosFancyInc = -99
            
            # all the lists to be used for associated lines
            lyaVList.append(float(lyaV))
            lyaWList.append(float(lyaW))
            naList.append(na)
            bList.append(float(b))
            impactList.append(float(impact))
            azList.append(az)
            incList.append(float(inc))
            fancyIncList.append(fancyInc)
            cosIncList.append(cosInc)
            cosFancyIncList.append(cosFancyInc)
            paList.append(pa)
            vcorrList.append(vcorr)
            majList.append(maj)
            difList.append(float(vel_diff))
            envList.append(float(env))
            morphList.append(morph)
            m15List.append(m15)
            virList.append(virialRadius)
            likeList.append(likelihood)
            likem15List.append(likelihoodm15)

    results.close()
        
    
##########################################################################################
##########################################################################################
    
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
        print 'KS for cosIncList vs all: ',ans 
        
        
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
        print 'KS for incList vs all: ',ans 
        
        
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
        print 'KS for fancyIncList vs all: ',ans 
        
        
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
        rvs1all = cosFancyIncList
        rvs1 = []
        rvs2 = allCosFancyInclinations
        
        # remove -99 'no-data' values
        for i in rvs1all:
            if float(i) >=0:
                rvs1.append(i)
                
                
        # perform the K-S test
        ans = stats.ks_2samp(rvs1, rvs2)
        print 'KS for cosFancyIncList vs all: ',ans
        
        
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
    
        blues = []
        reds = []
        all = allCosFancyInclinations
        
        # remove null "-99" values and split into red and blue groups
        for i,d in zip(cosFancyIncList,difList):
            # check for != -99
            if i>=0:
                # d = vel_galaxy - vel_absorber --> positive = blue shifted absorber (closer to us)
                if d>=0:
                    blues.append(i)
                if d<0:
                    reds.append(i)
    
                
        # perform the K-S test
        ans1 = stats.ks_2samp(blues, reds)
        ans1a = stats.anderson_ksamp([blues,reds])
        print 'KS for blue vs red: ',ans1
        print 'AD for blue vs red: ',ans1a
        
        ans2 = stats.ks_2samp(blues, all)
        print 'KS for blue vs all: ',ans2
        
        ans3 = stats.ks_2samp(reds, all)
        print 'KS for red vs all: ',ans3
        
        # plot the distributions 
        fig = figure()
        
        bins = 15
        
        ax1 = fig.add_subplot(311)
        plot1 = hist(blues,bins=bins)
        title('blueshifted Cos(fancy_inc)')
        
        ax2 = fig.add_subplot(312)
        plot2 = hist(reds,bins=bins)
        title('redshifted Cos(fancy_inc)')
        
        ax3 = fig.add_subplot(313)
        plot3 = hist(all,bins=bins)
        title('Full galaxy table Cos(fancy_inc)')
        
        show()



########################################################################################
########################################################################################

    # plot histograms of the inclinations for both associated galaxies and the 
    # full galaxy data set, combining both redshifted and blueshifted
    plot_dist_cosInc_red_blue = False
    
    if plot_dist_cosInc_red_blue:
    
        blues = []
        reds = []
        all = allCosInclinations
        
        # remove null "-99" values and split into red and blue groups
        for i,d in zip(cosIncList,difList):
            # check for != -99
            if i>=0:
                # d = vel_galaxy - vel_absorber --> positive = blue shifted absorber (closer to us)
                if d>=0:
                    blues.append(i)
                if d<0:
                    reds.append(i)
                
                
        # perform the K-S test
        ans1 = stats.ks_2samp(blues, reds)
        ans1a = stats.anderson_ksamp([blues,reds])
        print 'KS for blue vs red: ',ans1
        print 'AD for blue vs red: ',ans1a
        
        ans2 = stats.ks_2samp(blues, all)
        print 'KS for blue vs all: ',ans2
        
        ans3 = stats.ks_2samp(reds, all)
        print 'KS for red vs all: ',ans3
        
        # plot the distributions 
        fig = figure(figsize=(8,8))
        subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.4)

        bins = 15
        
        ax1 = fig.add_subplot(311)
        plot1 = hist(blues,bins=bins)
        title('blueshifted Cos(inc)')
        
        ax2 = fig.add_subplot(312)
        plot2 = hist(reds,bins=bins)
        title('redshifted Cos(inc)')
        
        ax3 = fig.add_subplot(313)
        plot3 = hist(all,bins=bins)
        title('Full galaxy table Cos(inc)')
        
        show()


########################################################################################
########################################################################################

    # plot histograms of the inclinations for both associated galaxies and the 
    # full galaxy data set, combining both redshifted and blueshifted
    plot_dist_fancy_inc_red_blue = True
    
    if plot_dist_fancy_inc_red_blue:
    
        blues = []
        reds = []
        all = allFancyInclinations
        
        # remove null "-99" values and split into red and blue groups
        for i,d in zip(fancyIncList,difList):
            # check for != -99
            if i>=0:
                # d = vel_galaxy - vel_absorber --> positive = blue shifted absorber (closer to us)
                if d>=0:
                    blues.append(i)
                if d<0:
                    reds.append(i)
                
                
        # perform the K-S test
        ans1 = stats.ks_2samp(blues, reds)
        ans1a = stats.anderson_ksamp([blues,reds])
        print 'KS for blue vs red: ',ans1
        print 'AD for blue vs red: ',ans1a
        
        ans2 = stats.ks_2samp(blues, all)
        print 'KS for blue vs all: ',ans2
        
        ans3 = stats.ks_2samp(reds, all)
        print 'KS for red vs all: ',ans3
        
        # plot the distributions 
        fig = figure(figsize=(8,8))
        subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.4)

        bins = 15
        
        ax1 = fig.add_subplot(311)
        plot1 = hist(blues,bins=bins)
        title('blueshifted fancy_inc')
        
        ax2 = fig.add_subplot(312)
        plot2 = hist(reds,bins=bins)
        title('redshifted fancy_inc')
        
        ax3 = fig.add_subplot(313)
        plot3 = hist(all,bins=bins)
        title('Full galaxy table fancy_inc')
        
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
    