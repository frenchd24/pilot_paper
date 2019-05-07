#!/usr/bin/env python

'''
By David French (frenchd@astro.wisc.edu)

$Id:  plotHist_randomDist.py, v 1.0 09/22/16

Plot histograms of the inclinations of the associated galaxies, all galaxies, and a 
random distribution of inclinations.


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
import random

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

fontScale = 24
rc('text', usetex=True)
rc('font', size=24, family='serif', weight='normal')
rc('xtick.major',size=8,width=0.6)
rc('xtick.minor',size=5,width=0.6)
rc('ytick.major',size=8,width=0.6)
rc('ytick.minor',size=5,width=0.6)
rc('xtick',labelsize = fontScale)
rc('ytick',labelsize = fontScale)
rc('axes',labelsize = fontScale)
rc('xtick', labelsize = fontScale)
rc('ytick',labelsize = fontScale)
# rc('font', weight = 450)
# rc('axes',labelweight = 'bold')
rc('axes',linewidth = 1)


###########################################################################

    
def main():
    # assuming 'theFile' contains one name per line, read the file
    
    
    if getpass.getuser() == 'David':
        pickleFilename = '/Users/David/Research_Documents/inclination/git_inclination/pilot_paper_code/pilotData2.p'
#         resultsFilename = '/Users/David/Research_Documents/inclination/git_inclination/LG_correlation_combined5_8_edit2.csv'
#         saveDirectory = '/Users/David/Research_Documents/inclination/git_inclination/pilot_paper_code/plots2/'
        resultsFilename = '/Users/David/Research_Documents/inclination/git_inclination/LG_correlation_combined5_11_25cut_edit4.csv'
        saveDirectory = '/Users/David/Research_Documents/inclination/git_inclination/pilot_paper_code/plots5/'

    elif getpass.getuser() == 'frenchd':
        pickleFilename = '/usr/users/frenchd/inclination/git_inclination/pilot_paper_code/pilotData2.p'
#         resultsFilename = '/usr/users/frenchd/inclination/git_inclination/LG_correlation_combined5_8_edit2.csv'
#         saveDirectory = '/usr/users/frenchd/inclination/git_inclination/pilot_paper_code/plots2/'
        resultsFilename = '/usr/users/frenchd/inclination/git_inclination/LG_correlation_combined5_11_25cut_edit4.csv'
        saveDirectory = '/usr/users/frenchd/inclination/git_inclination/pilot_paper_code/plots5/'

    else:
        print 'Could not determine username. Exiting.'
        sys.exit()
    
    # use the old pickle file to get the full galaxy dataset info
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
    
    maxEnv = 2000
    minL = 0.001
    
    # if match, then the includes in the file have to MATCH the includes above. e.g., if 
    # virInclude = False, cusInclude = True, finalInclude = False, then only systems
    # matching those three would be included. Otherwise, all cusInclude = True would be included
    # regardless of the others
    match = False
    
    # all the lists to be used for associated lines
    lyaVList = []
    lyaWList = []
    lyaErrList = []
    naList = []
    bList = []
    impactList = []
    azList = []
    incList = []
    fancyIncList = []
    spiralFancyIncList = []
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
            morph = l['final_morphology']
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
            if float(env) <= maxEnv and float(likelihood) >=minL:
                lyaVList.append(float(lyaV))
                lyaWList.append(float(lyaW))
                lyaErrList.append(float(lyaW_err))
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
                
#                 lmorph = morph.lower()
#                 if bfind(morph,'Spiral'):
#                     spiralFancyIncList.append(fancyInc)

    results.close()
    
        
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
    # fancyInclination histograms for redshifted vs blueshifted distributions of absorbers
    # plotted together and separately all in one, then the whole table as a subplot
    #
    
    plotFancyIncDifHist_full_all_random = True
    save = False
    
    if plotFancyIncDifHist_full_all_random:
    
        galaxyFile = open('/Users/David/Research_Documents/gt/NewGalaxyTable5.csv','rU')
        reader = csv.DictReader(galaxyFile)
        
        allDiameters = []
        incGT25diam = []
        spiralIncList = []
        
        
        q0 = 0.2
        for i in reader:
            major,minor = eval(i['linDiameters (kpc)'])
            morph = i['morphology'].lower()
            if bfind(morph,'s'):
                if not bfind(morph,'sph') and not bfind(morph,'s0'):
                
                    if isNumber(major):
                        if isNumber(minor):
                            if float(major) > float(minor):
                                fInc = calculateFancyInclination(major,minor,q0)
                                spiralIncList.append(fInc)
                                
                                if float(major) >=25.0:
                                    incGT25diam.append(fInc)
        
        galaxyFile.close()
        
    #########################################################################################
        # build a random sample of galaxies
        randMinorList = []
        randIncList = []
        randCosIncList = []
        randQList = []
        size = len(allFancyInclinations)
        minMin = 0.2
        maxMin = 1.0
        q0=0.2
        for i in range(0,size):
            minor = random.uniform(minMin,maxMin)
            randMinorList.append(minor)
            randQList.append(minor/maxMin)
            
            inc = calculateFancyInclination(maxMin,minor,q0)
            cosInc = math.cos(calculateFancyInclination(maxMin,minor,q0)*(math.pi/180.0))
    
            randIncList.append(inc)
            randCosIncList.append(cosInc)
            
        n, bins, patches = hist(randQList, bins = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
        print 'n: ',n
        print 'bins: ',bins
    #########################################################################################

    
        fig = figure(figsize=(10,6))
        subplots_adjust(hspace=0.210)

#         bins = [0,.10,.20,.30,.40,.50,.60,.70,.80,.90]
        bins = arange(0,100,10)
        blue = []
        red = []
        
        alpha = 0.6
        
        blueLya = []
        blueLyaErr = []

        redLya = []
        redLyaErr = []
        
        totalIncList = []
        
        for d,i,l,e,m in zip(difList,fancyIncList,lyaWList,lyaErrList,morphList):
            totalIncList.append(i)
            if d >=0:
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
                
                
        # just associated
        ax = fig.add_subplot(211)
        
        # x-axis
        majorLocator   = MultipleLocator(10)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(5)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        
        # y-axis
        majorLocator   = MultipleLocator(5)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(2.5)
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)

#         hist(red,bins=bins,histtype='bar',color='red',hatch='\\',lw=1.5,alpha = alpha,label=r'$\rm Redshifted$')        
#         hist(blue,bins=bins,histtype='bar',color='Blue',lw=1.5,alpha = alpha,label=r'$\rm Blueshifted$')
        hist(blue,bins=bins,histtype='bar',color='Blue',hatch='\\',lw=1.7,alpha = alpha+0.25,label=r'$\rm Blueshifted$')
        hist(red,bins=bins,histtype='bar',color='red',lw=1.7,alpha = alpha,label=r'$\rm Redshifted$')        
        hist(totalIncList,bins=bins,histtype='step',color='Black',lw=2.1,alpha = 0.9,label=r'$\rm All ~ Associated$')
        
        legend(scatterpoints=1,prop={'size':16},loc=2,fancybox=True)
        ylabel(r'$\rm Number$')
        

        # full table
        ax = fig.add_subplot(212)
        
        # x-axis
        majorLocator   = MultipleLocator(10)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(5)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        
        # y-axis
        majorLocator   = MultipleLocator(5000)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(2500)
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)
        
        hist(allFancyInclinations,bins=bins,histtype='bar',lw=1.5,color = 'green',alpha=0.85,label=r'$\rm Observed$')
#         hist(fIncList,bins=bins,histtype='bar',lw=1.5,color = 'green',alpha=0.85,label=r'$\rm Observed$')
        hist(randIncList,bins=bins,histtype='step',lw=1.5,ls='dashed',color = 'black',alpha=0.95,label=r'$\rm Random$')
        legend(scatterpoints=1,prop={'size':16},loc=2,fancybox=True)
        xlabel(r'$\rm Galaxy ~ Inclination ~ [deg]$')
        ylabel(r'$\rm Number$')
#         tight_layout()

        if save:
            savefig('{0}/hist(fancy_inclination)_red_blue_all_random2.pdf'.format(saveDirectory),format='pdf',bbox_inches='tight')
        else:
            show()


    
#########################################################################################
#########################################################################################
    # fancyInclination histograms for redshifted vs blueshifted distributions of absorbers
    # plotted together and separately all in one, then the whole table as a subplot
    #
    # ONLY SPIRAL GALAXIES
    
    plotFancyIncDifHist_full_all_random_spiral = False
    save = False
    
    if plotFancyIncDifHist_full_all_random_spiral:
    
        galaxyFile = open('/Users/David/Research_Documents/gt/NewGalaxyTable5.csv','rU')
        reader = csv.DictReader(galaxyFile)
        
        allDiameters = []
        incGT25diam = []
        spiralIncList = []
        
        
        q0 = 0.2
        for i in reader:
            major,minor = eval(i['linDiameters (kpc)'])
            morph = i['morphology'].lower()
            if bfind(morph,'s'):
                if not bfind(morph,'sph') and not bfind(morph,'s0'):
                
                    if isNumber(major):
                        if isNumber(minor):
                            if float(major) > float(minor):
                                fInc = calculateFancyInclination(major,minor,q0)
                                spiralIncList.append(fInc)
                                
                                if float(major) >=25.0:
                                    incGT25diam.append(fInc)
        
        galaxyFile.close()
        
    #########################################################################################
        # build a random sample of galaxies
        randMinorList = []
        randIncList = []
        randAcosIncList = []
        size = len(spiralIncList)
        minMin = 0.2
        maxMin = 1.0
        q0=0.2
        for i in range(0,size):
            minor = random.uniform(minMin,maxMin)
            randMinorList.append(minor)
        
            inc = calculateFancyInclination(maxMin,minor,q0)
            acosInc = calculateInclination(maxMin,minor)
    
            randIncList.append(inc)
            randAcosIncList.append(acosInc)
    #########################################################################################

    
        fig = figure(figsize=(10,6))
        subplots_adjust(hspace=0.200)

#         bins = [0,.10,.20,.30,.40,.50,.60,.70,.80,.90]
        bins = arange(0,100,10)
        blue = []
        red = []
        
        alpha = 0.6
        
        blueLya = []
        blueLyaErr = []

        redLya = []
        redLyaErr = []
        
        totalIncList = []
        
        for d,i,l,e,m in zip(difList,fancyIncList,lyaWList,lyaErrList,morphList):
            if m == 'Spiral':
                totalIncList.append(i)
                if d >=0:
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
                
                
        # just associated
        ax = fig.add_subplot(211)
        
        # x-axis
        majorLocator   = MultipleLocator(10)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(5)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        
        # y-axis
        majorLocator   = MultipleLocator(5)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(2.5)
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)

#         hist(red,bins=bins,histtype='bar',color='red',hatch='\\',lw=1.5,alpha = alpha,label=r'$\rm Redshifted$')        
#         hist(blue,bins=bins,histtype='bar',color='Blue',lw=1.5,alpha = alpha,label=r'$\rm Blueshifted$')
        hist(blue,bins=bins,histtype='bar',color='Blue',hatch='\\',lw=1.7,alpha = alpha+0.25,label=r'$\rm Blueshifted$')
        hist(red,bins=bins,histtype='bar',color='red',lw=1.7,alpha = alpha,label=r'$\rm Redshifted$')        
        hist(totalIncList,bins=bins,histtype='step',color='Black',lw=2.1,alpha = 0.9,label=r'$\rm All ~ Associated$')
        
        legend(scatterpoints=1,prop={'size':16},loc=2,fancybox=True)
        ylabel(r'$\rm Number$')
        

        # full table
        ax = fig.add_subplot(212)
        
        # x-axis
        majorLocator   = MultipleLocator(10)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(5)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        
        # y-axis
        majorLocator   = MultipleLocator(5000)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(2500)
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)
        
#         hist(allFancyInclinations,bins=bins,histtype='bar',lw=1.5,color = 'green',alpha=0.9,label=r'$\rm Observed$')
        hist(spiralIncList,bins=bins,histtype='bar',lw=1.5,color = 'green',alpha=0.85,label=r'$\rm Observed$')
        hist(randIncList,bins=bins,histtype='step',lw=1.5,ls='dashed',color = 'black',alpha=0.95,label=r'$\rm Random$')
        legend(scatterpoints=1,prop={'size':16},loc=2,fancybox=True)
        xlabel(r'$\rm Galaxy ~ Inclination ~ [deg]$')
        ylabel(r'$\rm Number$')
#         tight_layout()

        if save:
            savefig('{0}/hist(fancy_inclination)_red_blue_all_random_spiral.pdf'.format(saveDirectory),format='pdf',bbox_inches='tight')
        else:
            show()
            

#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
    # cos(fancyInclination) histograms for redshifted vs blueshifted distributions of absorbers
    # plotted together and separately all in one, then the whole table as a subplot
    #
    
    plotCosFancyIncDifHist_full_all = False
    save = False
    
    if plotCosFancyIncDifHist_full_all:
    
        fig = figure(figsize=(10,6))
        subplots_adjust(hspace=0.200)

#         bins = [0,.10,.20,.30,.40,.50,.60,.70,.80,.90]
        bins = arange(0,1.0,0.1)
        blue = []
        red = []
        
        blueLya = []
        blueLyaErr = []

        redLya = []
        redLyaErr = []
        
        for d,i,l,e in zip(difList,cosFancyIncList,lyaWList,lyaErrList):
            if d >=0:
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
                
                
        # just associated
        ax = fig.add_subplot(211)
        
        # x-axis
        majorLocator   = MultipleLocator(0.1)
        majorFormatter = FormatStrFormatter(r'$\rm %s$')
        minorLocator   = MultipleLocator(0.05)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        
        # y-axis
        majorLocator   = MultipleLocator(5)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(1)
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)

        hist(blue,bins=bins,histtype='bar',color='Blue',lw=1.5,alpha = 0.7,label=r'$\rm Blueshifted$')
        hist(red,bins=bins,histtype='bar',color='red',lw=1.5,alpha = 0.7,label=r'$\rm Redshifted$')        
        hist(cosFancyIncList,bins=bins,histtype='step',color='Black',lw=2.5,alpha = 0.9,label=r'$\rm All ~ Associated$')
        
        legend(scatterpoints=1,prop={'size':14},loc=2,fancybox=True)
        ylabel(r'$\rm Number$')
        

        # full table
        ax = fig.add_subplot(212)
        
        # x-axis
        majorLocator   = MultipleLocator(0.1)
        majorFormatter = FormatStrFormatter(r'$\rm %s$')
        minorLocator   = MultipleLocator(0.05)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        
        # y-axis
        majorLocator   = MultipleLocator(5000)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(1000)
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)
        
        hist(allCosFancyInclinations,bins=bins,histtype='bar',lw=1.5,color = 'green',alpha=0.9,label=r'$\rm All$')
        legend(scatterpoints=1,prop={'size':14},loc=2,fancybox=True)
        xlabel(r'$\rm Cos[inc]$')
        ylabel(r'$\rm Number$')
#         tight_layout()

        if save:
            savefig('{0}/hist(cos(fancy_inclination))_red_blue_full_all.pdf'.format(saveDirectory),format='pdf',bbox_inches='tight')
        else:
            show()
            
            
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
    # fancyInclination histograms for redshifted vs blueshifted distributions of absorbers
    # plotted together and separately all in one, then the whole table as a subplot, but
    # only those galaxies with D>25kpc
    #
    
    plotFancyIncDifHist_full_largeOnly = False
    save = False
    
    if plotFancyIncDifHist_full_largeOnly:
    
        galaxyFile = open('/Users/David/Research_Documents/gt/NewGalaxyTable5.csv','rU')
        reader = csv.DictReader(galaxyFile)
        
        allDiameters = []
        incGT25diam = []
        q0 = 0.2
        for i in reader:
            major,minor = eval(i['linDiameters (kpc)'])
        
            if isNumber(major):
                if isNumber(minor):
                    if float(major) > float(minor):
                        if float(major) >=25.0:
                            fInc = calculateFancyInclination(major,minor,q0)
                            incGT25diam.append(fInc)
        
        galaxyFile.close()
    
        fig = figure(figsize=(10,6))
        subplots_adjust(hspace=0.200)

#         bins = [0,.10,.20,.30,.40,.50,.60,.70,.80,.90]
        bins = arange(0,100,10)
        blue = []
        red = []
        
        alpha = 0.6
        
        blueLya = []
        blueLyaErr = []

        redLya = []
        redLyaErr = []
        
        for d,i,l,e in zip(difList,fancyIncList,lyaWList,lyaErrList):
            if d >=0:
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
                
                
        # just associated
        ax = fig.add_subplot(211)
        
        # x-axis
        majorLocator   = MultipleLocator(10)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(5)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        
        # y-axis
        majorLocator   = MultipleLocator(5)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(2.5)
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)

#         hist(red,bins=bins,histtype='bar',color='red',hatch='\\',lw=1.5,alpha = alpha,label=r'$\rm Redshifted$')        
#         hist(blue,bins=bins,histtype='bar',color='Blue',lw=1.5,alpha = alpha,label=r'$\rm Blueshifted$')
        hist(blue,bins=bins,histtype='bar',color='Blue',hatch='\\',lw=1.7,alpha = alpha+0.25,label=r'$\rm Blueshifted$')
        hist(red,bins=bins,histtype='bar',color='red',lw=1.7,alpha = alpha,label=r'$\rm Redshifted$')        
        hist(fancyIncList,bins=bins,histtype='step',color='Black',lw=2.1,alpha = 0.9,label=r'$\rm All ~ Associated$')
        
        legend(scatterpoints=1,prop={'size':15},loc=2,fancybox=True)
        ylabel(r'$\rm Number$')
        

        # full table
        ax = fig.add_subplot(212)
        
        # x-axis
        majorLocator   = MultipleLocator(10)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(5)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        
        # y-axis
        majorLocator   = MultipleLocator(5000)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(2500)
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)
        
        hist(incGT25diam,bins=bins,histtype='bar',lw=1.5,color = 'green',alpha=0.9,label=r'$\rm All ~(D>25~kpc)$')
        legend(scatterpoints=1,prop={'size':15},loc=2,fancybox=True)
        xlabel(r'$\rm Galaxy ~ Inclination ~ [deg]$')
        ylabel(r'$\rm Number$')
#         tight_layout()

        if save:
            savefig('{0}/hist(fancy_inclination)_red_blue_full_largeOnly.pdf'.format(saveDirectory),format='pdf',bbox_inches='tight')
        else:
            show()
            
            
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################


if __name__=="__main__":
    # do the work
    main()
    