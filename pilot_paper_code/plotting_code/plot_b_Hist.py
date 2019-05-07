#!/usr/bin/env python

'''
By David French (frenchd@astro.wisc.edu)

$Id:  plot_b_Hist.py, v 1.1 04/21/2016

Plots histograms of dopplar b parameters for red and blue shifted absorbers (03/02/2016)

v1.1: remake plots with v_hel instead of vcorr (4/21/16)


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



###########################################################################

    
def main():
    # assuming 'theFile' contains one name per line, read the file
    
    
    if getpass.getuser() == 'David':
#         resultsFilename = '/Users/David/Research_Documents/inclination/git_inclination/LG_correlation_combined5_8_edit2.csv'
#         saveDirectory = '/Users/David/Research_Documents/inclination/git_inclination/pilot_paper_code/plots2/'
        resultsFilename = '/Users/David/Research_Documents/inclination/git_inclination/LG_correlation_combined5_9_edit2.csv'
        saveDirectory = '/Users/David/Research_Documents/inclination/git_inclination/pilot_paper_code/plots3/'

    elif getpass.getuser() == 'frenchd':
        resultsFilename = '/usr/users/frenchd/inclination/git_inclination/LG_correlation_combined5_9_edit2.csv'
        saveDirectory = '/usr/users/frenchd/inclination/git_inclination/pilot_paper_code/plots3/'

    else:
        print 'Could not determine username. Exiting.'
        sys.exit()
    
    
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
    lyaErrList = []
    naList = []
    bList = []
    impactList = []
    azList = []
    incList = []
    cosIncList = []
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
            minor = l['minorAxis (kpc)']
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
                
                if isNumber(maj) and isNumber(minor):
                    q0 = 0.2
                    fancyInc = calculateFancyInclination(maj,minor,q0)
                    cosFancyInc = cos(fancyInc * pi/180)
                else:
                    fancyInc = -99
                    cosFancyInc = -99
            else:
                cosInc = -99
                inc = -99
                fancyInc = -99
                cosFancyInc = -99
                
            if isNumber(az):
                az = float(az)
            else:
                az = -99
            
            # all the lists to be used for associated lines
            lyaVList.append(float(lyaV))
            lyaWList.append(float(lyaW))
            lyaErrList.append(float(lyaW_err))
            naList.append(na)
            bList.append(float(b))
            impactList.append(float(impact))
            azList.append(az)
            incList.append(float(inc))
            cosIncList.append(cosInc)
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
    
    total = 0
    totalNo = 0
    totalYes = 0
    totalIsolated = 0
    totalGroup = 0
    

##########################################################################################
##########################################################################################
    # make a histogram of all dopplar b parameter
    #
    
    plot_b_Hist = False
    save = False
    
    if plot_b_Hist:
        fig = figure()
        ax = fig.add_subplot(111)

#         bins = arange(0,100,10)
        print 'bList: ',bList
        print
        
        binsize = 5
        bins = arange(min(bList),max(bList)+binsize,binsize)
        print 'bins: ',bins
        
        plot1 = hist(bList,bins=bins,histtype='bar')

        xlabel('Azimuth (deg)')
        ylabel('Number')
#         xlim(0,90)
#         ylim(0,10)
        
        if save:
            savefig('{0}/hist(b).pdf'.format(saveDirectory),format='pdf')
        else:
            show()
      
#########################################################################################
#########################################################################################
    # make histograms for red and blue shifted dopplar b parameters
    #
    
    plot_b_Hist_dif = False
    save = False
    
    if plot_b_Hist_dif:
    
        fig = figure(figsize=(10,6))
        subplots_adjust(hspace=0.200)

        binsize = 5
        bins = arange(min(bList),max(bList)+binsize,binsize)
        alpha = 0.85

        blue = []
        red = []
        blueLya = []
        blueLyaErr = []
        redLya = []
        redLyaErr = []
        
        for d,b,l,e in zip(difList,bList,lyaWList,lyaErrList):
            if b != -99 and isNumber(b):
                if d >=0:
                    # blue shifted absorber, but galaxy is REDSHIFTED
                    blue.append(b)
                    blueLya.append(l)
                    blueLyaErr.append(e)
                else:
                    # red shifted absorber, but galaxy is BLUESHIFTED
                    red.append(b)
                    redLya.append(l)
                    redLyaErr.append(e)
        
        print '--------stats--------'
        print 'max red: ',max(red)
        print 'min red: ',min(red)
        print 'max blue:' ,max(blue)
        print 'min blue: ',min(blue)
        
        ax = fig.add_subplot(211)        
        hist(red,bins=bins,histtype='bar',color='red',alpha = alpha,label='Redshifted absorbers')
#         title('Red shifted absorption: Galaxies')    
        ylabel("Number")
#         ylim(0,7)
#         xlim(0,90)
        legend()

        ax = fig.add_subplot(212)
        hist(blue,bins=bins,histtype='bar',color='Blue',alpha = alpha,label='Blueshifted absorbers')
#         title('Blue shifted absorption: Galaxies')
        ylabel('Number')
        xlabel(r'$b$ [km/s]')
#         xlim(0,90)
#         ylim(0,7)
        legend()

#         tight_layout()

        if save:
            savefig('{0}/hist(b)_dif_{1}.pdf'.format(saveDirectory,binsize),format='pdf')
        else:
            show()

#########################################################################################
#########################################################################################
    # make histograms for red and blue shifted dopplar b-parameters overlaid on each other
    #
    
    plot_b_Hist_over_dif = True
    save = False
    
    if plot_b_Hist_over_dif:
    
        fig = figure(figsize=(10,6))
        subplots_adjust(hspace=0.200)

        binsize = 5
        bins = arange(min(bList),max(bList)+binsize,binsize)
        alpha = 0.65
        
        blue = []
        red = []
        blueLya = []
        blueLyaErr = []
        redLya = []
        redLyaErr = []
        
        for d,b,l,e in zip(difList,bList,lyaWList,lyaErrList):
            if b != -99 and isNumber(b):
                if d >=0:
                    # blue shifted absorber, but galaxy is REDSHIFTED
                    blue.append(b)
                    blueLya.append(l)
                    blueLyaErr.append(e)
                else:
                    # red shifted absorber, but galaxy is BLUESHIFTED
                    red.append(b)
                    redLya.append(l)
                    redLyaErr.append(e)
        
        ax = fig.add_subplot(111)        
        hist(red,bins=bins,histtype='bar',color='red',alpha = alpha,hatch = '/',label='Redshifted absorbers')
        hist(blue,bins=bins,histtype='bar',color='Blue',alpha = alpha,label='Blueshifted absorbers')

        ylabel('Number')
        xlabel("Dopplar b-parameter (km/s)")
#         xlim(0,90)
#         ylim(0,7)
        legend()

#         tight_layout()

        if save:
            savefig('{0}/hist(b)_overlaid_dif.pdf'.format(saveDirectory),format='pdf')
        else:
            show()


##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################


if __name__=="__main__":
    # do the work
    main()
    