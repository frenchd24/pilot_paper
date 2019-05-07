#!/usr/bin/env python

'''
By David French (frenchd@astro.wisc.edu)

$Id:  plotAzHist2.py, v 5.5 9/26/16

This is the plotAzHist bit from histograms3.py. Now is separated, and loads in a pickle
file of the relevant data, as created by "buildDataLists.py"

Plots a histogram of azimuth angles for red and blue shifted absorbers

Previous (from histograms3.py):
    Plot some stuff for the 100largest initial results

    Make plots for AAS winter 2014 poster. Uses LG_correlation_combined2.csv file

    Updated for the pilot paper (05/06/15)

v5: updated to work with the new, automatically updated LG_correlation_combined5.csv
    (12/04/15) more on (12/28/2015)
    
v5.1: updated to work with LG_correlation_combined5_8_edit2.csv and l_min = 0.001
    (02/22/2016)

v5.2: remake plots with v_hel instead of vcorr (4/21/16)

v5.3: remake plots with newest large galaxy sample (07/13/16) -> /plots4/

v5.4: include ability to limit results based on environment number (7/14/16)

v5.5: update for LG_correlation_combined5_11_25cut_edit4.csv -> /plots5/ (9/26/16)

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

# from matplotlib import rc
# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
# ## for Palatino and other serif fonts use:
# #rc('font',**{'family':'serif','serif':['Palatino']})
# rc('text', usetex=True)

from matplotlib import rc
# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})

# plt.clf()
# plt.rcdefaults()
# fontScale = 18
# params = {'axes.labelsize': fontScale,
# 'axes.titlesize': fontScale,
# 'font.size': fontScale,
# 'xtick.labelsize': fontScale,
# 'ytick.labelsize': fontScale,
# 'font.weight': 450,
# 'axes.labelweight': 450,
# 'xtick.major.size': 4,
# 'xtick.major.width': 1.2,
# 'xtick.minor.size': 3,
# 'xtick.minor.width': 1.2,
# 'ytick.major.size': 4,
# 'ytick.major.width': 1.2,
# 'ytick.minor.size': 3,
# 'ytick.minor.width': 1.2
# }


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
#         resultsFilename = '/Users/David/Research_Documents/inclination/git_inclination/LG_correlation_combined5_8_edit2.csv'
#         saveDirectory = '/Users/David/Research_Documents/inclination/git_inclination/pilot_paper_code/plots2/'
#         resultsFilename = '/Users/David/Research_Documents/inclination/git_inclination/LG_correlation_combined5_9_edit2.csv'
        resultsFilename = '/Users/David/Research_Documents/inclination/git_inclination/LG_correlation_combined5_11_25cut_edit4.csv'
        saveDirectory = '/Users/David/Research_Documents/inclination/git_inclination/pilot_paper_code/plots5/'

    elif getpass.getuser() == 'frenchd':
#         resultsFilename = '/usr/users/frenchd/inclination/git_inclination/LG_correlation_combined5_8_edit2.csv'
#         saveDirectory = '/usr/users/frenchd/inclination/git_inclination/pilot_paper_code/plots2/'
#         resultsFilename = '/usr/users/frenchd/inclination/git_inclination/LG_correlation_combined5_9_edit2.csv'
        resultsFilename = '/usr/users/frenchd/inclination/git_inclination/LG_correlation_combined5_11_25cut_edit4.csv'
        saveDirectory = '/usr/users/frenchd/inclination/git_inclination/pilot_paper_code/plots5/'

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
    
    maxEnv = 3000
    
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
            if float(env) <= maxEnv:
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
    
#     for a,n,g in zip(azList,newAzList,galaxyNameList):
#         if a != -99 and n !=-99:
#             if abs(a-n) >5:
#                 print "{0} : old = {1}, new = {2}".format(g,a,n)
                
    

##########################################################################################
##########################################################################################
    # make a histogram of all the azimuth angles
    #
    
    plotAzHist = False
    save = False
    
    #font = {'family':'serif','size':16}
#     font = {'family':'serif','size':18, 'serif': ['computer modern roman']}
#     rc('font',**font)
#     rc('legend',**{'fontsize':18})
#     rcParams['text.latex.preamble']=[r'\usepackage{amsmath}']
    
    
    if plotAzHist:
        fig = figure(figsize=(10,8))
        ax = fig.add_subplot(111)
        alpha = 0.85
        
#         bins = [0,10,20,30,40,50,60,70,80,90]
#         bins = array([0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90])
#         bins +=4
#         bins = 45

        bins = arange(0,100,10)
        print 'bins:' ,bins
        print 'stdev: ',std(azList)
        
        bot = 0
        top = 45
        count = 0
        for a in azList:
            if a <=top and a>=bot:
                count+=1
                
        print 'between {0} and {1} = {2}'.format(top,bot,count)
        print 'total number: ',len(azList)
        print 'ratio = ',float(count)/len(azList)
        
        print 'mean: ',mean(azList)
        print 'median: ',median(azList)
        
#         rc('text', usetex=True)
#         rc('font', family='serif',weight='bold')
        
        plot1 = hist(azList,bins=bins,histtype='bar',alpha=alpha)
#         ax.tick_params(axis='both', which='major', labelsize=20)
        
#         print 'azList: ',azList
#         hist(azList)
#         title('Distribution of azimuths')

#         xlabel(r'Azimuth (deg)',fontsize=20,weight='bold')
#         ylabel(r'Number',fontsize=20,weight='bold')
        xlabel(r'Azimuth (deg)')
        ylabel(r'Number')
        xlim(0,90)
        ylim(0,10)
        tight_layout()
        
        if save:
            savefig('{0}/hist(azimuth).pdf'.format(saveDirectory),format='pdf')
        else:
            show()
      
#########################################################################################
#########################################################################################
    # make histograms for red and blue shifted azimuth angles
    #
    
    plotAzHist_dif = False
    save = False
    
    if plotAzHist_dif:
    
        fig = figure(figsize=(10,6))
        subplots_adjust(hspace=0.200)

#         bins = [0,15,30,45,60,75,90]
#         bins = arange(0,90,10)
        binsize = 10
        bins = arange(0,100,binsize)
        blue = []
        red = []
        alpha = 0.80
        
        blueLya = []
        blueLyaErr = []

        redLya = []
        redLyaErr = []
        
        for d,a,l,e in zip(difList,azList,lyaWList,lyaErrList):
            if a != -99:
                if d >=0:
                    # blue shifted absorber, but galaxy is REDSHIFTED
                    print 'd: ',d
                    blue.append(a)
                    blueLya.append(l)
                    blueLyaErr.append(e)
                else:
                    # red shifted absorber, but galaxy is BLUESHIFTED
                    red.append(a)
                    redLya.append(l)
                    redLyaErr.append(e)
        
        print 'stats-----'
        print 'max red: ',max(red)
        print 'min red: ',min(red)
        print 'max blue:' ,max(blue)
        print 'min blue: ',min(blue)
        
        ax = fig.add_subplot(211)        
        hist(red,bins=bins,histtype='bar',color='red',alpha = alpha,label='Redshifted absorbers')
#         title('Red shifted absorption: Galaxies')    
        ylabel("Number")
        ylim(0,7)
        xlim(0,90)
        legend()

        ax = fig.add_subplot(212)
        hist(blue,bins=bins,histtype='bar',color='Blue',alpha = alpha,label='Blueshifted absorbers')
#         title('Blue shifted absorption: Galaxies')
        ylabel('Number')
        xlabel("Azimuth (deg)")
        xlim(0,90)
        ylim(0,7)
        legend()

#         tight_layout()

        if save:
            savefig('{0}/hist(azimuth)_dif_{1}.pdf'.format(saveDirectory,binsize),format='pdf')
        else:
            show()

#########################################################################################
#########################################################################################
    # make histograms for red and blue shifted azimuth angles overlaid on each other
    #
    
    plotAzHist_over_dif = False
    save = False
    
    if plotAzHist_over_dif:
    
        fig = figure(figsize=(10,6))
        subplots_adjust(hspace=0.200)
        
        alpha = 0.75

        bins = arange(0,100,10)
        blue = []
        red = []
        
        blueLya = []
        blueLyaErr = []

        redLya = []
        redLyaErr = []
        
        for d,a,l,e in zip(difList,azList,lyaWList,lyaErrList):
            if a != -99:
                if d >=0:
                    # blue shifted absorber, but galaxy is REDSHIFTED
                    print 'd: ',d
                    blue.append(a)
                    blueLya.append(l)
                    blueLyaErr.append(e)
                else:
                    # red shifted absorber, but galaxy is BLUESHIFTED
                    red.append(a)
                    redLya.append(l)
                    redLyaErr.append(e)
        
        ax = fig.add_subplot(111)        
        hist(red,bins=bins,histtype='bar',color='red',alpha = alpha,hatch = '/',label='Redshifted absorbers')
        hist(blue,bins=bins,histtype='bar',color='Blue',alpha = alpha,label='Blueshifted absorbers')

        ylabel('Number')
        xlabel("Azimuth (deg)")
        xlim(0,90)
        ylim(0,7)
        legend()

#         tight_layout()

        if save:
            savefig('{0}/hist(azimuth)_overlaid_dif.pdf'.format(saveDirectory),format='pdf')
        else:
            show()
            
            
#########################################################################################
#########################################################################################
    # make histograms for all azimuth angles, with red and blue shifted overlaid
    #
    
    plotAzHist_all_over = True
    save = True
    
    if plotAzHist_all_over:
    
        fig = figure(figsize=(10,5))
#         subplots_adjust(hspace=0.200)
        
        alpha = 0.6

        bins = arange(0,100,10)
        blue = []
        red = []
        
        blueLya = []
        blueLyaErr = []

        redLya = []
        redLyaErr = []
        
        for d,a,l,e in zip(difList,azList,lyaWList,lyaErrList):
            if a != -99:
                if d >=0:
                    # blue shifted absorber, but galaxy is REDSHIFTED
                    print 'd: ',d
                    blue.append(a)
                    blueLya.append(l)
                    blueLyaErr.append(e)
                else:
                    # red shifted absorber, but galaxy is BLUESHIFTED
                    red.append(a)
                    redLya.append(l)
                    redLyaErr.append(e)
        
        
        ax = fig.add_subplot(111)
        
        # all first
        hist(blue,bins=bins,histtype='bar',color='Blue',alpha=alpha+0.25,lw=1.7,hatch='//',label=r'$\rm Blueshifted$')
        hist(red,bins=bins,histtype='bar',color='red',alpha=alpha,lw=1.7,label=r'$\rm Redshifted$')
        plot1 = hist(azList,bins=bins,histtype='step',lw=2.1,alpha=0.9,color='black',label=r'$\rm All$')


        # x-axis
        majorLocator   = MultipleLocator(10)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(5)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        
        # y axis
        majorLocator   = MultipleLocator(2)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(1)
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)

        ylabel(r'$\rm Number$')
        xlabel(r'$\rm Azimuth ~ [deg]$')
        xlim(0,90)
        ylim(0,10)
        legend(scatterpoints=1,prop={'size':16},loc=2,fancybox=True)


#         tight_layout()

        if save:
            savefig('{0}/hist(azimuth)_overlaid_all.pdf'.format(saveDirectory),format='pdf',bbox_inches='tight')
        else:
            show()

#########################################################################################
#########################################################################################
    # make histograms for red and blue shifted azimuth angles vs flat distributions
    #
    
    plotAzHist_dif_flat = False
    save = False
    
    if plotAzHist_dif_flat:
    
        fig = figure(figsize=(10,6))
        subplots_adjust(hspace=0.200)

#         bins = [0,15,30,45,60,75,90]
#         bins = arange(0,90,10)
        bins = arange(0,90,12)
        blue = []
        red = []
        
        blueLya = []
        blueLyaErr = []

        redLya = []
        redLyaErr = []
        
        for d,a,l,e in zip(difList,azList,lyaWList,lyaErrList):
            if a != -99:
                if d >=0:
                    # blue shifted absorber, but galaxy is REDSHIFTED
                    blue.append(a)
                    blueLya.append(l)
                    blueLyaErr.append(e)
                else:
                    # red shifted absorber, but galaxy is BLUESHIFTED
                    red.append(a)
                    redLya.append(l)
                    redLyaErr.append(e)
        
        
#         flatRedNum = len(red)/9.0
#         print 'flatRedNum: ',flatRedNum
#         flatRed = ones(9)*flatRedNum
#     
#         flatBlueNum = len(blue)/9.0
#         flatBlue = ones(9)*flatBlueNum

        flatRed = arange(0,90,5)
        flatBlue = arange(0,90,5)
        
        print 'flatRed: ',flatRed
        print
        print 'flatBlue: ',flatBlue
        print
        
        ax = fig.add_subplot(411)        
        hist(red,bins=bins,histtype='bar',color='red',alpha = 0.80,label='Redshifted absorbers')
#         title('Red shifted absorption: Galaxies')    
        ylabel("Number")
        ylim(0,6)
        xlim(0,90)
        legend()
        
        ax = fig.add_subplot(412)        
        hist(flatRed,bins=bins,histtype='bar',color='red',alpha = 0.80,label='Flat Redshifted absorbers')
#         title('Red shifted absorption: Galaxies')    
        ylabel("Number")
        ylim(0,6)
        xlim(0,90)
        legend()

        ax = fig.add_subplot(413)
        hist(blue,bins=bins,histtype='bar',color='Blue',alpha = 0.80,label='Blueshifted absorbers')
#         title('Blue shifted absorption: Galaxies')
        ylabel('Number')
        xlabel("Azimuth (deg)")
        xlim(0,90)
        ylim(0,6)
        legend()

        ax = fig.add_subplot(414)
        hist(flatBlue,bins=bins,histtype='bar',color='Blue',alpha = 0.80,label='Flat Blueshifted absorbers')
#         title('Blue shifted absorption: Galaxies')
        ylabel('Number')
        xlabel("Azimuth (deg)")
        xlim(0,90)
        ylim(0,6)
        legend()
        
#         tight_layout()

        if save:
            savefig('{0}/hist(azimuth)_dif_flat.pdf'.format(saveDirectory),format='pdf')
        else:
            show()

##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################


if __name__=="__main__":
    # do the work
    main()
    