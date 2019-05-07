#!/usr/bin/env python

'''
By David French (frenchd@astro.wisc.edu)

$Id:  plot_multiple_lines.py, v 1.0 08/03/2016

Investigate cases where either multiple sightlines probe a golaxy, or multiple lines in
one sightlines are associated with one galaxy

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

fontScale = 15
rc('text', usetex=True)
# rc('font', size=15, family='serif', weight=450)
rc('font', size=15, family='serif',weight='normal')
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
        resultsFilename = '/Users/David/Research_Documents/inclination/git_inclination/LG_correlation_combined5_11_25cut_edit.csv'
        saveDirectory = '/Users/David/Research_Documents/inclination/git_inclination/pilot_paper_code/plots4/'

    elif getpass.getuser() == 'frenchd':
        pickleFilename = '/usr/users/frenchd/inclination/git_inclination/pilot_paper_code/pilotData2.p'
        resultsFilename = '/usr/users/frenchd/inclination/git_inclination/LG_correlation_combined5_11_25cut_edit.csv'
        saveDirectory = '/usr/users/frenchd/inclination/git_inclination/pilot_paper_code/plots4/'

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
    
    # for ambiguous lines
    lyaVAmbList = []
    lyaWAmbList = []
    envAmbList = []
    
    galaxyDict = {}
    
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
            
            # check if this is one of multiple lines probing the same galaxy
            if galaxyDict.has_key(galaxyName):
                i = galaxyDict[galaxyName]
                i.append((galaxyName,lyaW,vel_diff))
                galaxyDict[galaxyName] = i
            
            else:
                galaxyDict[galaxyName] = [(galaxyName,float(lyaW),float(vel_diff))]
            
        else:
            lyaV = l['Lya_v']
            lyaW = l['Lya_W'].partition('pm')[0]
            lyaW_err = l['Lya_W'].partition('pm')[2]
            env = l['environment']
            
            lyaVAmbList.append(float(lyaV))
            lyaWAmbList.append(float(lyaW))
            envAmbList.append(float(env))

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
    # plot the EW as a function of vel_dif for galaxies with multiple lines associated
    #
    
    plotW_dif_multiple = True
    save = False
    
    if plotW_dif_multiple:
#         fig = figure(figsize=(2,8))
        fig = figure()
        
        # sort out the multiple galaxy lines
        vals = galaxyDict.values()
        keys = galaxyDict.keys()
        multipleLines = []
        
        blues = []
        reds = []
        names = []
        namesR = []
        namesB = []
        i = -1
        
        for v,k in zip(vals,keys):
            # does this galaxy have multiple lines associated with it?
            if len(v) >1:
                # for each galaxy, go through the associated lines and add EW to either
                # the blue or red shifted bin, and the name to the names list
                names.append(k)
                i +=1
                
                for l in v:
                    EW = l[1]
                    vel = float(l[2])
                    if vel >= 0:
                        blues.append(EW)
                        namesB.append(i)
                        print i,' = ',vel
                    else:
                        reds.append(EW)
                        namesR.append(i)
                        print 'red'
                    
        x = range(len(names))
        
        print 'x: ',x
        print ' names: ',names
        print 'blues: ',blues
        print 'reds: ',reds
                
        ax = fig.add_subplot(111)
        plot1 = plot(namesB, blues, lw=0,marker='o',label=r'$\rm Blueshifted$',color="Blue",alpha = 0.85)
        plot2 = plot(namesR, reds, lw=0,marker='o',label=r'$\rm Redshifted$',color="Red",alpha = 0.85)
        
        plt.xticks(x, names, rotation='vertical')

        legend(scatterpoints=1,prop={'size':12},loc=1)
        
#         ylabel(r'$\rm Number$')
#         ylim(0,10)
        
        # X-axis
#         majorLocator   = MultipleLocator(1)
#         majorFormatter = FormatStrFormatter(r'$\rm %d$')
#         minorLocator   = MultipleLocator(1)
# 
#         ax.xaxis.set_major_locator(majorLocator)
#         ax.xaxis.set_major_formatter(majorFormatter)
#         ax.xaxis.set_minor_locator(minorLocator)
        
        # Y-axis
        majorLocator   = MultipleLocator(200)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(100)
        
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)
        
        ylabel(r'$\rm Equivalent ~ Width ~ [m\AA]$')
        tight_layout()
        
        if save:
            savefig('{0}/W(vel_diff_multiple).pdf'.format(saveDirectory),format='pdf')
        else:
            show()
            

    
#########################################################################################
#########################################################################################
    # make a CDF of the distribution of Ly-alpha equivalent widths for the 
    # associated sample, splitting on red vs blue shifted absorbers
    #
    # USE THE BUILT IN CDF FUNCTION INSTEAD
    
    plotWCDF_alt_dif = False
    save = False
    
    if plotWCDF_alt_dif:
#         fig = figure(figsize=(2,8))
        fig = figure()
        ax = fig.add_subplot(111)
        
        lyaWArray = array(lyaWList)
        lyaWAmbArray = array(lyaWAmbList)
        
        envArray = array(envList)
        envAmbArray = array(envAmbList)
        
        bins = arange(0,max(lyaWArray),5)

        
        blues = []
        reds = []
        for d,l in zip(difList,lyaWList):
            if float(d) >0:
                # blueshifted ABSORBER: dif = v_gal - v_absorber
                blues.append(float(l))
            else:
                reds.append(float(l))
        
        n_reds, bins_reds, patches_reds = hist(reds, bins, normed=1, histtype="step",\
        cumulative=True,color='red',lw=1,label=r'$\rm Redshifted ~CDF$')
        n_blues, bins_bluess, patches_bluess = hist(blues, bins, normed=1, histtype="step",\
        cumulative=True,lw=1,label=r'$\rm Blueshifted ~CDF$',color='blue')

    
#         # plot the y-values against the sorted data
#         plt.plot(sorted_reds,yvals_reds,color='red',lw=3,label=r'$\rm Redshifted ~CDF$')
#         plt.plot(sorted_blues,yvals_blues,color='blue',lw=3,label=r'$\rm Blueshifted ~CDF$')
        
        xlabel(r'$\rm Equivalent ~ Width ~ [m\AA]$')
        ylabel(r'$\rm Number$')
        
        # format all the axis and stuff
        # X-axis
        majorLocator   = MultipleLocator(100)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(50)

        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        
        # Y-axis
        majorLocator   = MultipleLocator(0.1)
        majorFormatter = FormatStrFormatter(r'$\rm %s$ ')
        minorLocator   = MultipleLocator(0.05)
        
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)

        tight_layout()
        legend(scatterpoints=1,prop={'size':12},loc=2)
        xlim(0,max(lyaWArray)+1)
    
        if save:
            savefig('{0}/CDF_alt(lyaW_blue_vs_red).pdf'.format(saveDirectory),format='pdf')
        else:
            show()
    
###############################################################################
###############################################################################
###############################################################################
###############################################################################


if __name__=="__main__":
    # do the work
    main()
    