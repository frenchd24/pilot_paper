#!/usr/bin/env python

'''
By David French (frenchd@astro.wisc.edu)

$Id:  plotW_impact.py, v 5.6 9/26/16

Plot EW as a function of impact parameter, and impact parameter/diameter and /R_vir
    (01/04/2016)


This is the plotW_b_diam bit from histograms3.py. Now is separated, and loads in a pickle
file of the relevant data, as created by "buildDataLists.py"

Previous (from histograms3.py):
    Plot some stuff for the 100largest initial results

    Make plots for AAS winter 2014 poster. Uses LG_correlation_combined2.csv file

    Updated for the pilot paper (05/06/15)


v5.1: updated for LG_correlation_combined5_8_edit2.csv for l_min = 0.001 (02/24/2016)

v5.2: remake plots with v_hel instead of vcorr (4/21/16)

v5.3: remake plots with new large galaxy sample (7/13/16) -> /plots4/

v5.4: add the ability to limit results based on 'environment' number (7/14/16)
        also add a likelihood limit 


v5.5: major edits to structure and functions included. Same ideas, but better formatting
    and removed some duplicate functions. Made plots4/ for new pilot paper (8/05/16)
    
v5.6: update with LG_correlation_combined5_11_25cut_edit4.csv and /plots5/
    (9/26/16)

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
from scipy import stats


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

fontScale = 18
rc('text', usetex=True)
rc('font', size=18, family='serif', weight='normal')
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

def perc90(a):
    if len(a)>0:
        return percentile(a,90)
    else:
        return 0
        
def perc10(a):
    if len(a)>0:
        return percentile(a,10)
    else:
        return 0
        
def perc70(a):
    if len(a)>0:
        return percentile(a,70)
    else:
        return 0

    
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
    
    maxEnv = 3000
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
        
            if isNumber(pa):
                pa = float(pa)
            elif isNumber(RC3pa):
                pa = float(RC3pa)
            else:
                pa = -99
                
            if isNumber(az):
                az = float(az)
            else:
                az = -99
                
            if isNumber(maj):
                maj = float(maj)
                virialRadius = float(virialRadius)
            else:
                maj = -99
                virialRadius = -99
            
            # all the lists to be used for associated lines
            if float(env) <= maxEnv and float(likelihood) >= minL:
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
    
    
########################################################################################
########################################################################################

    # plot equivalent width as a function of impact parameter/diameter, split between
    # red and blue shifted absorption
    #
    
    plotW_b_diam= False
    save = False
    
    if plotW_b_diam:
        fig = figure()
        ax = fig.add_subplot(111)
        countb = 0
        countr = 0
        count = -1
        labelr = 'Redshifted Absorber'
        labelb = "Blueshifted Absorber"
        alpha = 0.7
        
        for d,i,w,m in zip(difList,impactList,lyaWList,majList):
            # check if all the values are good
            if isNumber(d) and isNumber(i) and isNumber(w) and isNumber(m):
                if d !=-99 and i !=-99 and w!=-99 and m!=-99:
                    count +=1
                    if d>0:
                        # galaxy is behind absorber, so gas is blue shifted
                        color = 'Blue'
                        if countb == 0:
                            countb +=1
                            plotb = ax.scatter(i/m,w,c='Blue',s=50,label= labelb,alpha=alpha)
                    if d<0:
                        # gas is red shifted compared to galaxy
                        color = 'Red'
                        if countr == 0:
                            countr +=1
                            plotr = ax.scatter(i/m,w,c='Red',s=50,label= labelr,alpha=alpha)
                    
                    plot1 = scatter(i/m,w,c=color,s = 50,alpha=alpha)
            
        # make the legend work properly
#         labelr = 'Red Shifted Absorber'
#         labelb = "Blue Shifted Absorber"
#         plotb = scatter(i[countb]/m[countb],w[countb],c='Blue',s=50,label= labelb)
#         plotr = scatter(i[countr]/m[countr],w[countr],c='Red',s=50,label= labelr)
        
#         title('W(impact/diameter) for red and blue shifted absorption')
        xlabel(r'$\rm \rho / D$')
        ylabel(r'$\rm Equivalent ~Width [m\AA]$')
        ax.grid(b=None,which='major',axis='both')
        ylim(-1,1200)
#         xlim(-1,150)

        ax.legend(scatterpoints=1)
        
        if save:
            savefig('{0}/W(impact_diam)_dif_cut.pdf'.format(saveDirectory),format='pdf')
        else:
            show()
    
    

##########################################################################################
##########################################################################################
    # plot equivalent width as a function of impact parameter, split between
    # red and blue shifted absorption, overplot median histograms for red and blue
    #
    
    plotW_impact_difhist = True
    save = True
    
    if plotW_impact_difhist:
        fig = figure(figsize=(7.7,5.7))
        ax = fig.add_subplot(111)
        
        countb = 0
        countr = 0
        count = -1
        alpha = 0.7
        
        binSize = 125
        bins = arange(0,625,binSize)
        
        labelr = r'$\rm Redshifted ~Absorber$'
        labelb = r'$\rm Blueshifted ~Absorber$'
        bSymbol = 'D'
        rSymbol = 'o'
        
        
        xVals = []
        redX = []
        redY = []
        blueX = []
        blueY = []
        
        for d,i,w,v in zip(difList,impactList,lyaWList,virList):
            # check if all the values are good
            if isNumber(d) and isNumber(i) and isNumber(w) and isNumber(v):
                if d!=-99 and i!=-99 and w!=-99 and v!=-99:
                    xVal = float(i)
                    yVal = float(w)
                    
                    xVals.append(xVal)
                    
                    if d>0:
                        # galaxy is behind absorber, so gas is blue shifted
                        color = 'Blue'
                        symbol = bSymbol
                        
                        blueX.append(xVal)
                        blueY.append(yVal)
                        
                        if countb == 0:
                            countb +=1
#                             plotb = ax.scatter(v,w,c='Blue',s=50,label= labelb)
                            plotb = ax.scatter(xVal,yVal,marker=symbol,c='Blue',s=50,\
                            alpha=alpha,label=labelb)

                    if d<0:
                        # gas is red shifted compared to galaxy
                        color = 'Red'
                        symbol = rSymbol
                        
                        redX.append(xVal)
                        redY.append(yVal)
                        
                        if countr == 0:
                            countr +=1
#                             plotr = ax.scatter(v,w,c='Red',s=50,label= labelr)
                            plotr = ax.scatter(xVal,yVal,marker=symbol,c='Red',s=50,\
                            alpha=alpha,label=labelr)

                    plot1 = scatter(xVal,yVal,marker=symbol,c=color,s=50,alpha=alpha)
    
        
#         # median red
#         bin_means,edges,binNumber = stats.binned_statistic(array(redX), array(redY), \
#         statistic='median', bins=bins)
#         left,right = edges[:-1],edges[1:]        
#         X = array([left,right]).T.flatten()
#         Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
#         plot(X,Y, ls='dotted',color='red',lw=1.7,alpha=alpha+0.1,label=r'$\rm Median~ Redshifted ~EW$')
#     
#     
#         # median blue
#         bin_means,edges,binNumber = stats.binned_statistic(array(blueX), array(blueY), \
#         statistic='median', bins=bins)
#         left,right = edges[:-1],edges[1:]        
#         X = array([left,right]).T.flatten()
#         Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
#         plot(X,Y, ls='dashed',color='blue',lw=1.7,alpha=alpha+0.1,label=r'$\rm Median~ Blueshifted ~EW$')        
        
        # avg red
        bin_means,edges,binNumber = stats.binned_statistic(array(redX), array(redY), \
        statistic='mean', bins=bins)
        left,right = edges[:-1],edges[1:]        
        X = array([left,right]).T.flatten()
        Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
        plot(X,Y, ls='dotted',color='red',lw=2.1,alpha=alpha+0.2,label=r'$\rm Mean~ Redshifted ~EW$')
    
    
        # avg blue
        bin_means,edges,binNumber = stats.binned_statistic(array(blueX), array(blueY), \
        statistic='mean', bins=bins)
        left,right = edges[:-1],edges[1:]        
        X = array([left,right]).T.flatten()
        Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
        plot(X,Y, ls='dashed',color='blue',lw=1.7,alpha=alpha+0.1,label=r'$\rm Mean~ Blueshifted ~EW$')    
        
        
        # x-axis
        majorLocator   = MultipleLocator(100)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(50)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        
        # y-axis
        majorLocator   = MultipleLocator(200)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(100)
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)
        
        xlabel(r'$\rm \rho ~[kpc]$')
        ylabel(r'$\rm Equivalent ~ Width ~ [m\AA]$')
        leg = ax.legend(scatterpoints=1,prop={'size':13},loc=1,fancybox=True)
#         leg.get_frame().set_alpha(0.5)
        ax.grid(b=None,which='major',axis='both')
        ylim(0,1300)
        xlim(0,500)

        if save:
            savefig('{0}/W(impact)_mean_{1}_difHistograms2.pdf'.format(saveDirectory,binSize),format='pdf',bbox_inches='tight')
        else:
            show()


##########################################################################################
##########################################################################################
    # plot equivalent width as a function of impact parameter/R_vir, split between
    # red and blue shifted absorption, overplot median histograms for red and blue
    #
    
    plotW_impact_vir_difhist = True
    save = True
    
    if plotW_impact_vir_difhist:
        fig = figure(figsize=(7.7,5.7))
        ax = fig.add_subplot(111)
        
        countb = 0
        countr = 0
        count = -1
        alpha = 0.7
        
        binSize = 0.5
        bins = arange(0,2.5,binSize)
        
        labelr = r'$\rm Redshifted ~Absorber$'
        labelb = r'$\rm Blueshifted ~Absorber$'
        bSymbol = 'D'
        rSymbol = 'o'
        
        
        xVals = []
        redX = []
        redY = []
        blueX = []
        blueY = []
        
        for d,i,w,v in zip(difList,impactList,lyaWList,virList):
            # check if all the values are good
            if isNumber(d) and isNumber(i) and isNumber(w) and isNumber(v):
                if d!=-99 and i!=-99 and w!=-99 and v!=-99:
                    xVal = float(i)/float(v)
                    yVal = float(w)
                    
                    xVals.append(xVal)
                    
                    if d>0:
                        # galaxy is behind absorber, so gas is blue shifted
                        color = 'Blue'
                        symbol = bSymbol
                        
                        blueX.append(xVal)
                        blueY.append(yVal)
                        
                        if countb == 0:
                            countb +=1
#                             plotb = ax.scatter(v,w,c='Blue',s=50,label= labelb)
                            plotb = ax.scatter(xVal,yVal,marker=symbol,c='Blue',s=50,\
                            alpha=alpha,label=labelb)

                    if d<0:
                        # gas is red shifted compared to galaxy
                        color = 'Red'
                        symbol = rSymbol
                        
                        redX.append(xVal)
                        redY.append(yVal)
                        
                        if countr == 0:
                            countr +=1
#                             plotr = ax.scatter(v,w,c='Red',s=50,label= labelr)
                            plotr = ax.scatter(xVal,yVal,marker=symbol,c='Red',s=50,\
                            alpha=alpha,label=labelr)

                    plot1 = scatter(xVal,yVal,marker=symbol,c=color,s=50,alpha=alpha)
        
#         # median red
#         bin_means,edges,binNumber = stats.binned_statistic(array(redX), array(redY), \
#         statistic='median', bins=bins)
#         left,right = edges[:-1],edges[1:]        
#         X = array([left,right]).T.flatten()
#         Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
#         plot(X,Y, ls='dotted',color='red',lw=1.7,alpha=alpha+0.1,label=r'$\rm Median~ Redshifted ~EW$')
#     
#     
#         # median blue
#         bin_means,edges,binNumber = stats.binned_statistic(array(blueX), array(blueY), \
#         statistic='median', bins=bins)
#         left,right = edges[:-1],edges[1:]        
#         X = array([left,right]).T.flatten()
#         Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
#         plot(X,Y, ls='dashed',color='blue',lw=1.7,alpha=alpha+0.1,label=r'$\rm Median~ Blueshifted ~EW$')        
        
        # avg red
        bin_means,edges,binNumber = stats.binned_statistic(array(redX), array(redY), \
        statistic='mean', bins=bins)
        left,right = edges[:-1],edges[1:]        
        X = array([left,right]).T.flatten()
        Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
        plot(X,Y, ls='dotted',color='red',lw=2.1,alpha=alpha+0.2,label=r'$\rm Mean~ Redshifted ~EW$')
    
    
        # avg blue
        bin_means,edges,binNumber = stats.binned_statistic(array(blueX), array(blueY), \
        statistic='mean', bins=bins)
        left,right = edges[:-1],edges[1:]        
        X = array([left,right]).T.flatten()
        Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
        plot(X,Y, ls='dashed',color='blue',lw=1.7,alpha=alpha+0.1,label=r'$\rm Mean~ Blueshifted ~EW$')    
        
        
        # x-axis
        majorLocator   = MultipleLocator(0.5)
        majorFormatter = FormatStrFormatter(r'$\rm %0.1f$')
        minorLocator   = MultipleLocator(0.25)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        
        # y-axis
        majorLocator   = MultipleLocator(200)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(100)
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)
        
        xlabel(r'$\rm \rho / R_{vir}$')
        ylabel(r'$\rm Equivalent ~ Width ~ [m\AA]$')
        leg = ax.legend(scatterpoints=1,prop={'size':13},loc=1,fancybox=True)
#         leg.get_frame().set_alpha(0.5)

        ax.grid(b=None,which='major',axis='both')
        ylim(0,1300)
        xlim(0,2.0)

        if save:
            savefig('{0}/W(impact_vir)_mean_{1}_difHistograms2.pdf'.format(saveDirectory,binSize),format='pdf',bbox_inches='tight')
        else:
            show()


##########################################################################################
##########################################################################################
    # plot equivalent width as a function of impact parameter/R_vir, split between
    # red and blue shifted absorption, overplot single median histogram (total EW)
    #
    
    plotW_impact_vir_medHist = False
    save = False
    
    if plotW_impact_vir_medHist:
        fig = figure()
        ax = fig.add_subplot(111)
        
        countb = 0
        countr = 0
        count = -1
        alpha = 0.7
        
        labelr = 'Redshifted Absorber'
        labelb = "Blueshifted Absorber"
        bSymbol = 'D'
        rSymbol = 'o'
        
        xVals = []
        yVals = []
        
        for d,i,w,v in zip(difList,impactList,lyaWList,virList):
            # check if all the values are good
            if isNumber(d) and isNumber(i) and isNumber(w) and isNumber(v):
                if d!=-99 and i!=-99 and w!=-99 and v!=-99:
                    xVal = float(i)/float(v)
                    yVal = float(w)
                    
                    xVals.append(xVal)
                    yVals.append(yVal)
                    
                    if d>0:
                        # galaxy is behind absorber, so gas is blue shifted
                        color = 'Blue'
                        symbol = bSymbol

                        
                        if countb == 0:
                            countb +=1
#                             plotb = ax.scatter(v,w,c='Blue',s=50,label= labelb)
                            plotb = ax.scatter(xVal,yVal,marker=bSymbol,c='Blue',s=50,alpha=alpha)

                    if d<0:
                        # gas is red shifted compared to galaxy
                        color = 'Red'
                        symbol = rSymbol
                        
                        if countr == 0:
                            countr +=1
#                             plotr = ax.scatter(v,w,c='Red',s=50,label= labelr)
                            plotr = ax.scatter(xVal,yVal,marker=rSymbol,c='Red',s=50,alpha=alpha)

                    plot1 = scatter(xVal,yVal,marker=symbol,c=color,s=50,alpha=alpha)
        
        
        # x-axis
        majorLocator   = MultipleLocator(0.5)
        majorFormatter = FormatStrFormatter(r'$\rm %0.1f$')
        minorLocator   = MultipleLocator(0.25)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        
        # y-axis
        majorLocator   = MultipleLocator(200)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(100)
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)
    
        binSize = 0.5
        bins = arange(0,2.5,binSize)

        # 50% percentile
        bin_means,edges,binNumber = stats.binned_statistic(array(xVals), array(yVals), \
        statistic='mean', bins=bins)
        left,right = edges[:-1],edges[1:]        
        X = array([left,right]).T.flatten()
        Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
        plot(X,Y, ls='solid',color='black',lw=1.7,alpha=alpha+0.1,label=r'$\rm Mean ~EW$')
        
#         # 90% percentile
#         bin_means,edges,binNumber = stats.binned_statistic(array(xVals), array(yVals), \
#         statistic=lambda y: perc90(y), bins=bins)
#         left,right = edges[:-1],edges[1:]
#         X = array([left,right]).T.flatten()
#         Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
#         plot(X,Y, ls='dashed',color='dimgray',lw=1.7,alpha=alpha+0.1,label=r'$\rm 90th\% ~EW$')
        
#         # 10th percentile
#         bin_means,edges,binNumber = stats.binned_statistic(array(xVals), array(yVals), \
#         statistic=lambda y: perc10(y), bins=bins)
#         left,right = edges[:-1],edges[1:]        
#         X = array([left,right]).T.flatten()
#         Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
#         plot(X,Y, ls='dotted',color='green',lw=1.7,alpha=alpha+0.1,label=r'$\rm 10th\% ~EW$')
        
        
        xlabel(r'$\rm \rho / R_{vir}$')
        ylabel(r'$\rm Equivalent ~ Width ~ [m\AA]$')
        ax.legend(scatterpoints=1,prop={'size':15},loc=1,fancybox=True)
        ax.grid(b=None,which='major',axis='both')
        ylim(0,1200)
        xlim(0,2.0)

        if save:
            savefig('{0}/W(impact_vir)_mean_{1}_Histograms.pdf'.format(saveDirectory,binSize),format='pdf',bbox_inches='tight')
        else:
            show()
            

##########################################################################################
##########################################################################################
    # plot equivalent width as a function of impact parameter, split between
    # red and blue shifted absorption, overplot single median histogram (total EW)
    #
    
    plotW_impact_medHist = False
    save = False
    
    if plotW_impact_medHist:
        fig = figure()
        ax = fig.add_subplot(111)
        
        countb = 0
        countr = 0
        count = -1
        alpha = 0.7
        
        labelr = 'Redshifted Absorber'
        labelb = "Blueshifted Absorber"
        bSymbol = 'D'
        rSymbol = 'o'
        
        xVals = []
        yVals = []
        
        for d,i,w,v in zip(difList,impactList,lyaWList,virList):
            # check if all the values are good
            if isNumber(d) and isNumber(i) and isNumber(w) and isNumber(v):
                if d!=-99 and i!=-99 and w!=-99 and v!=-99:
                    xVal = float(i)
                    yVal = float(w)
                    
                    xVals.append(xVal)
                    yVals.append(yVal)
                    
                    if d>0:
                        # galaxy is behind absorber, so gas is blue shifted
                        color = 'Blue'
                        symbol = bSymbol
                        
                        if countb == 0:
                            countb +=1
#                             plotb = ax.scatter(v,w,c='Blue',s=50,label= labelb)
                            plotb = ax.scatter(xVal,yVal,marker=bSymbol,c='Blue',s=50,alpha=alpha)

                    if d<0:
                        # gas is red shifted compared to galaxy
                        color = 'Red'
                        symbol = rSymbol
                        
                        if countr == 0:
                            countr +=1
#                             plotr = ax.scatter(v,w,c='Red',s=50,label= labelr)
                            plotr = ax.scatter(xVal,yVal,marker=rSymbol,c='Red',s=50,alpha=alpha)

                    plot1 = scatter(xVal,yVal,marker=symbol,c=color,s=50,alpha=alpha)
        
        
        # x-axis
        majorLocator   = MultipleLocator(100)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(50)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        
        # y-axis
        majorLocator   = MultipleLocator(200)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(100)
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)
    
        binSize = 100
        bins = arange(0,600,binSize)
        
        # 50% percentile
        bin_means,edges,binNumber = stats.binned_statistic(array(xVals), array(yVals), \
        statistic='mean', bins=bins)
        left,right = edges[:-1],edges[1:]        
        X = array([left,right]).T.flatten()
        Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
        plot(X,Y, ls='solid',color='black',lw=1.5,alpha=alpha,label=r'$\rm Mean ~EW$')
        
#         # 90% percentile
#         bin_means,edges,binNumber = stats.binned_statistic(array(xVals), array(yVals), \
#         statistic=lambda y: perc90(y), bins=bins)
#         left,right = edges[:-1],edges[1:]        
#         X = array([left,right]).T.flatten()
#         Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
#         plot(X,Y, ls='dashed',color='dimgray',lw=1.5,alpha=alpha,label=r'$\rm 90th\% ~EW$')
        
#         # 10th percentile
#         bin_means,edges,binNumber = stats.binned_statistic(array(xVals), array(yVals), \
#         statistic=lambda y: perc10(y), bins=bins)
#         left,right = edges[:-1],edges[1:]        
#         X = array([left,right]).T.flatten()
#         Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
#         plot(X,Y, ls='dotted',color='green',lw=1.5,alpha=alpha,label=r'$\rm 10th\% ~EW$')
        
        xlabel(r'$\rm \rho ~[kpc]$')
        ylabel(r'$\rm Equivalent ~ Width ~ [m\AA]$')
        ax.legend(scatterpoints=1,prop={'size':15},loc=1,fancybox=True)
        ax.grid(b=None,which='major',axis='both')
        ylim(0,1200)
        xlim(0,500)

        if save:
            savefig('{0}/W(impact)_mean_{1}_Histograms.pdf'.format(saveDirectory,binSize),format='pdf',bbox_inches='tight')
        else:
            show()



##########################################################################################
##########################################################################################
    # plot equivalent width as a function of impact parameter, split between
    # absorbers at azimuth >45 and <45, overplot median histograms for each
    #
    
    plotW_impact_az = False
    save = False
    
    if plotW_impact_az:
        fig = figure()
        ax = fig.add_subplot(111)
        
        countb = 0
        countr = 0
        count = -1
        alpha = 0.7
        binSize = 100
        bins = arange(0,600,binSize)
        
        labelr = 'Az < 45 Absorber'
        labelb = "Az > 45 Absorber"
        bSymbol = 'D'
        rSymbol = 'o'
        
        xVals = []
        redX = []
        redY = []
        blueX = []
        blueY = []
        
        for d,i,w,v,a in zip(difList,impactList,lyaWList,virList,azList):
            # check if all the values are good
            if isNumber(d) and isNumber(i) and isNumber(w) and isNumber(v) and isNumber(a):
                if d!=-99 and i!=-99 and w!=-99 and v!=-99 and a!=-99:
#                     xVal = float(i)/float(v)
                    xVal = float(i)
                    yVal = float(w)
                    
                    xVals.append(xVal)
                    
                    if float(a)>=45:
                        # galaxy is behind absorber, so gas is blue shifted
                        color = 'Blue'
                        symbol = bSymbol
                        
                        blueX.append(xVal)
                        blueY.append(yVal)
                        
                        if countb == 0:
                            countb +=1
                            plotb = ax.scatter(xVal,yVal,marker=symbol,c='Blue',s=50,alpha=alpha)

                    if float(a)<45:
                        # gas is red shifted compared to galaxy
                        color = 'Red'
                        symbol = rSymbol
                        
                        redX.append(xVal)
                        redY.append(yVal)
                        
                        if countr == 0:
                            countr +=1
                            plotr = ax.scatter(xVal,yVal,marker=symbol,c='Red',s=50,alpha=alpha)

                    plot1 = scatter(xVal,yVal,marker=symbol,c=color,s=50,alpha=alpha)
     
        
        # avg AZ < 45 = RED
        bin_means,edges,binNumber = stats.binned_statistic(array(redX), array(redY), \
        statistic='mean', bins=bins)
        left,right = edges[:-1],edges[1:]        
        X = array([left,right]).T.flatten()
        Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
        plot(X,Y, ls='dotted',color='red',lw=1.7,alpha=alpha+0.1,label=r'$\rm Az < 45$')
    
    
        # avg AZ >= 45 = BLUE
        bin_means,edges,binNumber = stats.binned_statistic(array(blueX), array(blueY), \
        statistic='mean', bins=bins)
        left,right = edges[:-1],edges[1:]        
        X = array([left,right]).T.flatten()
        Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
        plot(X,Y, ls='dashed',color='blue',lw=1.7,alpha=alpha+0.1,label=r'$\rm Az \geq 45$')
        

        # x-axis
        majorLocator   = MultipleLocator(100)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(50)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        
        # y-axis
        majorLocator   = MultipleLocator(200)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(100)
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)
        
        xlabel(r'$\rm Impact ~ Parameter ~ [kpc]$')
        ylabel(r'$\rm Equivalent ~ Width ~ [m\AA]$')
        ax.legend(scatterpoints=1,prop={'size':15},loc=1,fancybox=True)
        ax.grid(b=None,which='major',axis='both')
        ylim(0,1000)
#         xlim(0,2.0)
        xlim(0,500)


        if save:
            savefig('{0}/W(impact)_mean_{1}_az.pdf'.format(saveDirectory,binSize),format='pdf',bbox_inches='tight')
        else:
            show()


##########################################################################################
##########################################################################################
    # plot equivalent width as a function of impact parameter/R_vir, split between
    # absorbers at inc >55 and <55, overplot median histograms for each
    #
    
    plotW_impact_vir_inc = False
    save = False
    
    if plotW_impact_vir_inc:
        fig = figure()
        ax = fig.add_subplot(111)
        
        countb = 0
        countr = 0
        count = -1
        alpha = 0.7
        binSize = 0.5
        bins = arange(0,2.5,binSize)
        
        labelr = 'Inc < 55 Absorber'
        labelb = "Inc > 55 Absorber"
        bSymbol = 'D'
        rSymbol = 'o'
        
        xVals = []
        redX = []
        redY = []
        blueX = []
        blueY = []
        
        for d,i,w,v,inc in zip(difList,impactList,lyaWList,virList,incList):
            # check if all the values are good
            if isNumber(d) and isNumber(i) and isNumber(w) and isNumber(v) and isNumber(inc):
                if d!=-99 and i!=-99 and w!=-99 and v!=-99 and inc!=-99:
                    xVal = float(i)/float(v)
                    yVal = float(w)
                    
                    xVals.append(xVal)
                    
                    if float(inc)>55:
                        # galaxy inc > 55
                        color = 'Blue'
                        symbol = bSymbol
                        
                        blueX.append(xVal)
                        blueY.append(yVal)
                        
                        if countb == 0:
                            countb +=1
                            plotb = ax.scatter(xVal,yVal,marker=symbol,c='Blue',s=50,alpha=alpha)

                    if float(inc)<55:
                        # galaxy in < 55
                        color = 'Red'
                        symbol = rSymbol
                        
                        redX.append(xVal)
                        redY.append(yVal)
                        
                        if countr == 0:
                            countr +=1
                            plotr = ax.scatter(xVal,yVal,marker=symbol,c='Red',s=50,alpha=alpha)

                    plot1 = scatter(xVal,yVal,marker=symbol,c=color,s=50,alpha=alpha)
     
        
        # avg Inc < 55 = RED
        bin_means,edges,binNumber = stats.binned_statistic(array(redX), array(redY), \
        statistic='mean', bins=bins)
        left,right = edges[:-1],edges[1:]        
        X = array([left,right]).T.flatten()
        Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
        plot(X,Y, ls='dotted',color='red',lw=1.7,alpha=alpha+0.1,label=r'$\rm Inc < 55$')
    
    
        # avg Inc > 55 = BLUE
        bin_means,edges,binNumber = stats.binned_statistic(array(blueX), array(blueY), \
        statistic='mean', bins=bins)
        left,right = edges[:-1],edges[1:]        
        X = array([left,right]).T.flatten()
        Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
        plot(X,Y, ls='dashed',color='blue',lw=1.7,alpha=alpha+0.1,label=r'$\rm Inc > 55$')
        
        
        # x-axis
        majorLocator   = MultipleLocator(0.5)
        majorFormatter = FormatStrFormatter(r'$\rm %0.1f$')
        minorLocator   = MultipleLocator(0.25)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        
        # y-axis
        majorLocator   = MultipleLocator(200)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(100)
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)
        
        xlabel(r'$\rm \rho / R_{vir}$')
        ylabel(r'$\rm Equivalent ~ Width ~ [m\AA]$')
        ax.legend(scatterpoints=1,prop={'size':15},loc=1,fancybox=True)
        ax.grid(b=None,which='major',axis='both')
        ylim(0,1200)
        xlim(0,2.0)

        if save:
            savefig('{0}/W(impact_vir)_mean_{1}_inc.pdf'.format(saveDirectory,binSize),format='pdf',bbox_inches='tight')
        else:
            show()


##########################################################################################
##########################################################################################
    # plot equivalent width as a function of impact parameter, split between
    # absorbers at inc >55 and <55, overplot median histograms for each
    #
    
    plotW_impact_inc = False
    save = False
    
    if plotW_impact_inc:
        fig = figure()
        ax = fig.add_subplot(111)
        
        countb = 0
        countr = 0
        count = -1
        alpha = 0.7
        binSize = 125
        bins = arange(0,625,binSize)
        
        labelr = 'Inc < 55 Absorber'
        labelb = "Inc > 55 Absorber"
        bSymbol = 'D'
        rSymbol = 'o'
        
        xVals = []
        redX = []
        redY = []
        blueX = []
        blueY = []
        
        for d,i,w,v,inc in zip(difList,impactList,lyaWList,virList,incList):
            # check if all the values are good
            if isNumber(d) and isNumber(i) and isNumber(w) and isNumber(v) and isNumber(inc):
                if d!=-99 and i!=-99 and w!=-99 and v!=-99 and inc!=-99:
                    xVal = float(i)
                    yVal = float(w)
                    
                    xVals.append(xVal)
                    
                    if float(inc)>55:
                        # galaxy inc > 55
                        color = 'Blue'
                        symbol = bSymbol
                        
                        blueX.append(xVal)
                        blueY.append(yVal)
                        
                        if countb == 0:
                            countb +=1
                            plotb = ax.scatter(xVal,yVal,marker=symbol,c='Blue',s=50,alpha=alpha)

                    if float(inc)<55:
                        # galaxy inc < 55
                        color = 'Red'
                        symbol = rSymbol
                        
                        redX.append(xVal)
                        redY.append(yVal)
                        
                        if countr == 0:
                            countr +=1
                            plotr = ax.scatter(xVal,yVal,marker=symbol,c='Red',s=50,alpha=alpha)

                    plot1 = scatter(xVal,yVal,marker=symbol,c=color,s=50,alpha=alpha)
     
        
        # avg inc < 55 = RED
        bin_means,edges,binNumber = stats.binned_statistic(array(redX), array(redY), \
        statistic='mean', bins=bins)
        left,right = edges[:-1],edges[1:]        
        X = array([left,right]).T.flatten()
        Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
        plot(X,Y, ls='dotted',color='red',lw=1.7,alpha=alpha+0.1,label=r'$\rm Inc < 55$')
    
    
        # avg inc > 55 = BLUE
        bin_means,edges,binNumber = stats.binned_statistic(array(blueX), array(blueY), \
        statistic='mean', bins=bins)
        left,right = edges[:-1],edges[1:]        
        X = array([left,right]).T.flatten()
        Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
        plot(X,Y, ls='dashed',color='blue',lw=1.7,alpha=alpha+0.1,label=r'$\rm Inc > 55$')
        
        
        # x-axis
        majorLocator   = MultipleLocator(100)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(50)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        
        # y-axis
        majorLocator   = MultipleLocator(200)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(100)
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)
        
        xlabel(r'$\rm \rho ~[kpc]$')
        ylabel(r'$\rm Equivalent ~ Width ~ [m\AA]$')
        ax.legend(scatterpoints=1,prop={'size':15},loc=1,fancybox=True)
        ax.grid(b=None,which='major',axis='both')
        ylim(0,1200)
        xlim(0,500)

        if save:
            savefig('{0}/W(impact)_mean_{1}_inc.pdf'.format(saveDirectory,binSize),format='pdf',bbox_inches='tight')
        else:
            show()


##########################################################################################
##########################################################################################
    # plot equivalent width as a function of impact parameter, split between
    # absorbers at R_vir <200, R_vir >=200, overplot median histograms for each
    #
    
    plotW_impact_vir_separate = False
    save = False
    
    if plotW_impact_vir_separate:
        fig = figure()
        ax = fig.add_subplot(111)
        
        countb = 0
        countr = 0
        count = -1
        alpha = 0.7
        binSize = 125
        bins = arange(0,625,binSize)
        
        labelr = r'$\rm R_{vir} ~ < ~ 200~kpc$'
        labelb = r'$\rm R_{vir} ~ >= ~ 200~kpc$'
        bSymbol = 'D'
        rSymbol = 'o'
        
        xVals = []
        redX = []
        redY = []
        blueX = []
        blueY = []
        
        for d,i,w,v,inc in zip(difList,impactList,lyaWList,virList,incList):
            # check if all the values are good
            if isNumber(d) and isNumber(i) and isNumber(w) and isNumber(v) and isNumber(inc):
                if d!=-99 and i!=-99 and w!=-99 and v!=-99 and inc!=-99:
                    xVal = float(i)
                    yVal = float(w)
                    
                    xVals.append(xVal)
                    
                    if float(v)>200:
                        # galaxy inc > 55
                        color = 'Blue'
                        symbol = bSymbol
                        
                        blueX.append(xVal)
                        blueY.append(yVal)
                        
                        if countb == 0:
                            countb +=1
                            plotb = ax.scatter(xVal,yVal,marker=symbol,c='Blue',s=50,alpha=alpha)

                    if float(v)<=200:
                        # galaxy inc < 55
                        color = 'Red'
                        symbol = rSymbol
                        
                        redX.append(xVal)
                        redY.append(yVal)
                        
                        if countr == 0:
                            countr +=1
                            plotr = ax.scatter(xVal,yVal,marker=symbol,c='Red',s=50,alpha=alpha)

                    plot1 = scatter(xVal,yVal,marker=symbol,c=color,s=50,alpha=alpha)
     
        
        # avg inc < 55 = RED
        bin_means,edges,binNumber = stats.binned_statistic(array(redX), array(redY), \
        statistic='mean', bins=bins)
        left,right = edges[:-1],edges[1:]        
        X = array([left,right]).T.flatten()
        Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
        plot(X,Y, ls='dotted',color='red',lw=1.7,alpha=alpha+0.1,label=r'$\rm R_{vir} < 200~kpc$')
    
    
        # avg inc > 55 = BLUE
        bin_means,edges,binNumber = stats.binned_statistic(array(blueX), array(blueY), \
        statistic='mean', bins=bins)
        left,right = edges[:-1],edges[1:]        
        X = array([left,right]).T.flatten()
        Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
        plot(X,Y, ls='dashed',color='blue',lw=1.7,alpha=alpha+0.1,label=r'$\rm R_{vir} >= 200~kpc$')
        
        
        # x-axis
        majorLocator   = MultipleLocator(100)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(50)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        
        # y-axis
        majorLocator   = MultipleLocator(200)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(100)
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)
        
        xlabel(r'$\rm \rho ~[kpc]$')
        ylabel(r'$\rm Equivalent ~ Width ~ [m\AA]$')
        ax.legend(scatterpoints=1,prop={'size':15},loc=1,fancybox=True)
        ax.grid(b=None,which='major',axis='both')
        ylim(0,1200)
        xlim(0,500)

        if save:
            savefig('{0}/W(impact)_mean_{1}_vir_sep.pdf'.format(saveDirectory,binSize),format='pdf',bbox_inches='tight')
        else:
            show()



#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################


if __name__=="__main__":
    # do the work
    main()
    