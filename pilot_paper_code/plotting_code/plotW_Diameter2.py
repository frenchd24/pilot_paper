#!/usr/bin/env python

'''
By David French (frenchd@astro.wisc.edu)

$Id:  plotW_Diameter2.py, v 5.6 9/26/16

Plot equivalent width, NaV and Doppler parameter as a function of galaxy diameter and R_vir


This is the plotW_Diameter bit from histograms3.py. Now is separated, and loads in a pickle
file of the relevant data, as created by "buildDataLists.py"


Previous (from histograms3.py):
    Plot some stuff for the 100largest initial results

    Make plots for AAS winter 2014 poster. Uses LG_correlation_combined2.csv file

    Updated for the pilot paper (05/06/15)

v5: updated for pilot paper.
    (1/5/16)
    
v5.1: updated for LG_correlation_combined5_8_edit2.csv with l_min = 0.001 (02/24/2016)

v5.2: remake plots with v_hel instead of vcorr (4/22/16)

v5.3: remake plots with new large galaxy sample (7/13/16) -> /plots4/

v5.4: add the ability to limit results based on 'environment' number (7/14/16)

v5.5: minor formatting updates for pilot paper (8/08/16)

v5.6: update for LG_correlation_combined5_11_25cut_edit4.csv -> /plots5/ (9/26/16)
    
'''

import sys
import os
import csv
from scipy import stats

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
        resultsFilename = '/usr/users/frenchd/inclination/git_inclination/LG_correlation_combined5_11_25_cut_edit4.csv'
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
            major = l['majorAxis (kpc)']
            minor = l['minorAxis (kpc)']
            inc = l['inclination (deg)']
            az = l['azimuth (deg)']
            b = l['b'].partition('pm')[0]
            b_err = l['b'].partition('pm')[2]
            na = eval(l['Na'].partition(' pm ')[0])
#             print "l['Na'].partition(' pm ')[2] : ",l['Na'].partition(' pm ')
            na_err = eval(l['Na'].partition(' pm ')[2])
            likelihood = l['likelihood']
            likelihoodm15 = l['likelihood_1.5']
            virialRadius = l['virialRadius']
            m15 = l['d^1.5']
            vel_diff = l['vel_diff']
            
            if isNumber(inc):
                cosInc = cos(float(inc) * pi/180.)
                
                if isNumber(major) and isNumber(minor):
                    q0 = 0.2
                    fancyInc = calculateFancyInclination(major,minor,q0)
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
                
            if isNumber(major):
                major = float(major)
                virialRadius = float(virialRadius)
            else:
                major = -99
                virialRadius = -99
                
            if isNumber(b):
                b = float(b)
            else:
                b = -99
            
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
                fancyIncList.append(fancyInc)
                cosIncList.append(cosInc)
                cosFancyIncList.append(cosFancyInc)
                paList.append(pa)
                vcorrList.append(vcorr)
                majList.append(major)
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
    

##########################################################################################
##########################################################################################
    # plot equivalent width as a function of galaxy diameter
    #
    
    plotW_Diameter = False
    save = False
    
    if plotW_Diameter:
        fig = figure()
        ax = fig.add_subplot(111)
        countb = 0
        countr = 0
        count = -1
        labelr = 'Red Shifted Absorber'
        labelb = "Blue Shifted Absorber"
        for d,m,w in zip(difList,majList,lyaWList):
            if isNumber(d) and isNumber(m) and isNumber(w):
                if m != -99:
                    count +=1
                    if d>0:
                        # galaxy is behind absorber, so gas is blue shifted
                        color = 'Blue'
                        if countb == 0:
                            countb +=1
                            plotb = ax.scatter(m,w,c='Blue',s=50,label= labelb)
                    if d<0:
                        # gas is red shifted compared to galaxy
                        color = 'Red'
                        if countr == 0:
                            countr +=1
                            plotr = ax.scatter(m,w,c='Red',s=50,label= labelr)
                    
                    plot1 = scatter(m,w,c=color,s = 50)
        
#         title('Equivalent width vs. galaxy diameter')
        xlabel('Major Axis (kpc)')
        ylabel(r'Equivalent Width ($\rm m\AA$)')
        ax.grid(b=None,which='major',axis='both')
        ylim(0,1000)
        xlim(0,max(majList)+1)
        ax.legend(scatterpoints=1,prop={'size':15},loc=2)
        
        if save:
            savefig('{0}/W(diameter).pdf'.format(saveDirectory),format='pdf')
        else:
            show()
    
    
##########################################################################################
##########################################################################################
    # plot equivalent width as a function of galaxy R_vir
    #
    
    plotW_vir = False
    save = False
    
    if plotW_vir:
        fig = figure()
        ax = fig.add_subplot(111)
        countb = 0
        countr = 0
        count = -1
        labelr = 'Red Shifted Absorber'
        labelb = "Blue Shifted Absorber"
        for d,v,w in zip(difList,virList,lyaWList):
            if isNumber(d) and isNumber(v) and isNumber(w):
                if v != -99:
                    count +=1
                    if d>0:
                        # galaxy is behind absorber, so gas is blue shifted
                        color = 'Blue'
                        if countb == 0:
                            countb +=1
                            plotb = ax.scatter(v,w,c='Blue',s=50,label= labelb)
                    if d<0:
                        # gas is red shifted compared to galaxy
                        color = 'Red'
                        if countr == 0:
                            countr +=1
                            plotr = ax.scatter(v,w,c='Red',s=50,label= labelr)
                    
                    plot1 = scatter(v,w,c=color,s = 50)
        
#         title(r'Equivalent width vs. $\rm R_{vir}$')
        xlabel(r'$\rm R_{vir}$ (kpc)')
        ylabel(r'Equivalent Width ($\rm m\AA$)')
        ax.grid(b=None,which='major',axis='both')
        ylim(0,1200)
        xlim(0,max(virList)+5)
        ax.legend(scatterpoints=1,prop={'size':12},loc=2)
        
        if save:
            savefig('{0}/W(vir).pdf'.format(saveDirectory),format='pdf')
        else:
            show()
    
    
##########################################################################################
##########################################################################################
    # plot doppler parameter as a function of impact / R_vir
    #
    
    plotB_vir = False
    save = False
    
    if plotB_vir:
        fig = figure()
        ax = fig.add_subplot(111)
        countb = 0
        countr = 0
        count = -1

        bSymbol = 'D'
        rSymbol = 'o'
        alpha = 0.7
        
        labelr = 'Red Shifted Absorber'
        labelb = "Blue Shifted Absorber"
        for d,v,b,i in zip(difList,virList,bList,impactList):
            if isNumber(d) and isNumber(v) and isNumber(b) and isNumber(i):
                xVal = float(i)/float(v)
                yVal = float(b)
                if v != -99:
                    count +=1
                    if d>0:
                        # galaxy is behind absorber, so gas is blue shifted
                        color = 'Blue'
                        symbol = bSymbol
                        if countb == 0:
                            countb +=1
                            plotb = ax.scatter(xVal,yVal,c='Blue',s=50,label= labelb,marker=symbol,alpha=alpha)
                    if d<0:
                        # gas is red shifted compared to galaxy
                        color = 'Red'
                        symbol = rSymbol
                        if countr == 0:
                            countr +=1
                            plotr = ax.scatter(xVal,yVal,c='Red',s=50,label= labelr,marker=symbol,alpha=alpha)
                    
                    plot1 = scatter(xVal,yVal,c=color,s = 50,marker=symbol,alpha=alpha)
                    
        # x-axis
        majorLocator   = MultipleLocator(0.5)
        majorFormatter = FormatStrFormatter(r'$\rm %0.1f$')
        minorLocator   = MultipleLocator(0.25)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        
        # y-axis
        majorLocator   = MultipleLocator(10)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(5)
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)
        
        xlabel(r'$\rm \rho / R_{vir}$')
        ylabel(r'$\rm Doppler ~b ~Parameter ~[km/s]$')
        ax.grid(b=None,which='major',axis='both')
        ylim(min(bList)-5,max(bList)+5)
        xlim(0,2.0)
        ax.legend(scatterpoints=1,prop={'size':14},loc=2)
        
        if save:
            savefig('{0}/B(vir).pdf'.format(saveDirectory),format='pdf')
        else:
            show()



##########################################################################################
##########################################################################################
    # plot equivalent width as a function of virial radius for red vs blue
    # shifted absorption, include average histograms
    #
    
    plotW_vir_avg = True
    save = True
    
    if plotW_vir_avg:
        fig = figure(figsize=(7.7,5.7))
        ax = fig.add_subplot(111)
        
        countb = 0
        countr = 0
        count = -1
        alpha = 0.7
        binSize = 50

        
        bSymbol = 'D'
        rSymbol = 'o'
        labelr = 'Redshifted Absorber'
        labelb = "Blueshifted Absorber"
        
        placeArrayr = zeros(7)
        placeCountr = zeros(7)
        placeArrayb = zeros(7)
        placeCountb = zeros(7)
        redW = []
        blueW = []
        redVir = []
        blueVir = []
        
        
        for d,v,w,m in zip(difList,virList,lyaWList,majList):
            # check if all the values are good
            if isNumber(d) and isNumber(v) and isNumber(w) and isNumber(m):
                if d!=-99 and v!=-99 and w!=-99 and m!=-99:
                    
                    if d>0:
                        # galaxy is behind absorber, so gas is blue shifted
                        color = 'Blue'
                        symbol = bSymbol
                        
                        blueW.append(w)
                        blueVir.append(v)
                        
                        # which bin does it belong too?
                        place = v/binSize
                        print 'place: ',place
                        placeArrayb[place] += float(w)
                        print 'placeArrayb: ',placeArrayb
                        placeCountb[place] +=1.
                        print 'placecountb: ',placeCountb
                        
                        if countb == 0:
                            countb +=1
                            plotb = ax.scatter(v,w,marker=symbol,alpha=alpha,c='Blue',s=50)

                    if d<0:
                        # gas is red shifted compared to galaxy
                        color = 'Red'
                        symbol = rSymbol
                        
                        redW.append(w)
                        redVir.append(v)
                        
                        # which bin does it belong too?
                        place = v/binSize
                        placeArrayr[place] += float(w)
                        placeCountr[place] +=1.
                        
                        if countr == 0:
                            countr +=1
                            plotr = ax.scatter(v,w,marker=symbol,alpha=alpha,c='Red',s=50)
            
                    plot1 = scatter(v,w,marker=symbol,alpha=alpha,c=color,s=50)
                    
        rHist = placeArrayr/placeCountr
        print 'rHist: ',rHist
        bHist = placeArrayb/placeCountb
        print 'bHist: ',bHist
        
        totalrHist = []
        totalrVir = []
        totalbHist = []
        totalbVir = []
        
        for r,v in zip(rHist,arange(0,max(virList),binSize)):
            if not isNumber(r):
                r = 0
            
            totalrHist.append(r)
            totalrHist.append(r)

            totalrVir.append(v)
            totalrVir.append(v+binSize)
            
        for b,v in zip(bHist,arange(0,max(virList),binSize)):
            if not isNumber(b):
                b = 0
            totalbHist.append(b)
            totalbHist.append(b)

            totalbVir.append(v)
            totalbVir.append(v+binSize)
        
        print 'totalrVir: ',totalrVir
        print 'totalrHist: ',totalrHist
        print
        print 'totalbVir: ',totalbVir
        print 'totalbHist: ',totalbHist
        print
        
        
        # x-axis
        majorLocator   = MultipleLocator(50)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(25)
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
        
        # bins are in Rvir
        bins = arange(150,400,binSize)
           
#         bin_means,edges,binNumber = stats.binned_statistic(array(redVir), array(redW), statistic='mean', bins=bins)
#         left,right = edges[:-1],edges[1:]
#         X = array([left,right]).T.flatten()
#         Y = array([bin_means,bin_means]).T.flatten()
#         plot(X,Y, c='red',ls='dotted',lw=2.5,alpha=alpha,label=r'\rm $Mean ~Redshifted ~EW$')
#             
#         bin_means,edges,binNumber = stats.binned_statistic(array(blueVir), array(blueW), statistic='mean', bins=bins)
#         left,right = edges[:-1],edges[1:]
#         X = array([left,right]).T.flatten()
#         Y = array([bin_means,bin_means]).T.flatten()
#         plot(X,Y, c='blue',ls='dashed',lw=1.5,alpha=alpha,label=r'\rm $Mean ~Blueshifted ~EW$')
        
        plot2 = ax.plot(totalrVir,totalrHist,c='Red',lw=2.5,ls='dotted',label=r'$\rm Mean ~ Redshifted ~ EW$')
        plot3 = ax.plot(totalbVir,totalbHist,c='Blue',lw=1.5,ls='dashed',label=r'$\rm Mean ~ Blueshifted ~ EW$')
        
        xlabel(r'$\rm R_{vir} ~ [kpc]$')
        ylabel(r'$\rm Equivalent ~ Width ~ [m\AA]$')
        ax.legend(scatterpoints=1,prop={'size':15},loc=2,fancybox=True)
        ax.grid(b=None,which='major',axis='both')
        ylim(0,1200)
        xlim(150,350)

        if save:
            savefig('{0}/W(vir)_avgHistograms2.pdf'.format(saveDirectory),format='pdf',bbox_inches='tight')
        else:
            show()


##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################

if __name__=="__main__":
    # do the work
    main()
    