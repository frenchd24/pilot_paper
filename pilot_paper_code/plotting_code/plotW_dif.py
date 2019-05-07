#!/usr/bin/env python

'''
By David French (frenchd@astro.wisc.edu)

$Id:  plotW_dif.py, v 1.6 9/23/16

Plot EW as a function of velocity difference



Previous (from histograms3.py):
    Plot some stuff for the 100largest initial results

    Make plots for AAS winter 2014 poster. Uses LG_correlation_combined2.csv file

    Updated for the pilot paper (05/06/15)
    
    
v1.1: updated for LG_correlation_combined5_8_edit2.csv for l_min 0.001 (02/24/2016)

v1.2: remake plots with v_hel instead of vcorr (4/21/16)

    
v1.3: remake plots with v_hel instead of vcorr (4/21/16)

v1.4: remake plots with newest large galaxy sample (7/13/16) -> /plots4/

v1.5: minor formatting updates for pilot paper (8/08/16)

v1.6 minor update for LG_correlation_combined5_11_25cut_edit4.csv (9/23/16) -/plots5/

'''

import sys
import os
import csv

from pylab import *
from scipy import stats
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
    
    
##########################################################################################
##########################################################################################
    # plot equivalent width as a function of vel_diff for red and blue shifted
    # absorption
    #
    
    plotW_diff = True
    save = True
    
    if plotW_diff:
        fig = figure(figsize=(7.7,5.7))
        ax = fig.add_subplot(111)
        countb = 0
        countr = 0
        count = -1
        labelr = 'Redshifted Absorber'
        labelb = "Blueshifted Absorber"
        alpha = 0.7
        bSymbol = 'D'
        rSymbol = 'o'
        
        rdif = []
        bdif = []
        rW = []
        bW = []
        
        for d,w,m in zip(difList,lyaWList,majList):
            # check if all the values are okay
            if isNumber(d) and isNumber(w) and isNumber(m):
                if d!=-99 and w!=-99 and m!=-99:
                    if d>0:
                        # galaxy is behind absorber, so gas is blue shifted
                        color = 'Blue'
                        symbol = bSymbol
                        
                        bdif.append(float(d))
                        bW.append(float(w))
                        
                        if countb == 0:
                            countb +=1
                            plotb = ax.scatter(d,w,marker=symbol,c='Blue',s=50,label= labelb,\
                            alpha = alpha)
                    if d<0:
                        # gas is red shifted compared to galaxy
                        color = 'Red'
                        symbol = rSymbol
                        
                        rdif.append(float(d))
                        rW.append(float(w))
                        
                        if countr == 0:
                            countr +=1
                            plotr = ax.scatter(d,w,marker=symbol,c='Red',s=50,label= labelr,\
                            alpha = alpha)
                
                    plot1 = scatter(d,w,marker=symbol,c=color,s=50,alpha = alpha)
        
        includeHist = False

        if includeHist:
            bins = arange(-400,100,100)
            bin_means,edges,binNumber = stats.binned_statistic(array(rdif), array(rW), statistic='median', bins=bins)
            left,right = edges[:-1],edges[1:]
            X = np.array([left,right]).T.flatten()
            Y = np.array([bin_means,bin_means]).T.flatten()
            plt.plot(X,Y, c='Red',ls='dashed',lw=2,alpha=alpha,label=r'$\rm Mean EW$')
        
        
            bins = arange(0,500,100)
            bin_means,edges,binNumber = stats.binned_statistic(array(bdif), array(bW), statistic='median', bins=bins)
            left,right = edges[:-1],edges[1:]
            X = np.array([left,right]).T.flatten()
            Y = np.array([bin_means,bin_means]).T.flatten()
            plt.plot(X,Y, c='Blue',ls='dashed',lw=2,alpha=alpha,label=r'$\rm Mean EW$')

        # x-axis
        majorLocator   = MultipleLocator(100)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(50)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        
        # y axis
        majorLocator   = MultipleLocator(200)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(100)
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)


        xlabel(r'$\rm \Delta v ~ [km ~ s^{-1}]$')
        ylabel(r'$\rm Equivalent ~ Width ~ [m\AA]$')
#         legend(scatterpoints=1,prop={'size':12},loc=1)
        ax.grid(b=None,which='major',axis='both')
        ylim(0,1000)
        xlim(-400,400)
        
        if save:
            savefig('{0}/W(vel_diff)2.pdf'.format(saveDirectory),format='pdf',bbox_inches='tight')
        else:
            show()
    
    
##########################################################################################
##########################################################################################
    # plot equivalent width as a function of vel_diff * impact for red and blue shifted
    # absorption
    #
    
    plotW_impact_diff = False
    save = False
    
    if plotW_impact_diff:
        fig = figure()
        ax = fig.add_subplot(111)
        countb = 0
        countr = 0
        count = -1
        labelr = 'Redshifted Absorber'
        labelb = "Blueshifted Absorber"
        alpha = 0.7
        bSymbol = 'D'
        rSymbol = 'o'
        
        rdif = []
        bdif = []
        rW = []
        bW = []
        
        for d,w,i in zip(difList,lyaWList,impactList):
            # check if all the values are okay
            if isNumber(d) and isNumber(w) and isNumber(i):
                if d!=-99 and w!=-99 and i!=-99:
                
                    xVal = float(i) * float(d)
                    yVal = float(w)
                    if d>0:
                        # galaxy is behind absorber, so gas is blue shifted
                        color = 'Blue'
                        symbol = bSymbol
                        
                        bdif.append(float(d))
                        bW.append(float(w))
                        
                        if countb == 0:
                            countb +=1
                            plotb = ax.scatter(xVal,yVal,marker=symbol,c='Blue',s=50,label= labelb,\
                            alpha = alpha)
                    if d<0:
                        # gas is red shifted compared to galaxy
                        color = 'Red'
                        symbol = rSymbol
                        
                        rdif.append(float(d))
                        rW.append(float(w))
                        
                        if countr == 0:
                            countr +=1
                            plotr = ax.scatter(xVal,yVal,marker=symbol,c='Red',s=50,label= labelr,\
                            alpha = alpha)
                
                    plot1 = scatter(xVal,yVal,marker=symbol,c=color,s=50,alpha = alpha)
        
        includeHist = False

        if includeHist:
            bins = arange(-400,100,100)
            bin_means,edges,binNumber = stats.binned_statistic(array(rdif), array(rW), statistic='median', bins=bins)
            left,right = edges[:-1],edges[1:]
            X = np.array([left,right]).T.flatten()
            Y = np.array([bin_means,bin_means]).T.flatten()
            plt.plot(X,Y, c='Red',ls='dashed',lw=2,alpha=alpha,label=r'$\rm Mean EW$')
        
        
            bins = arange(0,500,100)
            bin_means,edges,binNumber = stats.binned_statistic(array(bdif), array(bW), statistic='median', bins=bins)
            left,right = edges[:-1],edges[1:]
            X = np.array([left,right]).T.flatten()
            Y = np.array([bin_means,bin_means]).T.flatten()
            plt.plot(X,Y, c='Blue',ls='dashed',lw=2,alpha=alpha,label=r'$\rm Mean EW$')

        # x-axis
        majorLocator   = MultipleLocator(50000)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(25000)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        
        # y axis
        majorLocator   = MultipleLocator(200)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(100)
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)


        xlabel(r'$\rm \rho * \Delta v ~ [kpc ~ km ~ s^{-1}]$')
        ylabel(r'$\rm Equivalent ~ Width ~ [m\AA]$')
#         legend(scatterpoints=1,prop={'size':12},loc=1)
        ax.grid(b=None,which='major',axis='both')
        ylim(0,1000)
#         xlim(-400,400)
        
        if save:
            savefig('{0}/W(impact_vel_diff).pdf'.format(saveDirectory),format='pdf',bbox_inches='tight')
        else:
            show()    

    
##########################################################################################
##########################################################################################
    # plot equivalent width as a function of likelihood for red and blue shifted
    # absorption
    #
    
    plotW_likelihood = False
    save = False
    
    if plotW_likelihood:
        fig = figure()
        ax = fig.add_subplot(111)
        countb = 0
        countr = 0
        count = -1
        labelr = 'Redshifted Absorber'
        labelb = "Blueshifted Absorber"
        alpha = 0.7
        bSymbol = 'D'
        rSymbol = 'o'
        
        rlike = []
        blike = []
        rW = []
        bW = []
        xVals = []
        yVals = []
        
        for d,w,l in zip(difList,lyaWList,likeList):
            # check if all the values are okay
            if isNumber(d) and isNumber(w) and isNumber(l):
                if d!=-99 and w!=-99 and l!=-99:
                
                    xVal = float(l)
                    yVal = float(w)
                    
                    xVals.append(xVal)
                    yVals.append(yVal)
                    
                    if d>0:
                        # galaxy is behind absorber, so gas is blue shifted
                        color = 'Blue'
                        symbol = bSymbol
                        
                        blike.append(xVal)
                        bW.append(yVal)
                        
                        if countb == 0:
                            countb +=1
                            plotb = ax.scatter(xVal,yVal,marker=symbol,c='Blue',s=50,label= labelb,\
                            alpha = alpha)
                    if d<0:
                        # gas is red shifted compared to galaxy
                        color = 'Red'
                        symbol = rSymbol
                        
                        rlike.append(xVal)
                        rW.append(yVal)
                        
                        if countr == 0:
                            countr +=1
                            plotr = ax.scatter(xVal,yVal,marker=symbol,c='Red',s=50,label= labelr,\
                            alpha = alpha)
                
                    plot1 = scatter(xVal,yVal,marker=symbol,c=color,s=50,alpha = alpha)
        
        includeHist = True

        if includeHist:
            bins = arange(0,2.0,0.4)
            bin_means,edges,binNumber = stats.binned_statistic(array(xVals), array(yVals), statistic='mean', bins=bins)
            left,right = edges[:-1],edges[1:]
            X = np.array([left,right]).T.flatten()
            Y = np.array([bin_means,bin_means]).T.flatten()
            plt.plot(X,Y, c='black',ls='solid',lw=2,alpha=alpha,label=r'$\rm Mean EW$')
    

        # x-axis
        majorLocator   = MultipleLocator(0.2)
        majorFormatter = FormatStrFormatter(r'$\rm %0.1f$')
        minorLocator   = MultipleLocator(0.1)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        
        # y axis
        majorLocator   = MultipleLocator(200)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(100)
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)


        xlabel(r'$\rm \mathcal{L}$')
        ylabel(r'$\rm Equivalent ~ Width ~ [m\AA]$')
#         legend(scatterpoints=1,prop={'size':12},loc=1)
        ax.grid(b=None,which='major',axis='both')
        ylim(0,1000)
        xlim(0,2)
#         xlim(-400,400)
        
        if save:
            savefig('{0}/W(likelihood).pdf'.format(saveDirectory),format='pdf',bbox_inches='tight')
        else:
            show()
            
    
##########################################################################################
##########################################################################################

    # plot equivalent width as a function of cos(inclination) with red and blue shifted
    # absorption represented by a color bar
    #
    
    plotW_CosInc_colorbar= False
    save = False
    
    if plotW_CosInc_colorbar:
        # colormap the velocity difference of the absorber
        averaged =[]

        blueMap = cm.Blues
        redMap = cm.Reds
        
        fig = figure()
        ax = fig.add_subplot(111)
        
        countb = 0
        countr = 0
        count = -1

        rdif = []
        rcosInc = []
        rlyaW = []
        rMaj = []
        bdif = []
        bcosInc = []
        blyaW = []
        bMaj = []
        
        for d,i,w,m in zip(difList,cosIncList,lyaWList,majList):
            # check if all the values are okay
            if isNumber(d) and isNumber(i) and isNumber(w) and isNumber(m):
                if d!=-99 and i!=-99 and w!=-99 and m!=-99:
                    if d>0:
                        # galaxy is behind absorber, so gas is blue shifted
                        color = 'Blue'
                        if countb == 0:
                            bdif.append(d)
                            bcosInc.append(i)
                            blyaW.append(w)
                            bMaj.append(m)
        #                     plotb = ax.scatter(a, w, cmap=blueMap, c=d, s=50, vmin=0, vmax=400)

                    if d<0:
                        # gas is red shifted compared to galaxy
                        color = 'Red'
                        if countr == 0:
        #                     plotr = ax.scatter(a, w, cmap=redMap, c=d, s=50, vmin=0, vmax=400)
                            rdif.append(d)
                            rcosInc.append(i)
                            rlyaW.append(w)
                            rMaj.append(m)


        print
        print 'rdif: ',rdif
        print 'average, median redshifts: ',average(rdif),', ',median(rdif)
        print 'average, median blueshifts: ',average(bdif),', ',median(bdif)

        print 'max, min red dif: ',max(rdif), ', ',min(rdif)
        print 'max, min blue dif: ',max(bdif), ', ', min(bdif)
        
        plotr = ax.scatter(rcosInc, rlyaW, cmap=redMap, c=rdif, s=50, vmin=-400, vmax=0)
        norm = matplotlib.colors.Normalize(vmin = 0, vmax = 400)
        cbarRed = plt.colorbar(plotr,cmap=redMap,orientation='vertical')
        cbarRed.set_label('galaxy - absorber velocity (km/s)')
        
        plotb = ax.scatter(bcosInc, blyaW, cmap=blueMap, c=bdif, s=50, vmin=0, vmax=400)
        cbarBlue = plt.colorbar(plotb,cmap=blueMap,orientation='vertical')
        cbarBlue.set_label('galaxy - absorber velocity (km/s)')
#       cbar.ax.set_yticklabels(ticks)
                    
        title("Equivalent Width vs Cos(inclination)")
        xlabel(r'Cos(inclination) = b/a')
        ylabel(r'Equivalent Width ($\rm m\AA$)')
#         legend(scatterpoints=1)
        ax.grid(b=None,which='major',axis='both')
        ylim(-1,1200)
        xlim(-0.02,1)

        if save:
            savefig('{0}/W(cos(inclination))_colorbar.pdf'.format(saveDirectory),format='pdf')
        else:
            show()
        
        
        
#         fig = figure()
#         ax = fig.add_subplot(211)
#         hist(rdif,color='red')
#        
#         ax = fig.add_subplot(212)
#         hist(bdif,color='blue')
#         show()


##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################


if __name__=="__main__":
    # do the work
    main()
    