#!/usr/bin/env python

'''
By David French (frenchd@astro.wisc.edu)

$Id:  plotW_Inc2.py, v 5.4 9/20/16

This is the plotW_Inc bit from histograms3.py. Now is separated, and loads in a pickle
file of the relevant data, as created by "buildDataLists.py"

Plot EW as a function of inclination, fancy(inc), cos(inc) and cos(fancy(inc)). Also plot
with a colormap.


Previous (from histograms3.py):
    Plot some stuff for the 100largest initial results

    Make plots for AAS winter 2014 poster. Uses LG_correlation_combined2.csv file

    Updated for the pilot paper (05/06/15)


v5: Updated for the final pilot paper results (12/04/15)
    No longer uses pilotData.p or whatever, just loads the files here
    
    - all the versions of W_vs_inc are here now: plotW_Inc2.py, plotW_CosInc_colorbar2.py,
    plotW_fancyInc2.py, plotW_FancyCosInc2.py, plotW_FancyCosInc_colorbar2.py,
    plotW_CosInc2.py
    
v5.1: updated for LG_correlation_combined5_8_edit2.csv and l_min = 0.001
    (02/18/2016)
    
v5.2: remake plots with v_hel instead of vcorr (4/21/16)

v5.3: remake plots with new large galaxy sample (7/13/16) -> /plots4/

v5.4: remake plots after adjusting some targets (9/20/16)
        -> LG_correlation_combined5_11_25cut_edit2.csv now ->/plots5/
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

def perc50(a):
    if len(a)>0:
        return percentile(a,50)
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

    # plot equivalent width as a function of inclination for red and blue shifted
    # absorption
    #
    
    plotW_Inc = False
    save = False
    
    if plotW_Inc:
        fig = figure()
        ax = fig.add_subplot(111)
        countb = 0
        countr = 0
        count = -1
        labelr = 'Red Shifted Absorber'
        labelb = "Blue Shifted Absorber"
        
        for d,i,w,m in zip(difList,incList,lyaWList,majList):
            # check if all the values are okay
            if isNumber(d) and isNumber(i) and isNumber(w) and isNumber(m):
                if d!=-99 and i!=-99 and w!=-99 and m!=-99:
                    if d>0:
                        # galaxy is behind absorber, so gas is blue shifted
                        color = 'Blue'
                        if countb == 0:
                            countb +=1
                            plotb = ax.scatter(i,w,c='Blue',s=50,label= labelb)
                    if d<0:
                        # gas is red shifted compared to galaxy
                        color = 'Red'
                        if countr == 0:
                            countr +=1
                            plotr = ax.scatter(i,w,c='Red',s=50,label= labelr)
                
                    plot1 = scatter(i,w,c=color,s=50)
            
#         title('W(inclination) for red vs blue shifted absorption')
        xlabel(r'Inclination (deg)')
        ylabel(r'Equivalent Width ($\rm m\AA$)')
        legend(scatterpoints=1,prop={'size':12},loc=2)
        ax.grid(b=None,which='major',axis='both')
        ylim(-1,1200)
        xlim(0,90)
        
        if save:
            savefig('{0}/W(inclination)_dif.pdf'.format(saveDirectory),format='pdf')
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
        
        allLyaW = []
        allInc = []
        allDif = []
        
        for d,i,w,m in zip(difList,cosIncList,lyaWList,majList):
            # check if all the values are okay
            if isNumber(d) and isNumber(i) and isNumber(w) and isNumber(m):
                if d!=-99 and i!=-99 and w!=-99 and m!=-99:
                    allLyaW.append(w)
                    allInc.append(i)
                    allDif.append(d)
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
        
#         plotr = ax.scatter(rcosInc, rlyaW, cmap=redMap, c=rdif, s=50, vmin=-400, vmax=0)
#         norm = matplotlib.colors.Normalize(vmin = 0, vmax = 400)
#         cbarRed = plt.colorbar(plotr,cmap=redMap,orientation='vertical')
#         cbarRed.set_label('galaxy - absorber velocity (km/s)')
#         
#         plotb = ax.scatter(bcosInc, blyaW, cmap=blueMap, c=bdif, s=50, vmin=0, vmax=400)
#         cbarBlue = plt.colorbar(plotb,cmap=blueMap,orientation='vertical')
#         cbarBlue.set_label('galaxy - absorber velocity (km/s)')
#       cbar.ax.set_yticklabels(ticks)
                    
        plot = ax.scatter(allInc, allLyaW, cmap=cm.RdBu, c=allDif, s=50, vmin=-400, vmax=400)
        norm = matplotlib.colors.Normalize(vmin = -400, vmax = 400)
        cbarRed = plt.colorbar(plot,cmap=cm.RdBu,orientation='vertical')
        cbarRed.set_label('galaxy - absorber velocity (km/s)')
        
#         plotb = ax.scatter(bcosInc, blyaW, cmap=cm.RdBu, c=bdif, s=50, vmin=0, vmax=400)
#         cbarBlue = plt.colorbar(plotb,cmap=cm.RdBu,orientation='vertical')
#         cbarBlue.set_label('galaxy - absorber velocity (km/s)')
                    

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
    # plot equivalent width as a function of fancy_inclination for red and blue shifted
    # absorption
    #
    
    plotW_fancyInc = False
    save = False
    alpha = 0.75
    
    if plotW_fancyInc:
        fig = figure()
        ax = fig.add_subplot(111)
        countb = 0
        countr = 0
        count = -1
        labelr = 'Red Shifted Absorber'
        labelb = "Blue Shifted Absorber"
        for d,i,w,m in zip(difList,fancyIncList,lyaWList,majList):
            # check if all the values are okay
            if isNumber(d) and isNumber(i) and isNumber(w) and isNumber(m):
                if d!=-99 and i!=-99 and w!=-99 and m!=-99:
                    if i ==90:
                        print d, i, w, m
                    if d>0:
                        # galaxy is behind absorber, so gas is blue shifted
                        color = 'Blue'
                        if countb == 0:
                            countb +=1
                            plotb = ax.scatter(i,w,c='Blue',s=50,label= labelb,alpha=alpha)
                    if d<0:
                        # gas is red shifted compared to galaxy
                        color = 'Red'
                        if countr == 0:
                            countr +=1
                            plotr = ax.scatter(i,w,c='Red',s=50,label= labelr,alpha=alpha)
                
                    plot1 = scatter(i,w,c=color,s=50,alpha=alpha)
            
#         title('W(fancy_inclination) for red vs blue shifted absorption')

        majorLocator   = MultipleLocator(10)
        majorFormatter = FormatStrFormatter('%d')
        minorLocator   = MultipleLocator(2)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        
        majorLocator   = MultipleLocator(200)
        majorFormatter = FormatStrFormatter('%d')
        minorLocator   = MultipleLocator(50)
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)

        xlabel(r'$\rm Inclination (deg)$')
        ylabel(r'$ \rm Equivalent Width ($\rm m\AA$)$')
        legend(scatterpoints=1,prop={'size':12},loc=2)
        ax.grid(b=None,which='major',axis='both')
        ylim(-1,1200)
        xlim(0,90)
        
        if save:
            savefig('{0}/W(fancy_inclination)_dif2.pdf'.format(saveDirectory),format='pdf')
        else:
            show()
            
##########################################################################################
##########################################################################################
    # plot equivalent width as a function of fancy_inclination for red and blue shifted
    # absorption and include the dashed line w/ equation
    #
    
    plotW_fancyInc_line = False
    save = False
    
    if plotW_fancyInc_line:
        fig = figure()
        ax = fig.add_subplot(111)
        countb = 0
        countr = 0
        count = -1
        labelr = 'Red Shifted Absorber'
        labelb = "Blue Shifted Absorber"
        alpha = 0.85

        for d,i,w,m in zip(difList,fancyIncList,lyaWList,majList):
            # check if all the values are okay
            if isNumber(d) and isNumber(i) and isNumber(w) and isNumber(m):
                if d!=-99 and i!=-99 and w!=-99 and m!=-99:
                    if i ==90:
                        print d, i, w, m
                    if d>0:
                        # galaxy is behind absorber, so gas is blue shifted
                        color = 'Blue'
                        if countb == 0:
                            countb +=1
                            plotb = ax.scatter(i,w,c='Blue',s=50,label= labelb,alpha=alpha)
                    if d<0:
                        # gas is red shifted compared to galaxy
                        color = 'Red'
                        if countr == 0:
                            countr +=1
                            plotr = ax.scatter(i,w,c='Red',s=50,label= labelr,alpha=alpha)
                
                    plot1 = scatter(i,w,c=color,s=50,alpha=alpha)
            
#         title('W(fancy_inclination) for red vs blue shifted absorption')

        def lineFunc(yInt,slope):
            
            ys = []
            xs = []
            for x in arange(0,91,1):
                y = slope* x + yInt
                ys.append(y)
                xs.append(x)
        
            return xs,ys
        
        # steeper line
#         yInt = -100
#         slope = round(500./79.,2)
        
        # shallower line
        yInt = -25
        slope = round(500./129.,2)
        
        lineX,lineY = lineFunc(yInt,slope)
        plot2 = plot(lineX,lineY,c='Black',ls='dashed',alpha=alpha,lw=2)
        
        print
        print 'y = {0}*inc + {1}'.format(slope,yInt)
        print
        
        xlabel(r'Inclination (deg)')
        ylabel(r'Equivalent Width ($\rm m\AA$)')
        legend(scatterpoints=1,prop={'size':12},loc=2)
        ax.grid(b=None,which='major',axis='both')
        ylim(-1,1200)
        xlim(0,90)
        
        if save:
            savefig('{0}/W(fancy_inclination)_line_shallow.pdf'.format(saveDirectory),format='pdf')
        else:
            show()
            
            
##########################################################################################
##########################################################################################
    # plot equivalent width as a function of cos(fancy_inclination) for red and blue shifted
    # absorption
    #
    
    plotW_CosFancyInc = False
    save = False
    
    if plotW_CosFancyInc:
        fig = figure()
        ax = fig.add_subplot(111)
        countb = 0
        countr = 0
        count = -1
        labelr = 'Red Shifted Absorber'
        labelb = "Blue Shifted Absorber"
        for d,i,w,m in zip(difList,cosFancyIncList,lyaWList,majList):
        
            # check if all the values are good
            if isNumber(d) and isNumber(i) and isNumber(w) and isNumber(m):
                if d!=-99 and i!=-99 and w!=-99 and m!=-99:
                    if d>0:
                        # galaxy is behind absorber, so gas is blue shifted
                        color = 'Blue'
                        if countb == 0:
                            countb +=1
                            plotb = ax.scatter(i,w,c='Blue',s=50,label= labelb)
                    if d<0:
                        # gas is red shifted compared to galaxy
                        color = 'Red'
                        if countr == 0:
                            countr +=1
                            plotr = ax.scatter(i,w,c='Red',s=50,label= labelr)
                
                    plot1 = scatter(i,w,c=color,s=50)
            
        title('W(cos(fancy_inclination)) for red vs blue shifted absorption')
        xlabel(r'Cos(fancy_inclination)')
        ylabel(r'Equivalent Width ($\rm m\AA$)')
        legend(scatterpoints=1)
        ax.grid(b=None,which='major',axis='both')
        ylim(-1,1200)
        xlim(-0.01,1)
        
        if save:
            savefig('{0}/W(cos(fancy_inclination))_dif.pdf'.format(saveDirectory),format='pdf')
        else:
            show()
            
            
##########################################################################################
##########################################################################################

    # plot equivalent width as a function of cos(fancy_inclination) with velocity dif
    # represented by a color bar
    #
    
    plotW_FancyCosInc_colorbar= False
    save = False
    
    if plotW_FancyCosInc_colorbar:
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
        
        for d,i,w,m in zip(difList,cosFancyIncList,lyaWList,majList):
            # check if all the values are okay
            if isNumber(d) and isNumber(i) and isNumber(w) and isNumber(m):
                if d !=-99 and i !=-99 and w!=-99 and m!=-99:
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
                    
        title('Equivalent width vs Cos(fancy_inclination) colormap')
        xlabel(r'Cos(fancy_inclination) = b/a')
        ylabel(r'Equivalent Width ($\rm m\AA$)')
#         legend(scatterpoints=1)
        ax.grid(b=None,which='major',axis='both')
        ylim(-1,1200)
        xlim(-0.02,1)

        if save:
            savefig('{0}/W(cos(fancy_inclination))_colorbar.pdf'.format(saveDirectory),format='pdf')
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

    # plot equivalent width as a function of cos(inclination) for red and blue shifted
    # absorption
    #
    
    plotW_CosInc = False
    save = False
    
    if plotW_CosInc:
    
        fig = figure()
        ax = fig.add_subplot(111)
        countb = 0
        countr = 0
        count = -1
        labelr = 'Red Shifted Absorber'
        labelb = "Blue Shifted Absorber"
        for d,i,w,m in zip(difList,cosIncList,lyaWList,majList):
        
            # check if all the values are good
            if isNumber(d) and isNumber(i) and isNumber(w) and isNumber(m):
                if d!=-99 and i!=-99 and w!=-99 and m!=-99:
                    if d>0:
                        # galaxy is behind absorber, so gas is blue shifted
                        color = 'Blue'
                        if countb == 0:
                            countb +=1
                            plotb = ax.scatter(i,w,c='Blue',s=50,label= labelb)
                    if d<0:
                        # gas is red shifted compared to galaxy
                        color = 'Red'
                        if countr == 0:
                            countr +=1
                            plotr = ax.scatter(i,w,c='Red',s=50,label= labelr)
                
                    plot1 = scatter(i,w,c=color,s=50)
            
        title('W(cos(inclination)) for red vs blue shifted absorption')
        xlabel(r'Cos(inclination) = b/a')
        ylabel(r'Equivalent Width ($\rm m\AA$)')
        legend(scatterpoints=1)
        ax.grid(b=None,which='major',axis='both')
        ylim(-1,1200)
        xlim(0,1)
        
        if save:
            savefig('{0}/W(cos(inclination))_dif.pdf'.format(saveDirectory),format='pdf')
        else:
            show()



##########################################################################################
##########################################################################################

    # plot equivalent width as a function of cos(inclination) for red and blue shifted
    # absorption, separating azimuth >45 and <45 absorbers
    #
    
    plotW_inc_az = False
    save = False
    
    if plotW_inc_az:
    
        fig = figure()
        ax = fig.add_subplot(111)
        countb = 0
        countr = 0
        count = -1
        labelr = 'Red Shifted Absorber'
        labelb = "Blue Shifted Absorber"
        labelrless = 'Red Shifted Absorber <45 az'
        labelbless = "Blue Shifted Absorber <45 az"
        
        # limit how far away the absorbers can be: impact / virial <= distLimit
        distLimit = 1.4
        
        rless = []
        rmore = []
        bless = []
        bmore = []
        
        for d,i,w,a,v,imp in zip(difList,fancyIncList,lyaWList,azList,virList,impactList):
        
            # check if all the values are good
            if isNumber(d) and isNumber(i) and isNumber(w) and isNumber(a) and isNumber(v) and isNumber(imp):
                if d!=-99 and i!=-99 and w!=-99 and a!=-99 and float(imp)/float(v)<= distLimit:
                    if d>0:
                        # galaxy is behind absorber, so gas is blue shifted
                        color = 'Blue'

                        if a >=45:
                            plotb = ax.scatter(i,w,c='Blue',s=50,label= labelb)
                            bmore.append(float(w))
                        else:
                            plotb = ax.scatter(i,w,c='Blue',s=50,label= labelbless,marker='*',lw=0)
                            bless.append(float(w))
                    if d<0:
                        # gas is red shifted compared to galaxy
                        color = 'Red'

                        if a >=45:
                            plotr = ax.scatter(i,w,c='Red',s=50,label= labelr)
                            rmore.append(float(w))
                        else:
                            plotr = ax.scatter(i,w,c='Red',s=50,label= labelrless,marker='*',lw=0)
                            rless.append(float(w))

                
#                     plot1 = scatter(i,w,c=color,s=10)
            
        print 'average bmore: ',mean(bmore)
        print 'average bless: ',mean(bless)
        print
        print 'average rmore: ',mean(rmore)
        print 'average rless: ',mean(rless)


        xlabel(r'Inclination')
        ylabel(r'Equivalent Width ($\rm m\AA$)')
#         legend(scatterpoints=1)
        ax.grid(b=None,which='major',axis='both')
#         ylim(-1,1200)
#         xlim(0,1)
        
        if save:
            savefig('{0}/W(inclination)_dif_az.pdf'.format(saveDirectory),format='pdf')
        else:
            show()

##########################################################################################
##########################################################################################
    # plot equivalent width as a function of fancy_inclination for red and blue shifted
    # absorption, include an overlying histogram
    #
    #
    # this is not very enlightening
    
    plotW_inc_hist = False
    save = False
    
    if plotW_inc_hist:
        fig = figure()
        ax = fig.add_subplot(111)
        
        countb = 0
        countr = 0
        count = -1
        binSize = 12
        alpha = 0.75
        
        labelr = 'Redshifted Absorber'
        labelb = "Blueshifted Absorber"
        
        placeArrayr = zeros(int(90/binSize)+1)
        placeCountr = zeros(int(90/binSize)+1)
        placeArrayb = zeros(int(90/binSize)+1)
        placeCountb = zeros(int(90/binSize)+1)
        
        placeArray = zeros(int(90/binSize)+1)
        placeCount = zeros(int(90/binSize)+1)
        
        for d,i,w,v in zip(difList,fancyIncList,lyaWList,virList):
            # check if all the values are good
            if isNumber(d) and isNumber(i) and isNumber(w) and isNumber(v):
                if d!=-99 and i!=-99 and w!=-99 and v!=-99:
                
                    # do totals first
                    place = i/binSize
                    placeArray[place] += float(w)
                    placeCount[place] +=1.
                        
                    if d>0:
                        # galaxy is behind absorber, so gas is blue shifted
                        color = 'Blue'
                        
                        # which bin does it belong too?
                        place = i/binSize
                        placeArrayb[place] += float(w)
                        placeCountb[place] +=1.
                        
                        if countb == 0:
                            countb +=1
#                             plotb = ax.scatter(v,w,c='Blue',s=50,label= labelb)
                            plotb = ax.scatter(i,w,c='Blue',s=50,alpha = alpha)

                    if d<0:
                        # gas is red shifted compared to galaxy
                        color = 'Red'
                        
                        # which bin does it belong too?
                        place = i/binSize
                        placeArrayr[place] += float(w)
                        placeCountr[place] +=1.
                        
                        if countr == 0:
                            countr +=1
#                             plotr = ax.scatter(v,w,c='Red',s=50,label= labelr)
                            plotr = ax.scatter(i,w,c='Red',s=50,alpha = alpha)

                    plot1 = scatter(i,w,c=color,s=50,alpha=alpha)
        
#         print
#         print 'placeArrayTotal: ',placeArrayTotal
#         print 'placeCountTotal: ',placeCountTotal
#         print
#         print 'placeArrayr: ',placeArrayr
#         print 'placeCountr: ',placeArrayr
#         print
        
        rHist = placeArrayr/placeCountr
        bHist = placeArrayb/placeCountb
        allHist = placeArray/placeCount
        
#         print 'totalHist: ',totalHist
        
        totalrHist = []
        totalrVir = []
        totalbHist = []
        totalbVir = []
        totalHist = []
        totalVir = []
        
        for r,v in zip(rHist,arange(0,100,binSize)):
            if not isNumber(r):
                r = 0
            totalrHist.append(r)
            totalrHist.append(r)

            totalrVir.append(v)
            totalrVir.append(v+binSize)
            
        for b,v in zip(bHist,arange(0,100,binSize)):
            if not isNumber(b):
                b = 0
            totalbHist.append(b)
            totalbHist.append(b)

            totalbVir.append(v)
            totalbVir.append(v+binSize)
            
            
        print 'hello'
        print 'allHist: ',allHist
        for h,v in zip(allHist,arange(0,100,binSize)):
            if not isNumber(h):
                h = 0
                
            print 'h: ',h
            print 'v: ',v
            print
            
            totalHist.append(h)
            totalHist.append(h)

            totalVir.append(v)
            totalVir.append(v+binSize)
            
        print 'hello 2'
        
#         print 'totalrVir: ',totalrVir
#         print 'totalrHist: ',totalrHist
#         print
#         print 'totalbVir: ',totalbVir
#         print 'totalbHist: ',totalbHist
#         print
#         print 'totalVir: ',totalVir
#         print 'totalHist: ',totalHist
        
        plot2 = ax.plot(totalrVir,totalrHist,c='Red',lw=2,ls='dashed',\
        label='Average Redshifted EW',alpha=alpha)
        
        plot3 = ax.plot(totalbVir,totalbHist,c='Blue',lw=2,ls='dotted',\
        label='Average Blueshifted EW',alpha=alpha)

        plot4 = ax.plot(totalVir,totalHist,c='Black',lw=2,ls='dashed',\
        label='Average Total EW',alpha=alpha)
        
        xlabel(r'$\rm Inclination$ (deg)')
        ylabel(r'Equivalent Width ($\rm m\AA$)')
        ax.legend(scatterpoints=1,prop={'size':12},loc=2)
        ax.grid(b=None,which='major',axis='both')
        ylim(-5,max(lyaWList)+100)
        xlim(0,100)

        if save:
            savefig('{0}/W(inc)_avgHistograms.pdf'.format(saveDirectory),format='pdf')
        else:
            show()



##########################################################################################
##########################################################################################
    # plot equivalent width as a function of fancy_inclination for red and blue shifted
    # absorption, include an overlying median histograms
    #
    #
    
    plotW_inc_hist_median = False
    save = False
    
    if plotW_inc_hist_median:
        fig = figure()
        ax = fig.add_subplot(111)
        
        countb = 0
        countr = 0
        count = -1
        binSize = 10
        numBins = 9
        alpha = 0.75
        bins = [0,15,30,45,60,75,90,105]
#         bins = arange(0,100,10)

        labelr = r'$\rm Redshifted$'
        labelb = r'$\rm Blueshifted$'
        
        allInc = []
        allW = []
        allVir = []
        
        rInc = []
        bInc = []
        rW = []
        bW = []
        rVir = []
        bVir = []
        
        for d,i,w,v in zip(difList,fancyIncList,lyaWList,virList):
            # check if all the values are good
            if isNumber(d) and isNumber(i) and isNumber(w) and isNumber(v):
                if d!=-99 and i!=-99 and w!=-99 and v!=-99:
                    
                    allInc.append(float(i))
                    allW.append(float(w))
                    allVir.append(float(v))
                        
                    if d>0:
                        # galaxy is behind absorber, so gas is blue shifted
                        color = 'Blue'
                        
                        bInc.append(float(i))
                        bW.append(float(w))
                        bVir.append(float(v))
                                            
                        if countb == 0:
                            countb +=1
#                             plotb = ax.scatter(v,w,c='Blue',s=50,label= labelb)
                            plotb = ax.scatter(i,w,c='Blue',s=50,alpha = alpha)

                    if d<0:
                        # gas is red shifted compared to galaxy
                        color = 'Red'
                        
                        rInc.append(float(i))
                        rW.append(float(w))
                        rVir.append(float(v))
                        
                        if countr == 0:
                            countr +=1
#                             plotr = ax.scatter(v,w,c='Red',s=50,label= labelr)
                            plotr = ax.scatter(i,w,c='Red',s=50,alpha = alpha)

                    plot1 = scatter(i,w,c=color,s=50,alpha=alpha)
                            
        print 'rInc: ,',rInc
        print
        print 'bInc: ',bInc

        # totals
        bin_means,edges,binNumber = stats.binned_statistic(array(allInc), array(allW), statistic='mean', bins=bins)
        left,right = edges[:-1],edges[1:]
        X = np.array([left,right]).T.flatten()
        Y = np.array([bin_means,bin_means]).T.flatten()
        plt.plot(X,Y, c='Black',ls='dashed',lw=2,alpha=alpha,label=r'$\rm Average ~ EW$')
        
        # red shifted
#         bin_means,edges,binNumber = stats.binned_statistic(array(rInc), array(rW), statistic='mean', bins=bins)
#         left,right = edges[:-1],edges[1:]
#         X = np.array([left,right]).T.flatten()
#         Y = np.array([bin_means,bin_means]).T.flatten()
#         plt.plot(X,Y, c='red',ls='dotted',lw=2,alpha=alpha,label="Mean redshifted EW")
        
        # blue shifted
#         bin_means,edges,binNumber = stats.binned_statistic(array(bInc), array(bW), statistic='mean', bins=bins)
#         left,right = edges[:-1],edges[1:]
#         X = np.array([left,right]).T.flatten()
#         Y = np.array([bin_means,bin_means]).T.flatten()
#         plt.plot(X,Y, c='blue',ls='dotted',lw=2,alpha=alpha,label="Mean blueshifted EW")
        
#         plot2 = ax.plot(totalrVir,totalrHist,c='Red',lw=2,ls='dashed',\
#         label='Average Redshifted EW',alpha=alpha)
#         
#         plot3 = ax.plot(totalbVir,totalbHist,c='Blue',lw=2,ls='dotted',\
#         label='Average Blueshifted EW',alpha=alpha)

#         plot4 = ax.plot(totalVir,totalHist,c='Black',lw=2,ls='dashed',\
#         label='Average Total EW',alpha=alpha)

        # x axis
        majorLocator   = MultipleLocator(10)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(5)
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
        
        
        xlabel(r'$\rm Galaxy ~ Inclination ~ [deg]$')
        ylabel(r'$\rm Equivalent ~ Width ~ [m\AA]$')
        ax.legend(scatterpoints=1,prop={'size':14},loc=2,fancybox=True)
        ax.grid(b=None,which='major',axis='both')
#         ylim(-5,max(lyaWList)+100)
        ylim(0,1200)
        xlim(0,91)

        if save:
            savefig('{0}/W(inc)_medHistogram.pdf'.format(saveDirectory),format='pdf',bbox_inches='tight')
        else:
            show()
            
            
        
        includeHist = False
        if includeHist:
        
            fig = figure(figsize=(10,4))
            ax = fig.add_subplot(111)
            bins = [0,15,30,45,60,75,90,105]

            bin_means,edges,binNumber = stats.binned_statistic(array(allInc), array(allW), statistic='count', bins=bins)
    
            left,right = edges[:-1],edges[1:]
            X = np.array([left,right]).T.flatten()
            Y = np.array([bin_means,bin_means]).T.flatten()

            plt.plot(X,Y, c='Black',ls='dashed',lw=2,alpha=alpha,label="Mean EW")
            show()




##########################################################################################
##########################################################################################
    # plot equivalent width as a function of fancy_inclination for red and blue shifted
    # absorption, include an overlying 10th, 50th and 90th percentile histograms
    #
    #
    
    plotW_inc_percentile = True
    save = True
    
    if plotW_inc_percentile:
        fig = figure(figsize=(7.7,5.7))
        ax = fig.add_subplot(111)
        
        countb = 0
        countr = 0
        count = -1
        binSize = 10
        numBins = 9
        alpha = 0.7

        binSize = 12
        bins = arange(0,90+binSize,binSize)

        labelr = r'$\rm Redshifted$'
        labelb = r'$\rm Blueshifted$'
        bSymbol = 'D'
        rSymbol = 'o'
        
        allInc = []
        allW = []
        allVir = []
        
        rInc = []
        bInc = []
        rW = []
        bW = []
        rVir = []
        bVir = []
        
        for d,i,w,v in zip(difList,fancyIncList,lyaWList,virList):
            # check if all the values are good
            if isNumber(d) and isNumber(i) and isNumber(w) and isNumber(v):
                if d!=-99 and i!=-99 and w!=-99 and v!=-99:
                    
                    allInc.append(float(i))
                    allW.append(float(w))
                    allVir.append(float(v))
                        
                    if d>0:
                        # galaxy is behind absorber, so gas is blue shifted
                        color = 'Blue'
                        symbol = bSymbol
                        
                        bInc.append(float(i))
                        bW.append(float(w))
                        bVir.append(float(v))
                                            
                        if countb == 0:
                            countb +=1
#                             plotb = ax.scatter(v,w,c='Blue',s=50,label= labelb)
                            plotb = ax.scatter(i,w,marker=bSymbol,c='Blue',s=50,alpha=alpha)

                    if d<0:
                        # gas is red shifted compared to galaxy
                        color = 'Red'
                        symbol = rSymbol
                        
                        rInc.append(float(i))
                        rW.append(float(w))
                        rVir.append(float(v))
                        
                        if countr == 0:
                            countr +=1
#                             plotr = ax.scatter(v,w,c='Red',s=50,label= labelr)
                            plotr = ax.scatter(i,w,marker=rSymbol,c='Red',s=50,alpha=alpha)

                    plot1 = scatter(i,w,marker=symbol,c=color,s=50,alpha=alpha)
                            

#         totals
#         bin_means,edges,binNumber = stats.binned_statistic(array(allInc), array(allW), statistic='mean', bins=bins)
#         left,right = edges[:-1],edges[1:]
#         X = np.array([left,right]).T.flatten()
#         Y = np.array([bin_means,bin_means]).T.flatten()
#         plt.plot(X,Y, c='Black',ls='dashed',lw=2,alpha=alpha,label=r'$\rm Average ~ EW$')

        # mean
        bin_means,edges,binNumber = stats.binned_statistic(array(allInc), array(allW), \
        statistic='mean', bins=bins)
        left,right = edges[:-1],edges[1:]        
        X = array([left,right]).T.flatten()
        Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
        plt.plot(X,Y, ls='solid',color='black',lw=2.0,alpha=alpha+0.1,label=r'$\rm Mean ~EW$')
        
        # 90% percentile
        bin_means,edges,binNumber = stats.binned_statistic(array(allInc), array(allW), \
        statistic=lambda y: perc90(y), bins=bins)
        left,right = edges[:-1],edges[1:]        
        X = array([left,right]).T.flatten()
        Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
        plt.plot(X,Y, ls='dashed',color='dimgrey',lw=2.0,alpha=alpha+0.1,label=r'$\rm 90th\% ~EW$')

#         # 50% percentile
#         bin_means,edges,binNumber = stats.binned_statistic(array(allInc), array(allW), \
#         statistic=lambda y: perc50(y), bins=bins)
#         left,right = edges[:-1],edges[1:]        
#         X = array([left,right]).T.flatten()
#         Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
#         plt.plot(X,Y, ls='solid',color='black',lw=2.0,alpha=alpha,label=r'$\rm Median ~EW$')
        
        
#         10th percentile
#         bin_means,edges,binNumber = stats.binned_statistic(array(allInc), array(allW), \
#         statistic=lambda y: perc10(y), bins=bins)
#         left,right = edges[:-1],edges[1:]        
#         X = array([left,right]).T.flatten()
#         Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
#         plt.plot(X,Y, ls='dashed',color='green',lw=1.5,alpha=alpha,label=r'$\rm 10th\% ~EW$')


        # x axis
        majorLocator   = MultipleLocator(10)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(5)
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
        
        
        xlabel(r'$\rm Galaxy ~ Inclination ~ [deg]$')
        ylabel(r'$\rm Equivalent ~ Width ~ [m\AA]$')
        ax.legend(scatterpoints=1,prop={'size':15},loc=2,fancybox=True)
        ax.grid(b=None,which='major',axis='both')
#         ylim(-5,max(lyaWList)+100)
        ylim(-2,1000)
        xlim(0,90.5)

        if save:
            savefig('{0}/W(fancy_inc)_mean_90_hist.pdf'.format(saveDirectory),format='pdf',bbox_inches='tight')
        else:
            show()
            
            
        
        includeHist = False
        if includeHist:
        
            fig = figure(figsize=(10,4))
            ax = fig.add_subplot(111)
            bins = [0,15,30,45,60,75,90,105]

            bin_means,edges,binNumber = stats.binned_statistic(array(allInc), array(allW), statistic='count', bins=bins)
    
            left,right = edges[:-1],edges[1:]
            X = np.array([left,right]).T.flatten()
            Y = np.array([bin_means,bin_means]).T.flatten()

            plt.plot(X,Y, c='Black',ls='dashed',lw=2,alpha=alpha,label="Mean EW")
            show()



##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################


if __name__=="__main__":
    # do the work
    main()
    