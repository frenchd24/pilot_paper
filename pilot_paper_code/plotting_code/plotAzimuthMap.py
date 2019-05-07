#!/usr/bin/env python

'''
By David French (frenchd@astro.wisc.edu)

$Id:  plotAzimuthMap.py, v 1.0 09/19/2016

Plot absorption properties on a map showing their distribution around a central galaxy


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
from matplotlib.patches import Ellipse


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
        WS09data = '/Users/David/Research_Documents/inclination/git_inclination/WS2009_lya_data.tsv'

    elif getpass.getuser() == 'frenchd':
        pickleFilename = '/usr/users/frenchd/inclination/git_inclination/pilot_paper_code/pilotData2.p'
#         resultsFilename = '/usr/users/frenchd/inclination/git_inclination/LG_correlation_combined5_8_edit2.csv'
#         saveDirectory = '/usr/users/frenchd/inclination/git_inclination/pilot_paper_code/plots2/'
        resultsFilename = '/usr/users/frenchd/inclination/git_inclination/LG_correlation_combined5_11_25cut_edit4.csv'
        saveDirectory = '/usr/users/frenchd/inclination/git_inclination/pilot_paper_code/plots5/'
        WS09data = '/usr/users/frenchd/inclination/git_inclination/WS2009_lya_data.tsv'

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
    
    WS = open(WS09data,'rU')
    WSreader = csv.DictReader(WS,delimiter=';')
    
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
    raList = []
    decList = []
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
    
    
    # WS lists
    WSvcorr = []
    WSdiam = []
    WSimpact =[]
    WSew = []
    WSvel = []
    WSlya = []
    WSvel_dif = []
    WSvir = []
    WSlike = []
    
    l_min = 0.001

    for w in WSreader:
        vcorr = w['HV']
        diam = w['Diam']
        rho = w['rho']
        ew = w['EWLya']
        vel = w['LyaVel']
        lya = w['Lya']
        
        if lya == 'Lya  ' and isNumber(diam) and isNumber(ew) and isNumber(rho):
            if float(rho) <=500.0:
                # this is a single galaxy association
                vir = calculateVirialRadius(float(diam))
                
                vel_dif = float(vcorr) - float(vel)
    
                # try this "sphere of influence" value instead
                m15 = float(diam)**1.5

                # first for the virial radius
                likelihood = math.exp(-(float(rho)/vir)**2) * math.exp(-(vel_dif/200.)**2)
                
                if vir>= float(rho):
                    likelihood = likelihood*2
                    
                # then for the second 'virial like' m15 radius
                likelihoodm15 = math.exp(-(float(rho)/m15)**2) * math.exp(-(vel_dif/200.)**2)
                
                if m15>= float(rho):
                    likelihoodm15 = likelihoodm15*2
                    
                if likelihood <= likelihoodm15:
                    likelihood = likelihoodm15
                    
                WSlike.append(likelihood)
                
#                 l_min=0
                
                if likelihood >= l_min:
                
                    WSvcorr.append(float(vcorr))
                    WSdiam.append(float(diam))
                    WSvir.append(vir)
                    WSimpact.append(float(rho))
                    WSew.append(float(ew))
                    WSvel.append(float(vel))
                    WSlya.append(lya)
                    WSvel_dif.append(vel_dif)
    
    
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
            
            # all the lists to be used for associated lines
            if float(env) <= maxEnv and float(likelihood) >= minL:
                raList.append(galaxyRA_Dec[0])
                decList.append(galaxyRA_Dec[1])
                lyaVList.append(float(lyaV))
                lyaWList.append(float(lyaW))
                lyaErrList.append(float(lyaW_err))
                naList.append(na)
                bList.append(float(b))
                impactList.append(float(impact))
                print az
                azList.append(float(az))
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
    WS.close()
            
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

    # plot histograms of the inclinations for both associated galaxies and the 
    # full galaxy data set, combining both redshifted and blueshifted
    plot_edge_on = True
    plotVertical = False
    save = True
    
    if plot_edge_on:
        
        azBlue_1 = []
        azRed_1 = []
        azBlue_2 = []
        azRed_2 = []
        azBlue_3 = []
        azRed_3 = []
        
        # new coordinates of absorption for plotting (w/ galaxy at (0,0))
        yList_1r = []
        xList_1r = []
        yList_2r = []
        xList_2r = []
        yList_3r = []
        xList_3r = []
        
        yList_1b = []
        xList_1b = []
        yList_2b = []
        xList_2b = []
        yList_3b = []
        xList_3b = []
        
        
        # calculate the position on the sky for each absorption feature wrt to the galaxy
        for r,d,i,a,fInc,maj,vir,dif,inc in zip(raList,decList,impactList,azList,fancyIncList,majList,virList,difList,incList):
            if float(a) >= 0 and float(fInc)>=0 and float(fInc)<=90:
                if inc <=40:
                    # y coordinate
                    y = (float(i)/float(vir)) * sin((a*pi)/180.)
                
                    # x coordinate
                    x = (float(i)/float(vir)) * cos(a*pi/180.)
                    
                    if dif >0:
                        #blue absorber
                        azBlue_1.append(a)
                        xList_1b.append(x)
                        yList_1b.append(y)

                    else:
                        azRed_1.append(a)
                        xList_1r.append(x)
                        yList_1r.append(y)

                if inc > 40 and inc <=65:
                    # y coordinate
                    y = (float(i)/float(vir)) * sin((a*pi)/180.)
                
                    # x coordinate
                    x = (float(i)/float(vir)) * cos(a*pi/180.)
                    
                    if dif >0:
                        #blue absorber
                        azBlue_2.append(a)
                        yList_2b.append(y)
                        xList_2b.append(x)

                    else:
                        azRed_2.append(a)
                        yList_2r.append(y)
                        xList_2r.append(x)

                if inc > 65:
                    # y coordinate
                    y = (float(i)/float(vir)) * sin((a*pi)/180.)
                
                    # x coordinate
                    x = (float(i)/float(vir)) * cos(a*pi/180.)
                    
                    if dif >0:
                        #blue absorber
                        azBlue_3.append(a)
                        yList_3b.append(y)
                        xList_3b.append(x)

                    else:
                        azRed_3.append(a)
                        yList_3r.append(y)
                        xList_3r.append(x)
                                        
            else:
                print 'float(a) <0: ',r,d,i,a,fInc

        # calculate the average red vs blue azimuth line
        blueAvg1 = mean(azBlue_1)
        print 'blueAvg1: ',blueAvg1
        redAvg1 = mean(azRed_1)
        print 'redAvg1: ',redAvg1
        print
        blueAvg2 = mean(azBlue_2)
        print 'blueAvg2: ',blueAvg2
        redAvg2 = mean(azRed_2)
        print 'redAvg2: ',redAvg2
        print
        blueAvg3 = mean(azBlue_3)
        print 'blueAvg3: ',blueAvg3
        redAvg3 = mean(azRed_3)
        print 'redAvg3: ',redAvg3
        print
        
        xyBlueAvg1 = (500.* cos(blueAvg1 * pi/180.), 500.* sin(blueAvg1 * pi/180.))
        xyRedAvg1 = (500.* cos(redAvg1 * pi/180.), 500.* sin(redAvg1 * pi/180.))
        print 'xyBlueAvg1: ',xyBlueAvg1
        print 'xyRedAvg1: ',xyRedAvg1
        xyBlueAvg2 = (500.* cos(blueAvg2 * pi/180.), 500.* sin(blueAvg2 * pi/180.))
        xyRedAvg2 = (500.* cos(redAvg2 * pi/180.), 500.* sin(redAvg2 * pi/180.))
        print 'xyBlueAvg2: ',xyBlueAvg2
        print 'xyRedAvg2: ',xyRedAvg2
        xyBlueAvg3 = (500.* cos(blueAvg3 * pi/180.), 500.* sin(blueAvg3 * pi/180.))
        xyRedAvg3 = (500.* cos(redAvg3 * pi/180.), 500.* sin(redAvg3 * pi/180.))
        print 'xyBlueAvg3: ',xyBlueAvg3
        print 'xyRedAvg3: ',xyRedAvg3


        if plotVertical:
        
            # plot the distributions 
            fig = figure(figsize=(4,10))
            alpha = 0.7
            bSymbol = 'D'
            rSymbol = 'o'
            subplots_adjust(left=None,bottom=None,right=None,top=None,wspace=None,hspace=0.001)
        
            ax1 = fig.add_subplot(311)
            # Low inclinations = circular mostly
            e = Ellipse(xy=(0,0), width=1.0, height=0.9, angle=0)
    #         legend(scatterpoints=1,prop={'size':15},loc=2,fancybox=True)
                        
            # no transparency
            e.set_alpha(0.3)
            ax1.add_artist(e)
            e.set_facecolor('black')
            e.set_edgecolor('black')
        
            maxW = 120
            minW = 5
            minLya = min(lyaWList)
            maxLya = max(lyaWList)
        
            label1 = r'$\rm 0 \leq Inc \leq 40$'
            plot1b = scatter(xList_1b,yList_1b,c='blue',marker=bSymbol,alpha=alpha,s=50)
            plot1r = scatter(xList_1r,yList_1r,c='red',marker=rSymbol,alpha=alpha,s=50)
        
            # these are matplotlib.patch.Patch properties
            props = dict(boxstyle='round', alpha=1, facecolor='none')

            # place a text box in upper right in axes coords
            ax1.text(0.59, 0.9, label1, transform=ax1.transAxes, fontsize=14,verticalalignment='top', bbox=props)

            # x-axis
            majorLocator   = MultipleLocator(0.5)
            majorFormatter = FormatStrFormatter(r'$\rm %0.1f$')
            minorLocator   = MultipleLocator(0.25)
            ax1.yaxis.set_major_locator(majorLocator)
            ax1.yaxis.set_major_formatter(majorFormatter)
            ax1.yaxis.set_minor_locator(minorLocator)
        
            # y axis
            majorLocator   = MultipleLocator(0.5)
            majorFormatter = FormatStrFormatter(r'$\rm %0.1f$')
            minorLocator   = MultipleLocator(0.25)
            ax1.yaxis.set_major_locator(majorLocator)
            ax1.yaxis.set_major_formatter(majorFormatter)
            ax1.yaxis.set_minor_locator(minorLocator)
        
            xlim(0,2.0)
            ylim(0,2.5)
            plt.setp(ax1.get_xticklabels(), visible=False)
    #         ax1.get_yaxis().set_tick_params(which='both', direction='out')

    ##########################################################################################
            ax2 = fig.add_subplot(312,sharex=ax1)
            # medium inclination = ellipse
            e = Ellipse(xy=(0,0), width=1.0, height=0.6, angle=0)
                        
            # no transparency
            e.set_alpha(0.3)
            ax2.add_artist(e)
            e.set_facecolor('black')
            e.set_edgecolor('black')
        
            label2 = r'$\rm 40 \le Inc \leq 65$'
            plot2b = scatter(xList_2b,yList_2b,c='blue',marker=bSymbol,alpha=alpha,s=50)
            plot2r = scatter(xList_2r,yList_2r,c='red',marker=rSymbol,alpha=alpha,s=50)
        
            # these are matplotlib.patch.Patch properties
            props = dict(boxstyle='round', alpha=1, facecolor='none')

            # place a text box in upper right in axes coords
            ax2.text(0.55, 0.9, label2, transform=ax2.transAxes, fontsize=14,verticalalignment='top', bbox=props)
        
            # x-axis
            majorLocator   = MultipleLocator(0.5)
            majorFormatter = FormatStrFormatter(r'$\rm %0.1f$')
            minorLocator   = MultipleLocator(0.25)
            ax2.yaxis.set_major_locator(majorLocator)
            ax2.yaxis.set_major_formatter(majorFormatter)
            ax2.yaxis.set_minor_locator(minorLocator)
        
            # y axis
            majorLocator   = MultipleLocator(0.5)
            majorFormatter = FormatStrFormatter(r'$\rm %0.1f$')
            minorLocator   = MultipleLocator(0.25)
            ax2.yaxis.set_major_locator(majorLocator)
            ax2.yaxis.set_major_formatter(majorFormatter)
            ax2.yaxis.set_minor_locator(minorLocator)
            ax2.set_yticks([0.0,0.5,1.0,1.5,2.0])


            xlim(0,2.0)
            ylim(0,2.5)
            plt.setp(ax2.get_xticklabels(), visible=False)
            ylabel(r'$\rm \rho / R_{vir}$')

    #         ax2.yaxis.tick_right()
    #         ax2.get_yaxis().set_tick_params(which='both', direction='out')

    ##########################################################################################
            ax3 = fig.add_subplot(313,sharex=ax1)
            # plot the flat galaxy line
            e = Ellipse(xy=(0,0), width=1.0, height=0.3, angle=0)
                        
            # no transparency
            e.set_alpha(0.3)
            ax3.add_artist(e)
            e.set_facecolor('black')
            e.set_edgecolor('black')
        
            label3 = r'$\rm Inc > 65$'
            plot3b = scatter(xList_3b,yList_3b,c='blue',marker=bSymbol,alpha=alpha,s=50)
            plot3r = scatter(xList_3r,yList_3r,c='red',marker=rSymbol,alpha=alpha,s=50)
        
            # these are matplotlib.patch.Patch properties
            props = dict(boxstyle='round', alpha=1, facecolor='none')

            # place a text box in upper right in axes coords
            ax3.text(0.695, 0.9, label3, transform=ax3.transAxes, fontsize=14,verticalalignment='top', bbox=props)
        
            # x-axis
            majorLocator   = MultipleLocator(0.5)
            majorFormatter = FormatStrFormatter(r'$\rm %0.1f$')
            minorLocator   = MultipleLocator(0.25)
            ax3.yaxis.set_major_locator(majorLocator)
            ax3.yaxis.set_major_formatter(majorFormatter)
            ax3.yaxis.set_minor_locator(minorLocator)
        
            # y axis
            majorLocator   = MultipleLocator(0.5)
            majorFormatter = FormatStrFormatter(r'$\rm %0.1f$')
            minorLocator   = MultipleLocator(0.25)
            ax3.yaxis.set_major_locator(majorLocator)
            ax3.yaxis.set_major_formatter(majorFormatter)
            ax3.yaxis.set_minor_locator(minorLocator)
            ax3.set_yticks([0.0,0.5,1.0,1.5,2.0])


            xlim(0,2.0)
            ylim(0,2.5)
            xlabel(r'$\rm \rho / R_{vir}$')
    #         ax3.get_yaxis().set_tick_params(which='both', direction='out')

            if save:
                savefig('{0}/azimuthMap.pdf'.format(saveDirectory),\
                format='pdf',bbox_inches='tight')
            else:
                show()


##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
        # horizontal plot
        
        else:
            # plot the distributions 
            fig = figure(figsize=(14.5,5.5))
            alpha = 0.7
            bSymbol = 'D'
            rSymbol = 'o'
        
            ax1 = fig.add_subplot(1,3,1)
#             f, (ax1, ax2, ax3) = plt.subplots(1, 3, sharey=True)
            subplots_adjust(left=None,bottom=None,right=None,top=None,wspace=0.001,hspace=None)

            # Low inclinations = circular mostly
            e = Ellipse(xy=(0,0), width=1.0, height=0.9, angle=0)
    #         legend(scatterpoints=1,prop={'size':15},loc=2,fancybox=True)
                        
            # no transparency
            e.set_alpha(0.3)
            ax1.add_artist(e)
            e.set_facecolor('black')
            e.set_edgecolor('black')
        
            maxW = 120
            minW = 5
            minLya = min(lyaWList)
            maxLya = max(lyaWList)
        
            label1 = r'$\rm Inc \leq 40$'
            ax1.scatter(xList_1b,yList_1b,c='blue',marker=bSymbol,alpha=alpha,s=50)
            ax1.scatter(xList_1r,yList_1r,c='red',marker=rSymbol,alpha=alpha,s=50)
        
            # these are matplotlib.patch.Patch properties
            props = dict(boxstyle='round', alpha=1, facecolor='none')

            # place a text box in upper right in axes coords
            ax1.text(0.715, 0.93, label1, transform=ax1.transAxes, fontsize=15,verticalalignment='top', bbox=props)

            # x-axis
            ax1.set_xticks([0.0,0.5,1.0,1.5])

            majorLocator   = MultipleLocator(0.5)
            majorFormatter = FormatStrFormatter(r'$\rm %0.1f$')
            minorLocator   = MultipleLocator(0.25)
            ax1.yaxis.set_major_locator(majorLocator)
            ax1.yaxis.set_major_formatter(majorFormatter)
            ax1.yaxis.set_minor_locator(minorLocator)

        
            # y axis
            majorLocator   = MultipleLocator(0.5)
            majorFormatter = FormatStrFormatter(r'$\rm %0.1f$')
            minorLocator   = MultipleLocator(0.25)
            ax1.yaxis.set_major_locator(majorLocator)
            ax1.yaxis.set_major_formatter(majorFormatter)
            ax1.yaxis.set_minor_locator(minorLocator)
        
            xlim(0,2.0)
            ylim(0,2.5)
            ylabel(r'$\rm \rho / R_{vir}$')

#             plt.setp(ax1.get_xticklabels(), visible=False)
    #         ax1.get_yaxis().set_tick_params(which='both', direction='out')

    ##########################################################################################
            ax2 = fig.add_subplot(1,3,2,sharey=ax1)
            # medium inclination = ellipse
            e = Ellipse(xy=(0,0), width=1.0, height=0.6, angle=0)
                        
            # no transparency
            e.set_alpha(0.3)
            ax2.add_artist(e)
            e.set_facecolor('black')
            e.set_edgecolor('black')
        
            label2 = r'$\rm 40 < Inc \leq 65$'
            ax2.scatter(xList_2b,yList_2b,c='blue',marker=bSymbol,alpha=alpha,s=50)
            ax2.scatter(xList_2r,yList_2r,c='red',marker=rSymbol,alpha=alpha,s=50)
        
            # these are matplotlib.patch.Patch properties
            props = dict(boxstyle='round', alpha=1, facecolor='none')

            # place a text box in upper right in axes coords
            ax2.text(0.58, 0.93, label2, transform=ax2.transAxes, fontsize=15,verticalalignment='top', bbox=props)
        
            # x-axis
            majorLocator   = MultipleLocator(0.5)
            majorFormatter = FormatStrFormatter(r'$\rm %0.1f$')
            minorLocator   = MultipleLocator(0.25)
            ax2.yaxis.set_major_locator(majorLocator)
            ax2.yaxis.set_major_formatter(majorFormatter)
            ax2.yaxis.set_minor_locator(minorLocator)
            ax2.set_xticks([0.0,0.5,1.0,1.5])

            # y axis
            majorLocator   = MultipleLocator(0.5)
            majorFormatter = FormatStrFormatter(r'$\rm %0.1f$')
            minorLocator   = MultipleLocator(0.25)
            ax2.yaxis.set_major_locator(majorLocator)
            ax2.yaxis.set_major_formatter(majorFormatter)
            ax2.yaxis.set_minor_locator(minorLocator)

            xlim(0,2.0)
            ylim(0,2.5)
            plt.setp(ax2.get_yticklabels(), visible=False)
            xlabel(r'$\rm \rho / R_{vir}$')

    #         ax2.yaxis.tick_right()
    #         ax2.get_yaxis().set_tick_params(which='both', direction='out')

    ##########################################################################################
            ax3 = fig.add_subplot(1,3,3,sharey=ax1)
            # plot the flat galaxy line
            e = Ellipse(xy=(0,0), width=1.0, height=0.3, angle=0)
                        
            # no transparency
            e.set_alpha(0.3)
            ax3.add_artist(e)
            e.set_facecolor('black')
            e.set_edgecolor('black')
        
            label3 = r'$\rm Inc > 65$'
            ax3.scatter(xList_3b,yList_3b,c='blue',marker=bSymbol,alpha=alpha,s=50)
            ax3.scatter(xList_3r,yList_3r,c='red',marker=rSymbol,alpha=alpha,s=50)
        
            # these are matplotlib.patch.Patch properties
            props = dict(boxstyle='round', alpha=1, facecolor='none')

            # place a text box in upper right in axes coords
            ax3.text(0.715, 0.93, label3, transform=ax3.transAxes, fontsize=15,verticalalignment='top', bbox=props)
        
            # x-axis
            majorLocator   = MultipleLocator(0.5)
            majorFormatter = FormatStrFormatter(r'$\rm %0.1f$')
            minorLocator   = MultipleLocator(0.25)
            ax3.yaxis.set_major_locator(majorLocator)
            ax3.yaxis.set_major_formatter(majorFormatter)
            ax3.yaxis.set_minor_locator(minorLocator)
        
            # y axis
            majorLocator   = MultipleLocator(0.5)
            majorFormatter = FormatStrFormatter(r'$\rm %0.1f$')
            minorLocator   = MultipleLocator(0.25)
            ax3.yaxis.set_major_locator(majorLocator)
            ax3.yaxis.set_major_formatter(majorFormatter)
            ax3.yaxis.set_minor_locator(minorLocator)
#             ax3.set_yticks([0.0,0.5,1.0,1.5,2.0])

            xlim(0,2.0)
            ylim(0,2.5)
#             xlabel(r'$\rm \rho / R_{vir}$')
            plt.setp(ax3.get_yticklabels(), visible=False)

    #         ax3.get_yaxis().set_tick_params(which='both', direction='out')

            if save:
                savefig('{0}/azimuthMap_horizontal.pdf'.format(saveDirectory),\
                format='pdf',bbox_inches='tight')
            else:
                show()



###############################################################################
###############################################################################
###############################################################################
###############################################################################


if __name__=="__main__":
    # do the work
    main()
    