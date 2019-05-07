#!/usr/bin/env python

'''
By David French (frenchd@astro.wisc.edu)

$Id:  plotLyaWHist_both2.py, v 5.2 04/21/2016

DEPRECIATED. SEE 'plotWHist.py' FOR CURRENT VERSION


This is the plotLyaWHist_both bit from histograms3.py. Now is separated, and loads in a pickle
file of the relevant data, as created by "buildDataLists.py"

Previous (from histograms3.py):
    Plot some stuff for the 100largest initial results

    Make plots for AAS winter 2014 poster. Uses LG_correlation_combined2.csv file

    Updated for the pilot paper (05/06/15)


v5: updated to work with the new, automatically updated LG_correlation_combined5.csv
    (12/04/15) - original updates to the individual files
    
v5.1: update for LG_correlation_combined5_8_edit2.csv and l_min = 0.001 (02/23/2016)

v5.2: remake plots with v_hel instead of vcorr (4/21/16)


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

fontScale = 16
rc('text', usetex=True)
rc('font',size=18)
rc('xtick.major',size=5,width=1.2)
rc('xtick.minor',size=3,width=1.2)
rc('ytick.major',size=5,width=1.2)
rc('ytick.minor',size=3,width=1.2)
rc('xtick',labelsize=16)
rc('ytick',labelsize=16)
rc('axes',labelsize=16)
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
#         resultsFilename = '/Users/David/Research_Documents/inclination/git_inclination/LG_correlation_combined5_9_edit2.csv'
        resultsFilename = '/Users/David/Research_Documents/inclination/git_inclination/LG_correlation_combined5_9_corrected.csv'
        saveDirectory = '/Users/David/Research_Documents/inclination/git_inclination/pilot_paper_code/plots3/'

    elif getpass.getuser() == 'frenchd':
        pickleFilename = '/usr/users/frenchd/inclination/git_inclination/pilot_paper_code/pilotData2.p'
#         resultsFilename = '/usr/users/frenchd/inclination/git_inclination/LG_correlation_combined5_8_edit2.csv'
#         saveDirectory = '/usr/users/frenchd/inclination/git_inclination/pilot_paper_code/plots2/'
#         resultsFilename = '/usr/users/frenchd/inclination/git_inclination/LG_correlation_combined5_9_edit2.csv'
        resultsFilename = '/usr/users/frenchd/inclination/git_inclination/LG_correlation_combined5_9_corrected.csv'
        saveDirectory = '/usr/users/frenchd/inclination/git_inclination/pilot_paper_code/plots3/'

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
            min = l['minorAxis (kpc)']
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
    # make a histogram of the distribution of Ly-alpha equivalent widths for both the 
    # associated and ambiguous samples
    #
    #
    # normByEnv doesn't work for the ambigous lines, because most of them have env = 0
    #
    
    plotLyaWHist_both = False
    save = False
    
    if plotLyaWHist_both:
#         fig = figure(figsize=(2,8))
        fig = figure()
        ax = fig.add_subplot(211)
#         bins = [0,.10,.20,.30,.40,.50,.60,.70,.80,.90]
#         bins = arange(0,max(max(lyaWList),max(lyaWAmbList)),20)
        bins = 15
        
        lyaWArray = array(lyaWList)
        lyaWAmbArray = array(lyaWAmbList)
        
        envArray = array(envList)
        envAmbArray = array(envAmbList)
    
        print 'lyaWAmbArray: ',lyaWAmbArray
        print 'envAmbArray: ',envAmbArray
    
        # see above, does not work the for the ambigous ones
        normByEnv = False
        
        if normByEnv:
            plot1 = hist(lyaWArray/envArray,bins=bins,histtype='bar',orientation = 'vertical')
            title('Distribution of Lya W - Associated')

            ax = fig.add_subplot(212)
            plot1 = hist(lyaWAmbArray/envAmbArray,bins=bins,histtype='bar',orientation = 'vertical')

        else:
            plot1 = hist(lyaWList,bins=bins,histtype='bar',orientation = 'vertical',label='Associated')
            legend(scatterpoints=1,prop={'size':12},loc=1)
            ylabel('Number')

            ax = fig.add_subplot(212)
            plot1 = hist(lyaWAmbList,bins=bins,histtype='bar',orientation = 'vertical',label='Ambiguous')
            legend(scatterpoints=1,prop={'size':12},loc=1)

        
        xlabel(r'Equivalent Width ($\rm m\AA$)')
        ylabel('Number')
        ax.tick_params(axis='x', labelsize=11)
        ax.tick_params(axis='y',labelsize=11)
#         xlim(0,11)
#         tight_layout()
        
        if save:
            savefig('{0}/hist(lyaW_assoc_vs_ambig).pdf'.format(saveDirectory),format='pdf')
        else:
            show()


#########################################################################################
#########################################################################################
    # make a histogram of the distribution of Ly-alpha equivalent widths for both the 
    # associated sample, splitting on red vs blue shifted absorbers
    #
    #
    # normByEnv doesn't work for the ambigous lines, because most of them have env = 0
    #
    
    plotLyaWHist_dif = True
    save = False
    
    if plotLyaWHist_dif:
#         fig = figure(figsize=(2,8))
        fig = figure()
        ax = fig.add_subplot(211)
#         bins = [0,.10,.20,.30,.40,.50,.60,.70,.80,.90]
#         bins = arange(0,max(max(lyaWList),max(lyaWAmbList)),20)
        bins = arange(0,1300,100)
        
        lyaWArray = array(lyaWList)
        lyaWAmbArray = array(lyaWAmbList)
        
        envArray = array(envList)
        envAmbArray = array(envAmbList)
    
        print 'lyaWAmbArray: ',lyaWAmbArray
        print 'envAmbArray: ',envAmbArray
        
        blues = []
        reds = []
        for d,l in zip(difList,lyaWList):
            if float(d) >0:
                # blueshifted
                blues.append(float(l))
            else:
                reds.append(float(l))
    
        # see above, does not work the for the ambigous ones
        normByEnv = False
        
        
        majorLocator   = MultipleLocator(10)
        majorFormatter = FormatStrFormatter('%d')
        minorLocator   = MultipleLocator(2)

        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        
        majorLocator   = MultipleLocator(2)
        majorFormatter = FormatStrFormatter('%d')
        minorLocator   = MultipleLocator(1)
        
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)
        
        
        if normByEnv:
            plot1 = hist(lyaWArray/envArray,bins=bins,histtype='bar',orientation = 'vertical')
            title('Distribution of Lya W - Associated')

            ax = fig.add_subplot(212)
            plot1 = hist(lyaWAmbArray/envAmbArray,bins=bins,histtype='bar',orientation = 'vertical')

        else:
            plot1 = hist(blues,bins=bins,histtype='bar',orientation = 'vertical',label='Blueshifted',color="Blue",alpha = 0.85)
            legend(scatterpoints=1,prop={'size':12},loc=1)
            ylabel(r'Number')
            ylim(0,10)

            ax = fig.add_subplot(212)
            plot1 = hist(reds,bins=bins,histtype='bar',orientation = 'vertical',label='Redshifted',color="Red",alpha = 0.85)
            legend(scatterpoints=1,prop={'size':12},loc=1)
            ylabel(r'Number')
            ylim(0,10)

        
        xlabel(r'Equivalent Width ($\rm m\AA$)')
        ax.tick_params(axis='x', labelsize=11)
        ax.tick_params(axis='y',labelsize=11)
        xlim(0,1200)
        
        if save:
            savefig('{0}/hist(lyaW_blue_vs_red).pdf'.format(saveDirectory),format='pdf')
        else:
            show()
            



#########################################################################################
#########################################################################################
    # make a histogram of the distribution of Ly-alpha equivalent widths for both the 
    # associated sample, splitting on red vs blue shifted absorbers
    #
    #
    # normByEnv doesn't work for the ambigous lines, because most of them have env = 0
    #
    
    plotLyaWHist_dif = True
    save = False
    
    if plotLyaWHist_dif:
#         fig = figure(figsize=(2,8))
        fig = figure()
        ax = fig.add_subplot(211)
#         bins = [0,.10,.20,.30,.40,.50,.60,.70,.80,.90]
#         bins = arange(0,max(max(lyaWList),max(lyaWAmbList)),20)
        bins = arange(0,1300,100)
        
        lyaWArray = array(lyaWList)
        lyaWAmbArray = array(lyaWAmbList)
        
        envArray = array(envList)
        envAmbArray = array(envAmbList)
        
        blues = []
        reds = []
        for d,l in zip(difList,lyaWList):
            if float(d) >0:
                # blueshifted
                blues.append(float(l))
            else:
                reds.append(float(l))

        sorted_data = sort(randImpact_nolike)
    
        # compute the CDF y-values
        yvals=np.arange(len(sorted_data)) / float(len(sorted_data))
        
        # format all the axis and stuff
        majorLocator   = MultipleLocator(10)
        majorFormatter = FormatStrFormatter('%d')
        minorLocator   = MultipleLocator(2)

        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        
        majorLocator   = MultipleLocator(2)
        majorFormatter = FormatStrFormatter('%d')
        minorLocator   = MultipleLocator(1)
        
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)
        
        plot1 = hist(blues,bins=bins,histtype='bar',orientation = 'vertical',label='Blueshifted',color="Blue",alpha = 0.85)
        legend(scatterpoints=1,prop={'size':12},loc=1)
        ylabel(r'Number')
        ylim(0,10)

        ax = fig.add_subplot(212)
        plot1 = hist(reds,bins=bins,histtype='bar',orientation = 'vertical',label='Redshifted',color="Red",alpha = 0.85)
        legend(scatterpoints=1,prop={'size':12},loc=1)
        ylabel(r'Number')
        ylim(0,10)
    
        # plot the y-values against the sorted data
        plt.plot(sorted_data,yvals)
                
        xlabel(r'$\rm Equivalent Width (m\AA)$')
        ylabel(r'$\rm Number$')

        tight_layout()
#         ylim(0,5)

        if save:
            savefig('{0}/CDF(lyaW_blue_red_amb).pdf'.format(saveDirectory),format='pdf')
        else:
            show()
    
###############################################################################
###############################################################################
###############################################################################
###############################################################################


if __name__=="__main__":
    # do the work
    main()
    