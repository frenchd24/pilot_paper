#!/usr/bin/env python

'''
By David French (frenchd@astro.wisc.edu)

$Id:  plotVirHist.py, v 1.1 4/21/16

Plots histograms of virial radii for associated galaxies

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
        pickleFilename = '/Users/David/Research_Documents/inclination/git_inclination/pilot_paper_code/pilotData2.p'
#         resultsFilename = '/Users/David/Research_Documents/inclination/git_inclination/LG_correlation_combined5_8_edit2.csv'
#         saveDirectory = '/Users/David/Research_Documents/inclination/git_inclination/pilot_paper_code/plots2/'
        resultsFilename = '/Users/David/Research_Documents/inclination/git_inclination/LG_correlation_combined5_9_edit2.csv'
        saveDirectory = '/Users/David/Research_Documents/inclination/git_inclination/pilot_paper_code/plots3/'

    elif getpass.getuser() == 'frenchd':
        pickleFilename = '/usr/users/frenchd/inclination/git_inclination/pilot_paper_code/pilotData2.p'
#         resultsFilename = '/usr/users/frenchd/inclination/git_inclination/LG_correlation_combined5_8_edit2.csv'
#         saveDirectory = '/usr/users/frenchd/inclination/git_inclination/pilot_paper_code/plots2/'
        resultsFilename = '/usr/users/frenchd/inclination/git_inclination/LG_correlation_combined5_9_edit2.csv'
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
            
            if isNumber(virialRadius):
                virialRadius = float(virialRadius)
            else:
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
    # make a histogram of the distribution of R_vir for associated galaxies 
    #
    
    plotR_vir= True
    save = True
    
    if plotR_vir:
        fig = figure(figsize=(10,2))
#         subplots_adjust(hspace=0.200)
        ax = fig.add_subplot(111)
        
#         bins = [0,10,20,30,40,50,60,70,80,90]
        print virList
        
        binsize = 25
        bins = arange(min(virList),max(virList)+binsize,binsize)   
        
        impactArray = np.array([])
        virArray = np.array([])
        
        # make sure all the values are okay
        for i,v in zip(impactList,virList):
            if isNumber(i) and isNumber(v):
                if i !=-99 and v !=-99:
                    impactArray = append(impactArray,float(i))
                    virArray = append(virArray,float(v))
        
        
        plot1 = hist(virArray,bins=bins,histtype='bar')
        
        xlabel(r'$\rm \rho$ (kpc)')
        ylabel('Number')
        ax.tick_params(axis='x', labelsize=0)
        ax.tick_params(axis='y',labelsize=8)
#         ylim(0,5)

        if save:
            savefig('{0}/hist(vir).pdf'.format(saveDirectory),format='pdf')
        else:
            show()


##########################################################################################
##########################################################################################
    # make a histogram of the distribution of R_vi for associated galaxies, split
    # into red and blue shifted bins
    #
    
    plotVirHist_dif = True
    save = True
    
    if plotVirHist_dif:
        fig = figure(figsize=(10,5))
        ax = fig.add_subplot(211)
        subplots_adjust(hspace=0.200)
        alpha = 0.85

        binsize = 25
        bins = arange(min(virList),max(virList)+binsize,binsize)
        
        reds = np.array([])
        blues = np.array([])
        
        # make sure all the values are okay
        for i,v,d in zip(impactList,virList,difList):
            if isNumber(i) and isNumber(v):
                if i !=-99 and v !=-99:
                    val = float(v)
                    
                    # if blueshifted
                    if d>=0:
                        blues = append(blues,val)
                    
                    # for redshifted
                    else:
                        reds = append(reds,val)
        
#         title('Distribution of impact parameters')

        plot1 = hist(blues,bins=bins,histtype='bar',color='blue',label='Blueshifted',alpha=alpha)
#         xlabel('Impact Parameter/R_vir (blueshifted)')
        ylabel('Number')
        legend(scatterpoints=1,prop={'size':12})
        ax.tick_params(axis='y',labelsize=11)
        ax.tick_params(axis='x', labelsize=0)

        ax = fig.add_subplot(212)
        plot2 = hist(reds,bins=bins,histtype='bar',color="red",label='Redshifted',alpha=alpha)
        xlabel(r'$\rm R_{vir}$')
        ylabel('Number')
        
        legend(scatterpoints=1,prop={'size':12})
        
        ax.tick_params(axis='x', labelsize=11)
        ax.tick_params(axis='y',labelsize=11)
#         ylim(0,5)

        if save:
            savefig('{0}/hist(vir)_dif_bin_{1}.pdf'.format(saveDirectory,binsize),format='pdf')
        else:
            show()


##########################################################################################
##########################################################################################
    # plot the velocity difference of absorbers as a function of impact/ virial radius of 
    # the associated galaxy
    #
    
    plotImpact_Vir_vs_vel_dif = False
    save = False
    
    if plotImpact_Vir_vs_vel_dif:
        fig = figure()
        ax = fig.add_subplot(111)
        countb = 0
        countr = 0
        count = -1
        labelr = 'Redshifted Absorber'
        labelb = "Blueshifted Absorber"
        alpha = 0.85

        for d,i,v in zip(difList,impactList,virList):
            # check if all the values are okay
            if isNumber(d) and isNumber(i) and isNumber(v):
                if d!=-99 and i!=-99 and v!=-99:
#                     print 'd: ',d
                    xVal = float(i)/float(v)
                    yVal = d
                    
                    if d>0:
                        # galaxy is 'behind' absorber, so GAS = blue shifted
                        color = 'Blue'
                        if countb == 0:
                            countb +=1
                            plotb = ax.scatter(xVal,yVal,c='Blue',s=50,label= labelb,alpha=alpha)
                    if d<0:
                        # gas is RED shifted compared to galaxy
                        color = 'Red'
                        if countr == 0:
                            countr +=1
                            plotr = ax.scatter(xVal,yVal,c='Red',s=50,label= labelr,alpha=alpha)
                
                    plot1 = scatter(xVal,yVal,c=color,s=50,alpha=alpha)
            
        xlabel(r'$\rm \rho/R_{vir}$ (kpc)')
        ylabel(r'$\rm \Delta v$ (km/s)')
        legend(scatterpoints=1,prop={'size':12},loc=1)
        ax.grid(b=None,which='major',axis='both')
        xlim(0,3)
        ylim(-400,400)

        if save:
            savefig('{0}/vel_dif(impact_vir).pdf'.format(saveDirectory),format='pdf')
        else:
            show()
            

##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################


if __name__=="__main__":
    # do the work
    main()
    