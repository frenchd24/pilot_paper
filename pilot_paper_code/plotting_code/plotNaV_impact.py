#!/usr/bin/env python

'''
By David French (frenchd@astro.wisc.edu)

$Id:  plotNaV_impact.py, v 5.0 01/04/2016

This is the plotNaV_b_diam bit from histograms3.py. Now is separated, and loads in a pickle
file of the relevant data, as created by "buildDataLists.py"

plots apparent column density vs impact parameter, and also vs impact / diameter and 
impact parameter / R_vir

Previous (from histograms3.py):
    Plot some stuff for the 100largest initial results

    Make plots for AAS winter 2014 poster. Uses LG_correlation_combined2.csv file

    Updated for the pilot paper (05/06/15)


v5: updated to work with the new, automatically updated LG_correlation_combined5.csv
    (12/04/15) - original updates to the individual files
    
    - added second function to normalize impact parameter by R_vir instead of maj axis
        (1/3/16)
        
    - combined plotNaV_b_diam2.py and plotNaV_b2.py into this single file that contains
        all the previous functions (they were all related)
        (1/3/16)

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
        resultsFilename = '/Users/David/Research_Documents/inclination/git_inclination/LG_correlation_combined5_3.csv'
        saveDirectory = '/Users/David/Research_Documents/inclination/git_inclination/pilot_paper_code/plots/'

    elif getpass.getuser() == 'frenchd':
        pickleFilename = '/usr/users/frenchd/inclination/git_inclination/pilot_paper_code/pilotData2.p'
        resultsFilename = '/usr/users/frenchd/inclination/git_inclination/LG_correlation_combined5_3.csv'
        saveDirectory = '/usr/users/frenchd/inclination/git_inclination/pilot_paper_code/plots/'

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
                
            if isNumber(maj):
                virialRadius = float(virialRadius)
                maj = float(maj)
            else:
                virialRadius = -99
                maj = -99
            
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
    
    
    
########################################################################################
########################################################################################

    # plot apparent column density as a function of impact parameter, split between red and
    # blue shifted absorption
    #
    
    plotNaV_b = True
    save = False
    
    if plotNaV_b:
        fig = figure()
        ax = fig.add_subplot(111)
        countb = 0
        countr = 0
        count = -1
        labelr = 'Red Shifted Absorber'
        labelb = "Blue Shifted Absorber"
        for d,i,n in zip(difList,impactList,naList):
        
            # make sure all the values are good
            if isNumber(d) and isNumber(i) and isNumber(n):
                if d !=-99 and i !=-99 and n!=-99:
                    if d>0:
                        # galaxy is behind absorber, so gas is blue shifted
                        color = 'Blue'
                        if countb == 0:
                            countb +=1
                            plotb = ax.scatter(i,n,c='Blue',s=50,label= labelb)
                    if d<0:
                        # gas is red shifted compared to galaxy
                        color = 'Red'
                        if countr == 0:
                            countr +=1
                            plotr = ax.scatter(i,n,c='Red',s=50,label= labelr)
                
                    plot1 = scatter(i,n,c=color,s=50)
            
        title('Apparent N(HI) vs impact parameter for red vs blue absorption')
        xlabel('Impact Parameter (kpc)')
        ylabel(r'Na(v) (cm$^{\rm -2}$)')
        legend(scatterpoints=1)
        ax.grid(b=None,which='major',axis='both')
        ylim(-1,5e14)
        xlim(0,500)
        
        if save:
            savefig('{0}/NaV(impact)_dif.pdf'.format(saveDirectory),format='pdf')
        else:
            show()
        

########################################################################################
########################################################################################        

    # plot apparent column density as a function of impact parameter/diameter for red
    # and blue shifted absorption
    #
    
    plotNaV_b_diam = False
    save = False
    
    if plotNaV_b_diam:
        fig = figure()
        ax = fig.add_subplot(111)
        countb = 0
        countr = 0
        count = -1
        labelr = 'Red Shifted Absorber'
        labelb = "Blue Shifted Absorber"
        for d,i,n,m in zip(difList,impactList,naList,majList):
        
            # check if all the values are good
            if isNumber(d) and isNumber(i) and isNumber(n) and isNumber(m):
                if d !=-99 and i !=-99 and n!=-99 and m!=-99:
                    if d>0:
                        # galaxy is behind absorber, so gas is blue shifted
                        color = 'Blue'
                        if countb == 0:
                            countb +=1
                            plotb = ax.scatter(i/m,n,c='Blue',s=50,label= labelb)
                    if d<0:
                        # gas is red shifted compared to galaxy
                        color = 'Red'
                        if countr == 0:
                            countr +=1
                            plotr = ax.scatter(i/m,n,c='Red',s=50,label= labelr)
                
                    plot1 = scatter(i/m,n,c=color,s=50)
            
        title('Apparent N(HI) vs impact/diameter for red vs blue absorption')
        xlabel(r'Impact Parameter / Diameter')
        ylabel(r'Na(v) (cm$^{\rm -2}$)')
        legend(scatterpoints=1)
        ax.grid(b=None,which='major',axis='both')
        ylim(-1,5e14)
        xlim(-1,70)
        
        if save:
            savefig('{0}/NaV(impact_diameter)_dif.pdf'.format(saveDirectory),format='pdf')
        else:
            show()
                 

########################################################################################
########################################################################################        

    # plot apparent column density as a function of impact parameter/R_virial for red
    # and blue shifted absorption
    #
    
    plotNaV_b_vir = False
    save = False
    
    if plotNaV_b_vir:
        fig = figure()
        ax = fig.add_subplot(111)
        countb = 0
        countr = 0
        count = -1
        labelr = 'Red Shifted Absorber'
        labelb = "Blue Shifted Absorber"
        for d,i,n,v in zip(difList,impactList,naList,virList):
        
            # check if all the values are good
            if isNumber(d) and isNumber(i) and isNumber(n) and isNumber(v):
                if d !=-99 and i !=-99 and n!=-99 and v!=-99:
                    if d>0:
                        # galaxy is behind absorber, so gas is blue shifted
                        color = 'Blue'
                        if countb == 0:
                            countb +=1
                            plotb = ax.scatter(i/v,n,c='Blue',s=50,label= labelb)
                    if d<0:
                        # gas is red shifted compared to galaxy
                        color = 'Red'
                        if countr == 0:
                            countr +=1
                            plotr = ax.scatter(i/v,n,c='Red',s=50,label= labelr)
                
                    plot1 = scatter(i/v,n,c=color,s=50)
            
        title('Apparent N(HI) vs impact/R_vir for red vs blue absorption')
        xlabel(r'Impact Parameter / R_vir')
        ylabel(r'Na(v) (cm$^{\rm -2}$)')
        legend(scatterpoints=1)
        ax.grid(b=None,which='major',axis='both')
#         ylim(-1,5e14)
#         xlim(-1,70)
        
        if save:
            savefig('{0}/NaV(impact_vir)_dif.pdf'.format(saveDirectory),format='pdf')
        else:
            show()
                 
    
###############################################################################
###############################################################################
###############################################################################
###############################################################################


if __name__=="__main__":
    # do the work
    main()
    