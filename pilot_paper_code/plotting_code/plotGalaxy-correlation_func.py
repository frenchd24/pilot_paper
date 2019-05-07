#!/usr/bin/env python

'''
By David French (frenchd@astro.wisc.edu)

$Id:  plotGalaxy-correlation_func.py, v 1.0 04/25/2016


Plot \rho/R_vir for a random selection of galaxies, and positions on the sky (4/25/16)

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
import random

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
# rc('text', usetex=True)
# rc('font',size=16,weight='bold')
    
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


##########################################################################################
def main():
    # assuming 'theFile' contains one name per line, read the file
    
    if getpass.getuser() == 'David':
        pickleFilename = '/Users/David/Research_Documents/inclination/git_inclination/pilot_paper_code/pilotData2.p'
        saveDirectory = '/Users/David/Research_Documents/inclination/git_inclination/pilot_paper_code/plots3/'
        galaxyFilename = '/Users/David/Research_Documents/gt/NewGalaxyTable5.csv'
        pickle_outFilename = '/Users/David/Research_Documents/inclination/rand_pickle2.p'

    elif getpass.getuser() == 'frenchd':
        pickleFilename = '/usr/users/frenchd/inclination/git_inclination/pilot_paper_code/pilotData2.p'
        saveDirectory = '/usr/users/frenchd/inclination/git_inclination/pilot_paper_code/plots3/'
        galaxyFilename = '/usr/users/frenchd/gt/NewGalaxyTable5.csv'
        pickle_outFilename = '/usr/users/frenchd/inclination/rand_pickle2.p'


    else:
        print 'Could not determine username. Exiting.'
        sys.exit()
    
    # use the old pickle file to get the full galaxy dataset info
    pickleFile = open(pickleFilename,'rU')
    fullDict = pickle.load(pickleFile)
    pickleFile.close()
    
    # open a new pickle file
    pickle_outFile = open(pickle_outFilename,'wt')
    
    galFile = open(galaxyFilename,'rU')
    reader = csv.DictReader(galFile)
    
    galRA = []
    galDec = []
    galVel = []
    galR = []
    galPA = []
    galInc = []
    galDist = []
    
    # grab galaxy locations, velocities, and sizes
    for l in reader:
        ra,dec = eval(l["degreesJ2000RA_Dec"])
        vel = l['radialVelocity (km/s)']
        major, minor = eval(l['linDiameters (kpc)'])
        pa = l['positionAngle (deg)']
        pa3 = l['RC3pa (deg)']
        inc = l['inclination (deg)']
        inc3 = l['RC3inc (deg)']
        dist = l['Best Distance (Mpc)']
        
        if isNumber(pa3):
            if not isNumber(pa):
                pa = pa3
        else:
            pa = -99
            
        if isNumber(inc3):
            if not isNumber(inc):
                inc = inc3
        else:
            inc = -99
        
        if isNumber(vel) and isNumber(major) and isNumber(dist):
            vel = float(vel)
            if vel >= 400 and float(major) >0:
                galRA.append(float(ra))
                galDec.append(float(dec))
                galVel.append(vel)
                
                R_vir = calculateVirialRadius(float(major))
                
                galR.append(R_vir)
                galInc.append(float(inc))
                galPA.append(float(pa))
                galDist.append(float(dist))
    
    # close the file
    galFile.close()
    
    # now randomly generate some locations in the sky
    
#     random.randrange(start, stop[, step])   
    num = 1000
    count = 0
    randLower = 400.
    randUpper = 10000.
    step = 1
    randVel = []
    randRA = []
    randDec = []
    randLike = []
    randRho = []
    randImpact_rho = []
    randAz = []
    randInc = []
    
    # ignore likelihood, just include the closest galaxy in all cases
    randImpact_nolike = []
    randImpact_rho_nolike = []
    randDelV_nolike = []
    
    
    while count <=num:
        count +=1
        sys.stdout.write('Finished {0}\r'.format(count))
        sys.stdout.flush()
        
        vel = random.randint(randLower,randUpper)
        ra = random.random()*360.
        s = random.random()      
          
        if s<=0.5:
            dec = -random.random()*90.
        else:
            dec = random.random()*90.
            
        possible = []
        
        sImp = 999999
        sImp_rho = 999999
        sDelV = 999999
        for gRA,gDec,gVel,r,pa,inc,dist in zip(galRA,galDec,galVel,galR,galPA,galInc,galDist):
            
            # calculate impact parameter to this random sightline
            imp = calculateImpactParameter(gRA,gDec,ra,dec,dist)
            
            imp_r = imp/r
            
            # calculate likelihood
            delV = gVel - vel
            likelihood = math.exp(-(imp/r)**2) * math.exp(-(delV/200.)**2)
            
            # is this the smallest impact parameter yet?
            if imp < sImp:
                # if yes, set it as sImp
                sImp = imp
            
            if imp_r < sImp_rho:
                sImp_rho = imp_r
            
            if delV < sDelV:
                sDelV = delV
                
            # now likelihood
            if abs(delV) <=400:
                if r>= imp:
                    likelihood = likelihood*2
            
                if likelihood>=0.001:
                    az = -99
                    if isNumber(galPA):
                        if float(galPA) >=0:
                            az = calculateAzimuth(galRA,galDec,ra,dec,galDist,galPA)
                        
                    possible.append([likelihood,[imp,delV,r,inc,az]])
        
        # append the minimum values to the nolike lists
        randImpact_nolike.append(sImp)
        randImpact_rho_nolike.append(sImp_rho)
        randDelV_nolike.append(sDelV)
        
        possible.sort(reverse=True)
        
        # top two most likely
        if len(possible) >1:
            best = possible[0]
            second = possible[1]
        
            if best[0] >= 5* second[0]:
                rest = best[1]
                randVel.append(rest[1])
                randLike.append(best[0])
                randRho.append(rest[0])
            
                impact_rho = rest[0]/rest[2]
                randImpact_rho.append(impact_rho)
                randAz.append(rest[4])
                randInc.append(rest[3])
                                
        elif len(possible) >0:
            best = possible[0]

            rest = best[1]
            randVel.append(rest[1])
            randLike.append(best[0])
            randRho.append(rest[0])
        
            impact_rho = rest[0]/rest[2]
            randImpact_rho.append(impact_rho)
            randAz.append(rest[4])
            randInc.append(rest[3])
                        

    pickle.dump([randVel,randLike,randRho,randImpact_rho,randAz,randInc, randImpact_nolike,\
    randImpact_rho_nolike,randDelV_nolike],pickle_outFile)
    pickle_outFile.close()

            
    # lists for the full galaxy dataset
#     allPA = fullDict['allPA']
#     allInclinations = fullDict['allInclinations']
#     allCosInclinations = fullDict['allCosInclinations']
#     allFancyInclinations = fullDict['allFancyInclinations']
#     allCosFancyInclinations = fullDict['allCosFancyInclinations']
    

##########################################################################################
##########################################################################################
    # make a histogram of the distribution of impact/R_vir for random locations
    #
    
    plotRandHist = True
    save = True
    
    if plotRandHist:
        fig = figure(figsize=(10,4))
        ax = fig.add_subplot(111)
        alpha = 0.75

        bins = arange(0,4,0.2)
        
#         impactArray = np.array([])
#         virArray = np.array([])
#         
#         # make sure all the values are okay
#         for i,v in zip(impactList,virList):
#             if isNumber(i) and isNumber(v):
#                 if i !=-99 and v !=-99:
#                     impactArray = append(impactArray,float(i))
#                     virArray = append(virArray,float(v))
                
        print 'median: ',median(randImpact_rho)
        print 'mean: ',mean(randImpact_rho)
        
        plot1 = hist(randImpact_rho,bins=bins,histtype='bar',alpha=alpha)
        
        xlabel(r'$\rm \rho / R_{vir}$')
        ylabel(r'$\rm Number$')
#         grid(True)
    
#         ax.tick_params(axis='x', labelsize=10)
#         ax.tick_params(axis='y',labelsize=10)
        tight_layout()
#         ylim(0,5)

        if save:
            savefig('{0}/hist(random_impact_vir)update_{1}.pdf'.format(saveDirectory,num),format='pdf')
        else:
            show()

##########################################################################################
##########################################################################################
    # make a histogram of the distribution of impact parameters for associated galaxies 
    # normalized by virial radius
    #
    
    plotImpactHist_Vir = False
    save = False
    
    if plotImpactHist_Vir:
        fig = figure(figsize=(10,4))
        ax = fig.add_subplot(111)
        alpha = 0.75
        
#         bins = [0,.10,.20,.30,.40,.50,.60,.70,.80,.90]
    #     bins = [5,15,25,35,45,55,65,75,85]
    #     bins = [0,15,30,45,60,75,90]
#         bins = [0,5,10,15,20,25,30,35,40]
        bins = arange(0,3.5,0.2)
        
        impactArray = np.array([])
        virArray = np.array([])
        
        # make sure all the values are okay
        for i,v in zip(impactList,virList):
            if isNumber(i) and isNumber(v):
                if i !=-99 and v !=-99:
                    impactArray = append(impactArray,float(i))
                    virArray = append(virArray,float(v))
        
        normalizedImpactArray = impactArray/virArray
        
        print 'median: ',median(normalizedImpactArray)
        print 'mean: ',mean(normalizedImpactArray)
        
        plot1 = hist(normalizedImpactArray,bins=bins,histtype='bar',alpha=alpha)
        
#         title('Distribution of impact parameters')
        xlabel(r'$\rho / R_{vir}$')
        ylabel(r'$\rm Number$')
#         grid(True)
    
#         ax.tick_params(axis='x', labelsize=10)
#         ax.tick_params(axis='y',labelsize=10)
        tight_layout()
#         ylim(0,5)

        if save:
            savefig('{0}/hist(Impact_vir).pdf'.format(saveDirectory),format='pdf')
        else:
            show()


##########################################################################################
##########################################################################################
    # make a histogram of the distribution of impact parameters for associated galaxies, 
    # and split into red and blue shifted bins
    #
    
    plotImpactHist_dif = False
    save = False
    
    if plotImpactHist_dif:
        fig = figure(figsize=(10,2))
        ax = fig.add_subplot(111)
        subplots_adjust(hspace=0.200)
        alpha = 0.75

        binsize = 50
        bins = arange(0,500+binsize,binsize)
        
        reds = np.array([])
        blues = np.array([])
        
        # make sure all the values are okay
        for i,v,d in zip(impactList,virList,difList):
            if isNumber(i) and isNumber(v):
                if i !=-99 and v !=-99:
                    val = float(i)
                    
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
        xlabel(r'$\rho$ (kpc)')
        legend(scatterpoints=1,prop={'size':12})
        ax.tick_params(axis='y',labelsize=11)
        ax.tick_params(axis='x', labelsize=0)
        ylim(0,10)
        tight_layout()
        
        if save:
            savefig('{0}/hist(Impact_blue)_bin_{1}.pdf'.format(saveDirectory,binsize),format='pdf')
        else:
            show()

        fig = figure(figsize=(10,2))
        ax = fig.add_subplot(111)
        plot2 = hist(reds,bins=bins,histtype='bar',color="red",label='Redshifted',alpha=alpha)
        xlabel(r'$\rho$ (kpc)')
        ylabel('Number')
        legend(scatterpoints=1,prop={'size':12})
        
        ax.tick_params(axis='x', labelsize=11)
        ax.tick_params(axis='y',labelsize=11)
        ylim(0,10)
        tight_layout()
        
        if save:
            savefig('{0}/hist(Impact_red)_bin_{1}.pdf'.format(saveDirectory,binsize),format='pdf')
        else:
            show()


##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################


if __name__=="__main__":
    # do the work
    main()
    