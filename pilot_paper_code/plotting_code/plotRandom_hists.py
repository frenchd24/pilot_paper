#!/usr/bin/env python

'''
By David French (frenchd@astro.wisc.edu)

$Id:  plotRandom_hists.py, v 1.0 04/27/2016


Plot correlation functions, etc for a random selection of Lya locations across the sky
    - (4/27/16)

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
        saveDirectory = '/Users/David/Research_Documents/inclination/git_inclination/pilot_paper_code/plots3/'
        galaxyFilename = '/Users/David/Research_Documents/gt/NewGalaxyTable5.csv'
        pickleFilename = '/Users/David/Research_Documents/inclination/rand_pickle2.p'

    elif getpass.getuser() == 'frenchd':
        saveDirectory = '/usr/users/frenchd/inclination/git_inclination/pilot_paper_code/plots3/'
        galaxyFilename = '/usr/users/frenchd/gt/NewGalaxyTable5.csv'
        pickleFilename = '/usr/users/frenchd/inclination/rand_pickle2.p'


    else:
        print 'Could not determine username. Exiting.'
        sys.exit()
    
    # use the old pickle file to get the full galaxy dataset info
    pickleFile = open(pickleFilename,'rU')
    data = pickle.load(pickleFile)
    pickleFile.close()
    
#     galFile = open(galaxyFilename,'rU')
#     reader = csv.DictReader(galFile)

    # now randomly generate some locations in the sky

    randVel,randLike,randRho,randImpact_rho,randAz,randInc,randImpact_nolike,randImpact_rho_nolike,randDelV_nolike = data
    
#     print
#     print 'randImpact_nolike: ',randImpact_nolike
#     print
#     print 'randImpact_rho: ',randImpact_rho
#     print
#     print 'randImpact_rho_nolike: ',randImpact_rho_nolike

##########################################################################################
##########################################################################################
    # make a histogram of the distribution of impact/R_vir for random locations
    #
    
    plotRandHist_imp_rho = True
    save = False
    
    if plotRandHist_imp_rho:
        fig = figure(figsize=(10,4))
        ax = fig.add_subplot(111)
        alpha = 0.75

        bins = arange(0,15,0.1)

        print 'median: ',median(randImpact_rho_nolike)
        print 'mean: ',mean(randImpact_rho_nolike)
        
        plot1 = hist(randImpact_rho_nolike,bins=bins,histtype='bar',alpha=alpha)
        
        xlabel(r'$\rm \rho / R_{vir}$')
        ylabel(r'$\rm Number$')
#         grid(True)
    
#         ax.tick_params(axis='x', labelsize=10)
#         ax.tick_params(axis='y',labelsize=10)
        tight_layout()
#         ylim(0,5)

        if save:
            savefig('{0}/hist(random_impact_vir)full_1000.pdf'.format(saveDirectory,num),format='pdf')
        else:
            show()
            
##########################################################################################
##########################################################################################
    # make a CDF of the distribution of impact/R_vir for random locations
    #
    
    plotRandCDF_imp_rho = True
    save = False
    
    if plotRandCDF_imp_rho:
        fig = figure(figsize=(8,8))
        ax = fig.add_subplot(111)
        alpha = 0.75

        bins = arange(0,1000,10)

        print 'median: ',median(randImpact_nolike)
        print 'mean: ',mean(randImpact_nolike)
        
        sorted_data = sort(randImpact_rho_nolike)
    
        # compute the CDF y-values
        yvals=np.arange(len(sorted_data)) / float(len(sorted_data))
    
        # plot the y-values against the sorted data
        plt.plot(sorted_data,yvals)
                
        xlabel(r'$\rm \rho / R_{vir}$')
        ylabel(r'$\rm Number$')
#         grid(True)
    
#         ax.tick_params(axis='x', labelsize=10)
#         ax.tick_params(axis='y',labelsize=10)
        tight_layout()
#         ylim(0,5)

        if save:
            savefig('{0}/CDF(random_impact_vir)full_1000.pdf'.format(saveDirectory,num),format='pdf')
        else:
            show()

##########################################################################################
##########################################################################################
    # make a histogram of the distribution of impact parameter for random locations
    #
    
    plotRandHist_imp = True
    save = False
    
    if plotRandHist_imp:
        fig = figure(figsize=(10,4))
        ax = fig.add_subplot(111)
        alpha = 0.75

        bins = arange(0,1000,10)

        print 'median: ',median(randImpact_nolike)
        print 'mean: ',mean(randImpact_nolike)
        
        plot1 = hist(randImpact_nolike,bins=bins,histtype='bar',alpha=alpha)
        
        xlabel(r'$\rm \rho$')
        ylabel(r'$\rm Number$')
#         grid(True)
    
#         ax.tick_params(axis='x', labelsize=10)
#         ax.tick_params(axis='y',labelsize=10)
        tight_layout()
#         ylim(0,5)

        if save:
            savefig('{0}/hist(random_impact)full_1000.pdf'.format(saveDirectory,num),format='pdf')
        else:
            show()


##########################################################################################
##########################################################################################
    # make a CDF of the distribution of impact parameter for random locations
    #
    
    plotRandCDF_imp = True
    save = False
    
    if plotRandCDF_imp:
        fig = figure(figsize=(8,8))
        ax = fig.add_subplot(111)
        alpha = 0.75

        bins = arange(0,1000,10)

        print 'median: ',median(randImpact_nolike)
        print 'mean: ',mean(randImpact_nolike)
        
        sorted_data = sort(randImpact_nolike)
    
        # compute the CDF y-values
        yvals=np.arange(len(sorted_data)) / float(len(sorted_data))
    
        # plot the y-values against the sorted data
        plt.plot(sorted_data,yvals)
                
        xlabel(r'$\rm \rho$')
        ylabel(r'$\rm Number$')
#         grid(True)
    
#         ax.tick_params(axis='x', labelsize=10)
#         ax.tick_params(axis='y',labelsize=10)
        tight_layout()
#         ylim(0,5)

        if save:
            savefig('{0}/CDF(random_impact)full_1000.pdf'.format(saveDirectory,num),format='pdf')
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
    