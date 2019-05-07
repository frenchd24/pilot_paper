#!/usr/bin/env python

'''
By David French (frenchd@astro.wisc.edu)

$Id:  plotW_az_hist.py, v 5.2 04/21/2016

Sum the EW into bins in azimuth and average, plotting the step histogram over the data

 - doesn't yield anything too interesting


This is adapted from the plotW_Az bit from histograms3.py. Now loads in a pickle
file of the relevant data, as created by "buildDataLists.py"

Previous (from histograms3.py):
    Plot some stuff for the 100largest initial results

    Make plots for AAS winter 2014 poster. Uses LG_correlation_combined2.csv file

    Updated for the pilot paper (05/06/15)

v5: Updated for the final pilot paper results (12/04/15)
    - changed name from plotW_newAz_hist2.py to plotW_az_hist.py, and combined with
    plotW_newAz_hist2.py because they were similar
    (1/4/16)


v5.1: updated for LG_correlation_combined5_8_edit2.csv and l_min = 0.001 (02/24/2016)

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
    # plot equivalent width as a function of azimuth angle for red vs blue
    # shifted absorption
    #
    
    plotW_Az_avg = False
    save = False
    
    if plotW_Az_avg:
        fig = figure()
        ax = fig.add_subplot(111)
        
        countb = 0
        countr = 0
        count = -1
        
        labelr = 'Red Shifted Absorber'
        labelb = "Blue Shifted Absorber"

        print 'len(azList): ',len(azList)
        print
        print azList
        
        placeArrayr = zeros(10)
        placeCountr = zeros(10)
        placeArrayb = zeros(10)
        placeCountb = zeros(10)
        
        for d,a,w,m in zip(difList,azList,lyaWList,majList):
            # check if all the values are good
            if isNumber(d) and isNumber(a) and isNumber(w) and isNumber(m):
                if d!=-99 and a!=-99 and w!=-99 and m!=-99:
                    if d>0:
                        # galaxy is behind absorber, so gas is blue shifted
                        color = 'Blue'
                        
                        # which bin does it belong too?
                        place = a/10
                        print 'place: ',place
                        placeArrayb[place] += float(w)
                        print 'placeArrayb: ',placeArrayb
                        placeCountb[place] +=1.
                        print 'placecountb: ',placeCountb
                        
                        if countb == 0:
                            countb +=1
                            plotb = ax.scatter(a,w,c='Blue',s=50,label= labelb)
                    if d<0:
                        # gas is red shifted compared to galaxy
                        color = 'Red'
                        
                        # which bin does it belong too?
                        place = a/10
                        placeArrayr[place] += float(w)
                        placeCountr[place] +=1.
                        
                        if countr == 0:
                            countr +=1
                            plotr = ax.scatter(a,w,c='Red',s=50,label= labelr)
            
                    plot1 = scatter(a,w,c=color,s=50)
                    
        rHist = placeArrayr/placeCountr
        print 'rHist: ',rHist
        bHist = placeArrayb/placeCountb
        print 'bHist: ',bHist
        
#         rHist = placeArrayr
#         bHist = placeArrayb
        
        totalrHist = []
        totalrAz = []
        totalbHist = []
        totalbAz = []
        
        for r,a in zip(rHist,range(0,100,10)):
            totalrHist.append(r)
            totalrHist.append(r)

            totalrAz.append(a)
            totalrAz.append(a+10)

#             plot2 = ax.plot([a,a+10],[r,r],c='Red')
            
        for b,a in zip(bHist,range(0,100,10)):
            totalbHist.append(b)
            totalbHist.append(b)

            totalbAz.append(a)
            totalbAz.append(a+10)

#             plot2 = ax.plot([a,a+10],[b,b],c='Blue')
            
        print 'totalbHist: ',totalbHist
        print
        print 'totalbAz: ',totalbAz
        
        
        print 'totalbHist: ',totalbHist
        print
        print 'totalbAz: ',totalbAz
        
        plot2 = ax.plot(totalrAz,totalrHist,c='Red',lw=5)
        plot3 = ax.plot(totalbAz,totalbHist,c='Blue',lw=5)
        
        title('W(azimuth) for red vs blue shifted absorption')
        xlabel(r'Azimuth (deg)')
        ylabel(r'Equivalent Width ($\rm m\AA$)')
        legend(scatterpoints=1)
        ax.grid(b=None,which='major',axis='both')
#         ylim(0,1300)
#         xlim(0,90)

        if save:
            savefig('{0}/W(azimuth)_dif_avgHistograms.pdf'.format(saveDirectory),format='pdf')
        else:
            show()
            
            
            
##########################################################################################
##########################################################################################
    # plot equivalent width as a function of azimuth angle for red vs blue
    # shifted absorption, include median histograms
    #
    # DOESN'T WORK YET
    
    plotW_Az_med = True
    save = False
    
    if plotW_Az_med:
        fig = figure()
        ax = fig.add_subplot(111)
        
        countb = 0
        countr = 0
        count = -1
        
        labelr = 'Redshifted Absorber'
        labelb = "Blueshifted Absorber"
        alpha = 0.75

        print 'len(azList): ',len(azList)
        print
        print azList
        
        placeArrayr = zeros(10)
        placeCountr = zeros(10)
        placeArrayb = zeros(10)
        placeCountb = zeros(10)
        
        for d,a,w,m in zip(difList,azList,lyaWList,majList):
            # check if all the values are good
            if isNumber(d) and isNumber(a) and isNumber(w) and isNumber(m):
                if d!=-99 and a!=-99 and w!=-99 and m!=-99:
                    if d>0:
                        # galaxy is behind absorber, so gas is blue shifted
                        color = 'Blue'
                        
                        # which bin does it belong too?
                        place = a/10
                        placeArrayb[place] += float(w)
                        placeCountb[place] +=1.
                        
                        if countb == 0:
                            countb +=1
                            plotb = ax.scatter(a,w,c='Blue',s=50,label= labelb,alpha=alpha)
                    if d<0:
                        # gas is red shifted compared to galaxy
                        color = 'Red'
                        
                        # which bin does it belong too?
                        place = a/10
                        placeArrayr[place] += float(w)
                        placeCountr[place] +=1.
                        
                        if countr == 0:
                            countr +=1
                            plotr = ax.scatter(a,w,c='Red',s=50,label= labelr,alpha=alpha)
            
                    plot1 = scatter(a,w,c=color,s=50,alpha=alpha)
                    
        rHist = placeArrayr/placeCountr
        print 'rHist: ',rHist
        bHist = placeArrayb/placeCountb
        print 'bHist: ',bHist
        
#         rHist = placeArrayr
#         bHist = placeArrayb
        
        totalrHist = []
        totalrAz = []
        totalbHist = []
        totalbAz = []
        
        for r,a in zip(rHist,range(0,100,10)):
            totalrHist.append(r)
            totalrHist.append(r)

            totalrAz.append(a)
            totalrAz.append(a+10)

#             plot2 = ax.plot([a,a+10],[r,r],c='Red')
            
        for b,a in zip(bHist,range(0,100,10)):
            totalbHist.append(b)
            totalbHist.append(b)

            totalbAz.append(a)
            totalbAz.append(a+10)

#             plot2 = ax.plot([a,a+10],[b,b],c='Blue')
            
        print 'totalbHist: ',totalbHist
        print
        print 'totalbAz: ',totalbAz
        
        
        print 'totalbHist: ',totalbHist
        print
        print 'totalbAz: ',totalbAz
        
#         plot2 = ax.plot(totalrAz,totalrHist,c='Red',lw=5)
#         plot3 = ax.plot(totalbAz,totalbHist,c='Blue',lw=5)
        
        title('W(azimuth) for red vs blue shifted absorption')
        xlabel(r'Azimuth (deg)')
        ylabel(r'Equivalent Width ($\rm m\AA$)')
        legend(scatterpoints=1)
        ax.grid(b=None,which='major',axis='both')
#         ylim(0,1300)
#         xlim(0,90)

        if save:
            savefig('{0}/W(azimuth)_dif_avgHistograms.pdf'.format(saveDirectory),format='pdf')
        else:
            show()

##########################################################################################
##########################################################################################
    # plot equivalent width as a function of azimuth angle for red vs blue
    # shifted absorption
    #
    
    plotW_Az_total = False
    save = False
    
    if plotW_Az_total:
        fig = figure()
        ax = fig.add_subplot(111)
        
        countb = 0
        countr = 0
        count = -1
        
        labelr = 'Red Shifted Absorber'
        labelb = "Blue Shifted Absorber"

        print 'len(azList): ',len(azList)
        print
        print azList
        
        placeArrayr = zeros(10)
        placeCountr = zeros(10)
        placeArrayb = zeros(10)
        placeCountb = zeros(10)
        
        for d,a,w,m in zip(difList,azList,lyaWList,majList):
            # check if all the values are good
            if isNumber(d) and isNumber(a) and isNumber(w) and isNumber(m):
                if d!=-99 and a!=-99 and w!=-99 and m!=-99:
                    if d>0:
                        # galaxy is behind absorber, so gas is blue shifted
                        color = 'Blue'
                        
                        # which bin does it belong too?
                        place = a/10
                        print 'place: ',place
                        placeArrayb[place] += float(w)
                        print 'placeArrayb: ',placeArrayb
                        placeCountb[place] +=1.
                        print 'placecountb: ',placeCountb
                        
                        if countb == 0:
                            countb +=1
                            plotb = ax.scatter(a,w,c='Blue',s=50,label= labelb)
                    if d<0:
                        # gas is red shifted compared to galaxy
                        color = 'Red'
                        
                        # which bin does it belong too?
                        place = a/10
                        placeArrayr[place] += float(w)
                        placeCountr[place] +=1.
                        
                        if countr == 0:
                            countr +=1
                            plotr = ax.scatter(a,w,c='Red',s=50,label= labelr)
            
                    plot1 = scatter(a,w,c=color,s=50)
                    
        rHist = placeArrayr
        print 'rHist: ',rHist
        bHist = placeArrayb
        print 'bHist: ',bHist
        
        totalrHist = []
        totalrAz = []
        totalbHist = []
        totalbAz = []
        
        for r,a in zip(rHist,range(0,100,10)):
            totalrHist.append(r)
            totalrHist.append(r)

            totalrAz.append(a)
            totalrAz.append(a+10)

#             plot2 = ax.plot([a,a+10],[r,r],c='Red')
            
        for b,a in zip(bHist,range(0,100,10)):
            totalbHist.append(b)
            totalbHist.append(b)

            totalbAz.append(a)
            totalbAz.append(a+10)

#             plot2 = ax.plot([a,a+10],[b,b],c='Blue')
            
        print 'totalbHist: ',totalbHist
        print
        print 'totalbAz: ',totalbAz
        
        
        print 'totalbHist: ',totalbHist
        print
        print 'totalbAz: ',totalbAz
        
        plot2 = ax.plot(totalrAz,totalrHist,c='Red',lw=5)
        plot3 = ax.plot(totalbAz,totalbHist,c='Blue',lw=5)
        
        title('W(azimuth) for red vs blue shifted absorption')
        xlabel(r'Azimuth (deg)')
        ylabel(r'Equivalent Width ($\rm m\AA$)')
        legend(scatterpoints=1)
        ax.grid(b=None,which='major',axis='both')
#         ylim(0,2100)
#         xlim(0,90)

        if save:
            savefig('{0}/W(azimuth)_dif_totalHistograms.pdf'.format(saveDirectory),format='pdf')
        else:
            show()

##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################


if __name__=="__main__":
    # do the work
    main()
    