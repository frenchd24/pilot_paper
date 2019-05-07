#!/usr/bin/env python

'''
By David French (frenchd@astro.wisc.edu)

$Id:  pilot_stats.py, v 1.4 09/23/16

Print out all the relevant stats on the dataset for the pilot paper (01/04/2016)

v1.1: updates for v_hel velocity and probably something else? (2/22/16)

v1.2: updates for the new large galaxy sample (07/14/16) -> /plots4/

v1.3: added ability to limit results by environment variable (7/14/16)
    also same for likelihood values, including some stats about them also
    
v1.4: updated to LG_correlation_combined5_11_25cut_edit4.csv and plots5 (9/23/16)

'''

import sys
import os
import csv

from pylab import *
from scipy import stats
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
#         resultsFilename = '/Users/David/Research_Documents/inclination/git_inclination/LG_correlation_combined5_3.csv'
#         resultsFilename = '/Users/David/Research_Documents/inclination/git_inclination/LG_correlation_combined5_8_edit2.csv'
        resultsFilename = '/Users/David/Research_Documents/inclination/git_inclination/LG_correlation_combined5_11_25cut_edit4.csv'
        saveDirectory = '/Users/David/Research_Documents/inclination/git_inclination/pilot_paper_code/plots5/'

    elif getpass.getuser() == 'frenchd':
        pickleFilename = '/usr/users/frenchd/inclination/git_inclination/pilot_paper_code/pilotData2.p'
#         resultsFilename = '/usr/users/frenchd/inclination/git_inclination/LG_correlation_combined5_3.csv'
#         resultsFilename = '/usr/users/frenchd/inclination/git_inclination/LG_correlation_combined5_8_edit2.csv'
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
    
    maxEnv = 3000
    minL = 0.001
    
    # if match, then the includes in the file have to MATCH the includes above. e.g., if 
    # virInclude = False, cusInclude = True, finalInclude = False, then only systems
    # matching those three would be included. Otherwise, all cusInclude = True would be included
    # regardless of the others
    match = False
    
    # all the lists to be used for associated lines
    nameList = []
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
    AGNnameList = []
    
    # for ambiguous lines
    lyaVAmbList = []
    lyaWAmbList = []
    envAmbList = []
    ambAGNnameList = []
    
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
            AGNname = l['AGNname']
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
            morph = l['final_morphology']
            vcorr = l['vcorrGalaxy (km/s)']
            maj = l['majorAxis (kpc)']
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
            if float(env) <= maxEnv and float(likelihood) >=minL:
                nameList.append(galaxyName)
                AGNnameList.append(AGNname)
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
            AGNname = l['AGNname']
             
        
            lyaVAmbList.append(float(lyaV))
            lyaWAmbList.append(float(lyaW))
            envAmbList.append(float(env))
            ambAGNnameList.append(AGNname)

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
#########################################################################################
    
    # print all the things
    #
    
    # absorber info lists
    blues = []
    reds = []
    blueAbs = []
    redAbs = []
    blueW = []
    redW = []
    blueB = []
    redB = []
    blueErr = []
    redErr = []
    blueV = []
    redV = []
    blueImpact = []
    redImpact = []
    
    # galaxy info lists
    blueInc = []
    redInc = []
    blueFancyInc = []
    redFancyInc = []
    blueAz = []
    redAz = []
    bluePA = []
    redPA = []
    blueVcorr = []
    redVcorr = []
    blueEnv = []
    redEnv = []
    blueVir = []
    redVir = []
    blueLike = []
    redLike = []
    

    # ambiguous stuff
    void = []
    ambig = []
    for v,w,e in zip(lyaVAmbList,lyaWAmbList,envAmbList):
        if e == 0:
            void.append(w)
        else:
            ambig.append(w)
    
    
    # for targets
    finalTargets = {}
    for a in AGNnameList:
        if finalTargets.has_key(a):
            i = finalTargets[a]
            i+=1
            finalTargets[a] = i
            
        else:
            finalTargets[a] = 1
            
    # for ambiguous targets
    ambTargets = {}
    for a in ambAGNnameList:
        if ambTargets.has_key(a):
            i = ambTargets[a]
            i+=1
            ambTargets[a] = i
            
        else:
            ambTargets[a] = 1
        
    
    # for absorbers
    for d,w,e,v,i,b in zip(difList,lyaWList,lyaErrList,lyaVList,impactList,bList):
        if d>=0:
            blues.append(float(d))
            blueW.append(float(w))
            blueErr.append(float(e))
            blueV.append(float(v))
            blueImpact.append(float(i))
            blueAbs.append(abs(d))
            blueB.append(float(b))
        else:
            reds.append(float(d))
            redW.append(float(w))
            redErr.append(float(e))
            redV.append(float(v))
            redImpact.append(float(i))
            redAbs.append(abs(d))
            redB.append(float(b))
            
            
##########################################################################################
    blueSpiralInc = []
    redSpiralInc = []
    spiralIncList = []
    # for spirals only
    for d,inc in zip(difList,fancyIncList):
        spiralIncList.append(float(inc))
        if d>=0:
            blueSpiralInc.append(float(inc))
        else:
            redSpiralInc.append(float(inc))
                
                
    # compile a list of only spiral galaxy inclinations from the full galaxy table
    if getpass.getuser() == 'David':
        galaxyFile = open('/Users/David/Research_Documents/gt/NewGalaxyTable5.csv','rU')
    else:
        print 'Not on laptop, exiting'
        sys.exit()
        
    reader = csv.DictReader(galaxyFile)
    
    allDiameters = []
    incGT25diam = []
    allSpiralIncList = []
    
    q0 = 0.2
    for i in reader:
        major,minor = eval(i['linDiameters (kpc)'])
        morph = i['morphology'].lower()
        if bfind(morph,'s'):
            if not bfind(morph,'sph') and not bfind(morph,'s0'):
            
                if isNumber(major):
                    if isNumber(minor):
                        if float(major) > float(minor):
                            fInc = calculateFancyInclination(major,minor,q0)
                            allSpiralIncList.append(fInc)
                            
                            if float(major) >=25.0:
                                incGT25diam.append(fInc)
    
    galaxyFile.close()
    
##########################################################################################
                   
            
    nameDict = {}
    # for galaxies
    for d,inc,finc,az,pa,vcorr,e,vir,l,name in zip(difList,incList,fancyIncList,azList,paList,vcorrList,envList,virList, likeList,nameList):
        if nameDict.has_key(name):
            i = nameDict[name]
            i+=1
            nameDict[name] = i
        else:
            nameDict[name] = 1
        
        if d>=0:
            if inc !=-99:
                blueInc.append(float(inc))
            if finc !=-99:
                blueFancyInc.append(float(finc))
            if az !=-99:
                blueAz.append(float(az))
            if pa !=-99:
                bluePA.append(float(pa))
            if vcorr !=-99:
                blueVcorr.append(float(vcorr))
            blueEnv.append(float(e))
            if vir !=-99:
                blueVir.append(float(vir))
            if l !=-99:
                blueLike.append(float(l))
        else:
            if inc !=-99:
                redInc.append(float(inc))
            if finc !=-99:
                redFancyInc.append(float(finc))
            if az !=-99:
                redAz.append(float(az))
            if pa !=-99:
                redPA.append(float(pa))
            if vcorr !=-99:
                redVcorr.append(float(vcorr))
            redEnv.append(float(e))
            if vir !=-99:
                redVir.append(float(vir))
            if l !=-99:
                redLike.append(float(l))
                
    galaxyNames = nameDict.keys()
                
    # how many absorbers above vs below vel_cut?
    redVelCount200 = 0
    redVelCount100 = 0
    blueVelCount200 = 0
    blueVelCount100 = 0
    
    for b in blues:
        if b >=200:
            blueVelCount200 +=1
        if b >= 100:
            blueVelCount100 +=1
        
    for r in reds:
        if abs(r) >=200:
            redVelCount200 +=1
        if abs(r) >=100:
            redVelCount100 +=1
    

    assocFancyInc = blueFancyInc + redFancyInc

    
    print
    print '------------------------ Pilot Data ------------------------------'
    print
    print ' FOR THE FOLLOWING INCLUDE SET:'
    print ' Virial radius include = ',virInclude
    print ' Custom include =        ',cusInclude
    print ' Final include =         ',finalInclude
    print ' Match =                 ',match
    print
    print 'total number of lines: ', len(lyaWList) + len(lyaWAmbList)
    print 'total number of unique galaxies matched: ',len(galaxyNames)
    print 'total number of associated lines: ',len(difList)
    print 'total number of ambiguous lines: ',len(ambig)
    print 'total number of void lines: ',len(void)
    print '# of redshifted lines: ',len(reds)
    print '# of blueshifted lines: ',len(blues)
    print
    print
    print ' ASSOCIATED TARGETS '
    print
    print 'final target number: ',len(finalTargets.keys())
    for i in finalTargets.keys():
        print i
    print
    print
    print ' AMBIGUOUS TARGTS '
    print
    print 'final ambiguous number: ',len(ambTargets.keys())
    for i in ambTargets.keys():
        print i
    print
    print
    print '----------------------- Absorber info ----------------------------'
    print
    print 'avg blueshifted EW: ',mean(blueW)
    print 'median blueshifted EW: ',median(blueW)
    print 'avg blue err: ',mean(blueErr)
    print 'median blue err: ',median(blueErr)
    print
    print 'std(blue EW): ',std(blueW)
    print 'stats.sem(blue EW): ',stats.sem(blueW)
    print 'stats.describe(blue EW): ',stats.describe(blueW)
    print
    print 'avg blueshifted vel_diff: ',mean(blues)
    print 'median blueshifted vel_diff: ',median(blues)
    print 'std(blueshifted vel_diff): ',std(blues)
    print 'stats.sem(blue vel_dif): ',stats.sem(blues)
    print 'stats.describe(blue vel_dif: ',stats.describe(blues)
    print
    print '% blueshifted which have vel_diff >= 200 km/s: {0}'.format(float(blueVelCount200)/len(blues))
    print 'total number with abs(vel_diff) >= 200 km/s: {0}'.format(blueVelCount200)
    print '% blueshifted which have vel_diff >= 100 km/s: {0}'.format(float(blueVelCount100)/len(blues))
    print 'total number with abs(vel_diff) >= 100 km/s: {0}'.format(blueVelCount100)
    print
    
    print 'avg blue velocity: ',mean(blueV)
    print 'median blue velocity: ',median(blueV)
    print 'std(blue Velocity): ',std(blueV)
    print 'avg blue impact: ',mean(blueImpact)
    print 'median blue impact: ',median(blueImpact)
    print 'stats.sem(blue impact): ',stats.sem(blueImpact)
    print 'stats.describe(blue impact): ',stats.describe(blueImpact)

    print
    
    print 'avg redshifted EW: ',mean(redW)
    print 'median redshifted EW: ',median(redW)
    print 'avg red err: ',mean(redErr)
    print 'median red err: ',median(redErr)
    print
    print 'std(red EW): ',std(redW)
    print 'stats.sem(red EW): ',stats.sem(redW)
    print 'stats.describe(red EW): ',stats.describe(redW)

    print
    print 'avg redshifted vel_diff: ',mean(reds)
    print 'median redshifted vel_diff: ',median(reds)
    print 'std(redshifted vel_dif): ',std(reds)
    print 'stats.sem(red vel_dif): ',stats.sem(reds)
    print 'stats.describe(red vel_dif): ',stats.describe(reds)
    print
    print '% redshifted which have abs(vel_diff) >= 200 km/s: {0}'.format(float(redVelCount200)/len(reds))
    print 'total number with abs(vel_diff) >= 200 km/s: {0}'.format(redVelCount200)
    print '% redshifted which have abs(vel_diff) >= 100 km/s: {0}'.format(float(redVelCount100)/len(reds))
    print 'total number with abs(vel_diff) >= 100 km/s: {0}'.format(redVelCount100)
    print

    print 'avg red velocity: ',mean(redV)
    print 'median red velocity: ',median(redV)
    print
    print 'avg red impact: ',mean(redImpact)
    print 'median red impact: ',median(redImpact)
    print 'stats.sem(red impact): ',stats.sem(redImpact)
    print 'stats.describe(red impact): ',stats.describe(redImpact)
    print 'std(red impact): ',std(redImpact)



    print
    print '----------------------- Galaxy info ----------------------------'
    print
    
    # regular inclinations
    incCut = 50
    totalBlueInc = len(blueInc)
    totalRedInc = len(redInc)
    
    blueIncCount = 0
    for i in blueInc:
        if i >= incCut:
            blueIncCount +=1
            
    redIncCount = 0
    for i in redInc:
        if i >= incCut:
            redIncCount +=1
            
    totalInc = len(allInclinations)
    totalCount = 0
    for i in allInclinations:
        if i >= incCut:
            totalCount +=1
            
            
    # fancy inclinations
    totalBlueFancyInc = len(blueFancyInc)
    totalRedFancyInc = len(redFancyInc)
    
    blueFancyIncCount = 0
    for i in blueFancyInc:
        if i >= incCut:
            blueFancyIncCount +=1
            
    redFancyIncCount = 0
    for i in redFancyInc:
        if i >= incCut:
            redFancyIncCount +=1
            
    combinedCount = redFancyIncCount + blueFancyIncCount
    totalCombinedCount = totalRedFancyInc + totalBlueFancyInc
            
    totalFancyInc = len(allFancyInclinations)
    totalFancyCount = 0
    for i in allFancyInclinations:
        if i >= incCut:
            totalFancyCount +=1
    
    print
    print ' INCLINATIONS: '
    print 
    print 'Blue: {0} % of associated galaxies have >={1}% inclination'.format(float(blueIncCount)/float(totalBlueInc),incCut)
    print 'Red: {0} % of associated galaxies have >={1}% inclination'.format(float(redIncCount)/float(totalRedInc),incCut)
    print 'All: {0} % of ALL galaxies have >={1}% inclination'.format(float(totalCount)/float(totalInc),incCut)
    print
    print ' FANCY INCLINATIONS: '
    print
    print 'Blue: {0} % of associated galaxies have >={1}% fancy inclination'.format(float(blueFancyIncCount)/float(totalBlueFancyInc),incCut)
    print 'Red: {0} % of associated galaxies have >={1}% fancy inclination'.format(float(redFancyIncCount)/float(totalRedFancyInc),incCut)
    print 'All: {0} % of ALL galaxies have >={1}% fancy inclination'.format(float(totalFancyCount)/float(totalFancyInc),incCut)
    print 'Combined: {0} % of associated galaxies have >= {1} fancy inclination'.format(float(combinedCount)/float(totalCombinedCount),incCut)
    print
    print 'Average all fancy inclination: ',mean(allFancyInclinations)
    print 'stats.sem(all): ',stats.sem(allFancyInclinations)
    print    
    print 'avg blue inclination: ',mean(blueInc)
    print 'median blue inclination: ',median(blueInc)
    print 'avg blue fancy inclination: ',mean(blueFancyInc)
    print 'median blue fancy inclination: ',median(blueFancyInc)
    print
    print 'avg red inclination: ',mean(redInc)
    print 'median red inclination: ',median(redInc)
    print 'avg red fancy inclination: ',mean(redFancyInc)
    print 'median red fancy inclination: ',median(redFancyInc)
    
    print
    print 'mean associated: ',mean(assocFancyInc)
    print 'stats.sem(associated): ',stats.sem(assocFancyInc)
    print 'stats.describe(associated): ',stats.describe(assocFancyInc)
    print 'stats.sem(blue): ',stats.sem(blueFancyInc)
    print 'stats.describe(blue): ',stats.describe(blueFancyInc)
    print
    print 'stats.sem(red): ',stats.sem(redFancyInc)
    print 'stats.describe(red): ',stats.describe(redFancyInc)
    
    print
    print "  AZIMUTHS and PA:  "
    print
    print 'avg blue azimuth: ',mean(blueAz)
    print 'median blue azimuth: ',median(blueAz)
    print 'stats.sem(blue az): ',stats.sem(blueAz)
    print 'stats.describe(blue az): ',stats.describe(blueAz)
    print
    print 'avg red azimuth: ',mean(redAz)
    print 'median red azimuth: ',median(redAz)
    print 'stats.sem(red az): ',stats.sem(redAz)
    print 'stats.describe(red az): ',stats.describe(redAz)
    print
    print 'avg blue PA: ',mean(bluePA)
    print 'median blue PA: ',median(bluePA)
    print
    print 'avg red PA: ',mean(redPA)
    print 'median red PA: ',median(redPA)
    
    print
    print ' VCORR : '
    print
    print 'avg blue vcorr: ',mean(blueVcorr)
    print 'median blue vcorr: ',median(blueVcorr)
    print
    print 'avg red vcorr: ',mean(redVcorr)
    print 'median red vcorr: ',median(redVcorr)
    
    print
    print ' ENVIRONMENT: '
    print
    print 'avg blue environment: ',mean(blueEnv)
    print 'median blue environment: ',median(blueEnv)
    print
    print 'avg red environment: ',mean(redEnv)
    print 'median red environment: ',median(redEnv)
    
    print
    print ' R_vir: '
    print
    print 'avg blue R_vir: ',mean(blueVir)
    print 'median blue R_vir: ',median(blueVir)
    print 'stats.sem(blue R_vir): ',stats.sem(blueVir)
    print 'stats.describe(blue R_vir): ',stats.describe(blueVir)
    print
    print 'avg red R_vir: ',mean(redVir)
    print 'median red R_vir: ',median(redVir)
    print 'stats.sem(red R_vir): ',stats.sem(redVir)
    print 'stats.describe(red R_vir): ',stats.describe(redVir)

    print
    print ' LIKELIHOOD: '
    print
    print 'avg blue likelihood: ',mean(blueLike)
    print 'median blue likelihood: ',median(blueLike)
    print
    print 'avg red likelihood: ',mean(redLike)
    print 'median red likelihood: ',median(redLike)
    
    print
    print
    print '-------------------- Distribution analysis ----------------------'
    print
    print
    
    print ' FANCY INCLINATIONS: '
    
    # perform the K-S and AD tests for inclination
    ans1 = stats.ks_2samp(blueFancyInc, redFancyInc)
    ans1a = stats.anderson_ksamp([blueFancyInc,redFancyInc])

    print 'KS for blue vs red fancy inclinations: ',ans1
    print 'AD for blue vs red fancy inclinations: ',ans1a
    
    ans2 = stats.ks_2samp(blueFancyInc, allFancyInclinations)
    print 'KS for blue vs all fancy inclinations: ',ans2
    
    ans3 = stats.ks_2samp(redFancyInc, allFancyInclinations)
    print 'KS for red vs all fancy inclinations: ',ans3
    
    print
    z_statrb, p_valrb = stats.ranksums(blueFancyInc, redFancyInc)
    z_statall, p_valall = stats.ranksums(assocFancyInc, allFancyInclinations)
    print 'ranksum red vs blue p-value: ',p_valrb
    print 'ranksum associated vs all: ',p_valall


#     ans4 = stats.ks_2samp(assocFancyInc, allFancyInclinations)
#     ans4a = stats.anderson_ksamp([assocFancyInc,allFancyInclinations])
# 
#     print 'KS for all associated vs all fancy inclinations: ',ans4
#     print 'AD for all associated vs all fancy inclinations: ',ans4a
#     
    print

#     ans5 = stats.ks_2samp(spiralIncList, allSpiralIncList)
#     ans5a = stats.anderson_ksamp([spiralIncList,allSpiralIncList])
# 
#     print 'KS for all spiral associated vs all spiral fancy inclinations: ',ans5
#     print 'AD for all spiral associated vs all spiral fancy inclinations: ',ans5a
    
    print
    print ' INCLINATIONS: '
    print
    
    # perform the K-S and AD tests for inclination
    ans1 = stats.ks_2samp(blueInc, redInc)
    ans1a = stats.anderson_ksamp([blueInc,redInc])

    print 'KS for blue vs red inclinations: ',ans1
    print 'AD for blue vs red inclinations: ',ans1a
    
    ans2 = stats.ks_2samp(blueInc, allInclinations)
    print 'KS for blue vs all inclinations: ',ans2
    
    ans3 = stats.ks_2samp(redInc, allInclinations)
    print 'KS for red vs all inclinations: ',ans3
    
    assocInc = blueInc + redInc
    ans4 = stats.ks_2samp(assocInc, allInclinations)
    print 'KS for associated vs all inclinations: ',ans4
    
    print
    print ' EW Distributions: '
    print
    
    # perform the K-S and AD tests for EW
    ans1 = stats.ks_2samp(blueW, redW)
    ans1a = stats.anderson_ksamp([blueW,redW])
    print 'KS for blue vs red EW: ',ans1
    print 'AD for blue vs red EW: ',ans1a
    

    print
    print ' Impact parameter Distributions: '
    print
    
    # perform the K-S and AD tests for impact parameter
    ans1 = stats.ks_2samp(blueImpact, redImpact)
    ans1a = stats.anderson_ksamp([blueImpact,redImpact])
    print 'KS for blue vs red impact parameters: ',ans1
    print 'AD for blue vs red impact parameters: ',ans1a
    
    print
    print ' \Delta v Distributions: '
    print
    
    # perform the K-S and AD tests for \delta v
    ans1 = stats.ks_2samp(blueAbs, redAbs)
    ans1a = stats.anderson_ksamp([blueAbs,redAbs])
    print 'KS for blue vs red \Delta v: ',ans1
    print 'AD for blue vs red \Delta v: ',ans1a
    
    print
    print ' Azimuth Distributions: '
    print
    
    # perform the K-S and AD tests for azimuth
    ans1 = stats.ks_2samp(blueAz, redAz)
    ans1a = stats.anderson_ksamp([blueAz,redAz])
    print 'KS for blue vs red azimuth: ',ans1
    print 'AD for blue vs red azimuth: ',ans1a
    print
    
    # now against a flat distribution
    flatRed = arange(0,90,1)
    flatBlue = arange(0,90,1)

    ans1 = stats.ks_2samp(blueAz, flatBlue)
    ans1a = stats.anderson_ksamp([blueAz,flatBlue])
    print 'KS for blue vs flat azimuth: ',ans1
    print 'AD for blue vs flat azimuth: ',ans1a
    print
    ans1 = stats.ks_2samp(redAz, flatRed)
    ans1a = stats.anderson_ksamp([redAz,flatRed])
    print 'KS for red vs flat azimuth: ',ans1
    print 'AD for erd vs flat azimuth: ',ans1a
    print
    
            
    print
    print ' Environment Distributions: '
    print
    
    # perform the K-S and AD tests for environment
    ans1 = stats.ks_2samp(blueEnv, redEnv)
    ans1a = stats.anderson_ksamp([blueEnv,redEnv])
    print 'KS for blue vs red environment: ',ans1
    print 'AD for blue vs red environment: ',ans1a
    
    print
    print ' R_vir Distributions: '
    print
    
    # perform the K-S and AD tests for r_vir
    ans1 = stats.ks_2samp(blueVir, redVir)
    ans1a = stats.anderson_ksamp([blueVir,redVir])
    print 'KS for blue vs red R_vir: ',ans1
    print 'AD for blue vs red R_vir: ',ans1a
    
    print
    print ' Doppler parameter Distributions: '
    print
    
    # perform the K-S and AD tests for doppler parameter
    ans1 = stats.ks_2samp(blueB, redB)
    ans1a = stats.anderson_ksamp([blueB,redB])
    print 'KS for blue vs red doppler parameter: ',ans1
    print 'AD for blue vs red doppler parameter: ',ans1a
    
    print
    print ' Likelihood Distributions: '
    print
    
    # perform the K-S and AD tests for doppler parameter
    ans1 = stats.ks_2samp(blueLike, redLike)
    ans1a = stats.anderson_ksamp([blueLike,redLike])
    print 'KS for blue vs red likelihood: ',ans1
    print 'AD for blue vs red likelihood: ',ans1a
    
    print
    print ' COMPLETED. '

    

###############################################################################
###############################################################################
###############################################################################
###############################################################################


if __name__=="__main__":
    # do the work
    main()
    