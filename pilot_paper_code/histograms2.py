#!/usr/bin/env python

'''
By David French (frenchd@astro.wisc.edu)

$Id:  histograms2.py, v 2.0 01/04/2015

Plot some stuff for the 100largest initial results

Make plots for AAS winter 2014 poster. Uses LG_correlation_combined2.csv file

'''

import sys
# sys.path.append("/usr/users/frenchd/local/lib/python2.7/site-packages/")
import os
import csv
# import string
# import warnings
# import urllib
# import numpy
from pylab import *
# import atpy
from math import *
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


def isNumber(s):
    try:
        float(s)
        return True
    except Exception:
        return False
        
def bfind(k,s):
    # better find (returns true/false)
    if k.find(s) == -1:
        return False
    else:
        return True

    
def parseGalaxyNames(nameList):
    # format galaxy names for the url
    
    newNameList = []
    for name in nameList:
        nname = name.strip()
        nname = urllib.quote_plus(nname)
        nname = nname.replace('\n','')
        newNameList.append(nname)
    return newNameList
    

def createCSVTable(outFile,fieldnames):
    # creates and returns a DictReader object populated with header names
    
    writer = csv.DictWriter(outFile, fieldnames=fieldnames)
    headers = dict((n,n) for n in fieldnames)
    writer.writerow(headers)
    return writer


def decidePreferredName(listOfAlternatives,oldName):
    # decides which name to use as the first name, returns that name
    
    print 
    bestName = None
    found = False
    
    # priority:
    first = ['NGC ','IC ','UGC ']
    second = ['MRK ','MCG ','ISO ','SBS ']
    third = ['CGCG ','IRAS ','RXJ ','FGC ','KUG ','PGC ','SDSS ','VCC ']
    fourth = ['2MASS ','2DF ','6DF ','HIPASS ','2MASX ']
    
    for f in first:
        if oldName.find(f) != -1:
            bestName = oldName
            found = True
            break
    
    if not found:
        for name in listOfAlternatives:
            for f in first:
                if name.find(f) != -1 and name.find('[') == -1:
                    bestName = name
                    found = True
                    break
                    
            if not found:
                for s in second:
                    if name.find(s) != -1 and name.find('[') == -1:
                        bestName = name
                        found = True            
                        break
                        
            if not found:
                for t in third:
                    if name.find(t) != -1 and name.find('[') == -1:
                        bestName = name
                        found = True
                        break
    
            if not found:
                for q in fourth:
                    if name.find(q) != -1 and name.find('[') == -1:
                        bestName = name
                        found = True
                        break
            if not found:
                bestName = oldName
    
    found2 = False
    if not found:
        for name in listOfAlternatives:
            if str(name) == oldName:
                bestName = name
                found2 = True
                found = True
                break
    
    if not found2:
        bestName = listOfAlternatives[0]
        
    print 'bestName: ',bestName
#   if bestName != oldName:
#       listOfAlternatives.remove(bestName)
#       listOfAlternatives.append(oldName)
    return bestName


def returnRAandDec(ra,dec):
    # this function converts ra and dec in degrees to HH:MM:SS and DD:MM:SS
    
    raHours,raMinutes = str(((float(ra)*24)/360)).split('.')
    raMinutes,raSeconds = str(float('0.'+raMinutes) *60).split('.')
    raSeconds = float('0.'+raSeconds) *60

    decDegrees,decMinutes = str(dec).split('.')
    decMinutes,decSeconds = str(float('0.'+decMinutes)*60).split('.')
    decSeconds = float('0.'+decSeconds)*60
    
    return (raHours,raMinutes,round(raSeconds,2)),(decDegrees,decMinutes,round(decSeconds,2))


def returnVcorr(ra,dec,velocity):
    rav = 186.7833
    decv = 12.9333
    vcorr = velocity + 300*(math.sin(math.radians(dec)) * math.sin(math.radians(decv)) + \
    math.cos(math.radians(dec))*math.cos(math.radians(decv)) * math.cos(math.radians(ra-rav)))
    return vcorr
    
    
def returnLinDiameters(major,minor,distance):
    # input major and minor in arcsec, distance in Mpc
    # outputs major and minor in kpc
    newMajor = math.tan(math.radians(float(major)))*(distance*1000)
    newMinor = math.tan(math.radians(float(minor)))*(distance*1000)
    return (newMajor,newMinor)
    
    
def returnInclination(major,minor):
    # outputs inclination in degress
    inclination = math.acos(float(minor)/float(major))*180/math.pi
    return inclination
    

def returnAngDiameters(major,minor,distance):
    # input distances in mpc, major and minor is in kpc
    # outputs angular diameters in arcsec
    newMajor = math.atan((float(major)/1000)/float(distance))*(1/3600)
    newMinor = math.atan((float(minor)/1000)/float(distance))*(1/3600)
    return (newMajor,newMinor)


def calculateAngularSeparation(ra1,dec1,ra2,dec2,dist):
    angSep = acos(cos(dec1)*cos(dec2) + sin(dec1)*sin(dec2)*cos(ra1-ra2))
    distSep = float(angSep)*float(dist)*1000
    return distSep
    

def calculateAzimuth(galaxyRA,galaxyDec,AGNra,AGNdec,galaxyDist,pa):
    # calculates the QSO-galaxy azimuth
    
    # calculate angular separations in ra and dec to determine positions on chart w.r.t. target AGN
    gRA,gDec = float(galaxyRA)*(pi/180),(90-float(galaxyDec))*(pi/180)
    agnRA,agnDec = float(AGNra)*(pi/180),(90-float(AGNdec))*(pi/180)
    
    dRA = calculateAngularSeparation(gRA,agnDec,agnRA,agnDec,galaxyDist)
    dDec = calculateAngularSeparation(agnRA,gDec,agnRA,agnDec,galaxyDist)
#     angSepAgain = calculateAngularSeparation(gRA,gDec,agnRA,agnDec,galaxyDist)
    
    if gRA >= agnRA:
        dRA = dRA
    else:
        dRA = -dRA
    if gDec >= agnDec:
        dDec = dDec
    else:
        dDec = -dDec

#     theta = atan(dRA/dDec)
#     azimuth = float(pa) - theta

    theta = atan(dRA/dDec)*180/pi
    azimuth = 180 - float(pa) - theta
    
    print 'dRA, dDec, azimuth: ', dRA,', ',dDec,', ',azimuth
    print
    return azimuth
    

###########################################################################

    
def main():
    # assuming 'theFile' contains one name per line, read the file
    

#     filename = '/Users/David/Research Documents/inclination/correlationUpdate2_first100plus.csv'
#     filename = '/Users/David/Research Documents/100largest_identified_with_galaxies2.csv'
    filename = '/Users/David/Research Documents/inclination/LG_correlation_combined2.csv'
    galaxyFilename = '/Users/David/Research Documents/NewGalaxyTable2.csv'
    theFile = open(filename,'rU')
    galaxyFile = open(galaxyFilename,'rU')
    
    reader = csv.DictReader(theFile)
    galaxyReader = csv.DictReader(galaxyFile)
    
    
    saveDirectory = '/Users/David/Research Documents/inclination/AAS_winter2014/AAS_figures/final_figures2'
    
    
    # create and write output table header and fieldnames
#     fieldnames = ('preferredName','oldName','J2000RA_Dec','alternativeNames')
#     writerOutFile = open(outname,'wt')
#     writer = createCSVTable(writerOutFile,fieldnames)
    
    
    lyaVList = []
    lyaWList = []
    lyaErrorList = []
    naList = []
    bList = []
    impactList = []
    azList = []
    incList = []
    cosIncList = []
    paList = []
    vcorrList = []
    majList = []
    difList = []
    combList = []
    envList = []
    morphList = []
    galaxyNameList = []
    
    total = 0
    totalNo = 0
    totalYes = 0
    totalIsolated = 0
    totalGroup = 0
    
       
#     associated = array()
#     ambiguous = array()
    
    for line in reader:
        lyaV = line['Lya_v']
        lyaW = line['Lya_W']
        env = line['environment']
        galaxyName = line['galaxyName']
        impact = line['impactParameter (kpc)']
        vcorr = line['vcorrGalaxy (km/s)']
        if isNumber(vcorr):
            vcorr = float(vcorr)
        maj = line['majorAxis (kpc)']
        inc = line['inclination (deg)']
        az = line['corrected_az (deg)']
#         RC3pa = line['RC3pa (deg)']
        b = line['b'].partition('pm')[0]
        include = line['include']

        if include == 'yes' or include == 'no':
            total +=1
            if isNumber(env):
                if int(env) == 0:
                    totalIsolated +=1
                if int(env) >=2:
                    totalGroup +=1
            
        if include == 'yes':
            totalYes +=1
        if include == 'no':
            totalNo +=1
        
        
        
        if isNumber(lyaV) and isNumber(b) and include == 'yes':
            b = float(b)
            
            if isNumber(az):
                az = abs(float(az))
#             elif isNumber(RC3pa):
#                 if str(RC3pa) != '-99':
#                     galaxyRA, galaxyDec = eval(line['degreesJ2000RA_DecGalaxy'])
#                     AGNra, AGNdec = eval(line['degreesJ2000RA_DecAGN'])
#                     galaxyDist = line['distGalaxy (Mpc)']
#                     if not isNumber(galaxyDist):
#                         galaxyDist = vcorr/71.0
#                     az = calculateAzimuth(galaxyRA,galaxyDec,AGNra,AGNdec,galaxyDist,RC3pa)
#                     az = abs(az)
            else:
                az = -99
            
            morph = line['morphology']
        
            localType = 'x'
            morphology = morph.lower()
            m = morphology[:3]
            if morphology != 'x':
                if bfind(m,'s'):
                    if not bfind(m,'s0'):
                        # straight spiral type
                        localType = 's'
                    else:
                        # lenticular or S0 type
                        localType = 'e'
                elif bfind(m,'e') or bfind(m,'dwarf') or bfind(m,'pec'):
                    localType = 'e'
                elif bfind(m,'len'):
                    localType = 'e'
            
                elif bfind(m,'Ir') or bfind(m,'Im') or bfind(m,'I ') or bfind(m,'IA'):
                    localType = 's'
            
                else:
                    localType = 'x' 
        
            pa = line['positionAngle (deg)']
#             RC3pa = line['RC3pa (deg)']
            if isNumber(pa):
                finalPA = float(pa)
#             elif isNumber(RC3pa):
#                 if float(RC3pa) != -99:
#                     finalPA = float(RC3pa)
            else:
                finalPA = -99
            
            if not isNumber(inc):
                inc = -99
                
            if not isNumber(maj):
                maj = -99
                
            na = str(line['Na']).split()[0]
            print 'na: ',na
            na = float(eval(na.replace('E','e')))
            i = lyaW.find('pm')
            lyaW2 = float(lyaW[:i])
            lyaErr = float(str(lyaW.split('pm')[1]))
            
            if galaxyName !='x' and isNumber(env) and isNumber(lyaV) and isNumber(na) and inc !=-99:
                lyaV = float(lyaV)
                lyaVList.append(lyaV)
                lyaWList.append(lyaW2)
                lyaErrorList.append(lyaErr)
                naList.append(float(na))
                bList.append(float(b))
                impactList.append(float(impact))
                azList.append(float(az))
                incList.append(float(inc))
#                 print 'inc: ', inc
                cosIncList.append(cos(pi/180 *float(inc)))
                paList.append(float(finalPA))
                vcorrList.append(float(vcorr))
                majList.append(float(maj))
                difList.append(vcorr - lyaV)
                envList.append(float(env))
                morphList.append(localType)
                galaxyNameList.append(galaxyName)
                print 'galaxyName: ',galaxyName


#             elif galaxyName == 'x' and isNumber(na):
#               lyaV = float(lyaV)
#                 lyaVList.append(lyaV)
#                 lyaWList.append(lyaW2)
#                 naList.append(float(na))
#                 bList.append(float(b))
#                 envList.append(float(env))

    # do some stats:
    incList2 = []
    azList2 = []
    for i,a in zip(incList,azList):
        if i>=50:
            incList2.append(i)
        if a>50:
            azList2.append(a)
    print '{0} of {1} are inclined >=50%'.format(len(incList2),len(incList))
    print '{0} of {1} have azimuth >=50%'.format(len(azList2),len(azList))
    print 'a: ',a
    print 'len: ',len(incList)
    print
    print 'total yes: ',totalYes
    print 'total no: ',totalNo
    print 'total: ',total
    print
    print 'totalIsolated: ',totalIsolated
    print 'totalGroup: ',totalGroup

        
########################################################################################
########################################################################################
########################################################################################
########################################################################################

    allInclinations = []
    allCosInclinations = []
    allPA = []
    for line in galaxyReader:
        diameters = eval(line['linDiameters (kpc)'])
        inc = line['inclination (deg)']
        pa = line['positionAngle (deg)']
        
        if isNumber(diameters[0]) and isNumber(inc):
            allCosInclinations.append(cos(pi/180 * round(float(inc),1)))
            allInclinations.append(round(float(inc),1))
            
        if isNumber(pa):
            allPA.append(float(pa))
    galaxyFile.close()

########################################################################################

    plotAzHist = False
    if plotAzHist:
        fig = figure()
        ax = fig.add_subplot(111)
        bins = [0,10,20,30,40,50,60,70,80,90]
    #     bins = [5,15,25,35,45,55,65,75,85]

    #     bins = [0,15,30,45,60,75,90]
        plot1 = hist(azList,bins=bins,histtype='bar')
        xlabel('Azimuth (deg)')
        ylabel('Number')
#         ylim(0,9)
        xlim(0,90)
        savefig('{0}/hist(azimuth)_final.pdf'.format(saveDirectory),format='pdf')
#         show()
      
########################################################################################
  
    plotPAHist = False
    if plotPAHist:
        fig = figure()
#         subplots_adjust(hspace=0.200)
        ax = fig.add_subplot(211)
        bins = [0,10,20,30,40,50,60,70,80,90]
#         fig.subplots_adjust(left=0.01, bottom=0.01, right=0.01, top=0.01, wspace=0.01, hspace=0.01)
        fig.subplots_adjust(hspace=0.4)

        plot1 = hist(paList,bins=bins,histtype='bar')
        title('Absorber-associated galaxies')
        xlabel('Position Angle (deg)')
        ylabel('Number')
#         ylim(0,5)
        
        fig.add_subplot(212)
        bins = [0,10,20,30,40,50,60,70,80,90]
    #     bins = [5,15,25,35,45,55,65,75,85]
    #     bins = [0,15,30,45,60,75,90]
        plot1 = hist(allPA,bins=bins,histtype='bar')
        title('Full Galaxy Sample')
        xlabel('Position Angle (deg)')
        ylabel('Number')
#         ylim(0,8000)
#         tight_layout()
        savefig('{0}/hist(PA)_final.pdf'.format(saveDirectory),format='pdf')
#         show()

#########################################################################################

    plotImpactHist1 = False
    if plotImpactHist1:
        fig = figure(figsize=(10,2))
#         subplots_adjust(hspace=0.200)
        ax = fig.add_subplot(111)
#         bins = [0,.10,.20,.30,.40,.50,.60,.70,.80,.90]
    #     bins = [5,15,25,35,45,55,65,75,85]
    #     bins = [0,15,30,45,60,75,90]
        plot1 = hist(array(impactList)/array(majList),bins=10,histtype='bar')
#         title('Impact Parameters')
#         xlabel('Impact Parameter (kpc)')
        ylabel('Number')
        ax.tick_params(axis='x', labelsize=0)
        ax.tick_params(axis='y',labelsize=8)
#         ylim(0,5)
#         show()
        savefig('{0}/hist(Impact)_final.pdf'.format(saveDirectory),format='pdf')

########################################################################################

    plotLyaWHist1 = False
    if plotLyaWHist1:
        fig = figure(figsize=(2,8))
        ax = fig.add_subplot(111)
        bins = [0,.10,.20,.30,.40,.50,.60,.70,.80,.90]
    #     bins = [5,15,25,35,45,55,65,75,85]
    #     bins = [0,15,30,45,60,75,90]
        plot1 = hist(lyaWList,bins=10,histtype='bar',orientation = 'horizontal')
#         title('Full Galaxy Sample')
#         xlabel('Lya W (AA)')
        xlabel('Number')
        ax.tick_params(axis='x', labelsize=8)
        ax.tick_params(axis='y',labelsize=0)
#         xlim(0,11)
#         tight_layout()
#         show()
        savefig('{0}/hist(lyaW)_final2.pdf'.format(saveDirectory),format='pdf')

########################################################################################

    plotInc_vs_dif = False
    if plotInc_vs_dif:
        fig = figure()
#         subplots_adjust(hspace=0.200)
        ax = fig.add_subplot(211)
#         bins = [0,10,20,30,40,50,60,70,80,90]
#         bins = [5,15,25,35,45,55,65,75,85,95]

    #     bins = [5,15,25,35,45,55,65,75,85]
#         bins = [0,15,30,45,60,75,90]
#         bins = [0,30,60,90]
        bins = [0,30,60,90]
        incBlue = []
        incRed = []
        for i,d,l in zip(incList,difList,lyaWList):
            if d>0:
                incBlue.append(i)
            else:
                incRed.append(i)
                
        n, bins, patches = hist(incBlue, bins)
        setp(patches, 'facecolor', 'blue', 'alpha', 0.5)               
#         plot1 = hist(incBlue,bins=bins,histtype='bar',c='blue',alpha=0.5)
        title('Blue Shifted Absorbers')
        xlabel('Inclination (deg) / W (AA)')
        ylabel('Number')
#         xlim(0,90)
        
        ax = fig.add_subplot(212)
#         bins = [5,15,25,35,45,55,65,75,85,95]
    #     bins = [5,15,25,35,45,55,65,75,85]
#         bins = [0,15,30,45,60,75,90]
    
        n, bins, patches = hist(incRed, bins)
        setp(patches, 'facecolor', 'red', 'alpha', 0.5)
#         plot2 = hist(incRed,bins=bins,histtype='bar',c='red',alpha=0.5)
        title('Red Shifted Absorbers')
        xlabel('Inclination (deg) / W (AA)')
        ylabel('Number')
#         ylim(0,1)
#         xlim(0,90)
        tight_layout()
#         savefig('/Users/David/Research Documents/inclination/hist(PA)_final.pdf',format='pdf')


        # give me the stats:
        incList2 = []
        azList2 = []
        for i,a in zip(incList,azList):
            if i>=50:
                incList2.append(i)
            if a>50:
                azList2.append(a)
        print '{0} of {1} are inclined >=50%'.format(len(incList2),len(incList))
        print '{0} of {1} have azimuth >=50%'.format(len(azList2),len(azList))
        print 'a: ',a
        print 'len: ',len(incList)
        show()

########################################################################################

    plotIncHist = False
    if plotIncHist:
        fig = figure()
#         subplots_adjust(hspace=0.200)
        ax = fig.add_subplot(211)
        bins = [0,10,20,30,40,50,60,70,80,90]
        subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.4)

    #     bins = [5,15,25,35,45,55,65,75,85]
    #     bins = [0,15,30,45,60,75,90]
        plot1 = hist(incList,bins=bins,histtype='bar')
        title('Absorber-associated galaxies')
        xlabel('Inclination (deg)')
        ylabel('Number')
#         ylim(0,5)
        
        ax = fig.add_subplot(212)
        bins = [0,10,20,30,40,50,60,70,80,90]
    #     bins = [5,15,25,35,45,55,65,75,85]
    #     bins = [0,15,30,45,60,75,90]
        plot1 = hist(allInclinations,bins=bins,histtype='bar')
        title('Full Galaxy Sample')
        xlabel('Inclination (deg)')
        ylabel('Number')
#         ylim(0,8000)
#         tight_layout()

        savefig('{0}/hist(inclination)_final.pdf'.format(saveDirectory),format='pdf')
#         show()

########################################################################################

    plotCosIncHist = False
    if plotCosIncHist:
        fig = figure()
#         subplots_adjust(hspace=0.200)
        ax = fig.add_subplot(211)
        bins = [0,.10,.20,.30,.40,.50,.60,.70,.80,.90]
        subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.4)

    #     bins = [5,15,25,35,45,55,65,75,85]
    #     bins = [0,15,30,45,60,75,90]
        plot1 = hist(cosIncList,bins=bins,histtype='bar')
        title('Absorber-associated galaxies')
        xlabel('Inclination (deg)')
        ylabel('Number')
#         ylim(0,5)

        ax = fig.add_subplot(212)
        bins = [0,.10,.20,.30,.40,.50,.60,.70,.80,.90]
    #     bins = [5,15,25,35,45,55,65,75,85]
    #     bins = [0,15,30,45,60,75,90]
        plot1 = hist(allCosInclinations,bins=bins,histtype='bar')
        title('Full Galaxy Sample')
        xlabel('Inclination (deg)')
        ylabel('Number')
#         ylim(0,8000)
#         tight_layout()
        savefig('{0}/hist(cos(inclination))_final.pdf'.format(saveDirectory),format='pdf')
#         show()

########################################################################################

    plotCosIncHist2 = True
    if plotCosIncHist2:
    
#         fig = figure(figsize=(12,10))
#         fig = figure()
#         subplots_adjust(hspace=0.200)

#         subplots_adjust(hspace=0.7)
#         ax = fig.add_subplot(311)
#         subplot(311)
        bins = [0,.10,.20,.30,.40,.50,.60,.70,.80,.90]
        blue = []
        red = []
        
        blueLya = []
        blueLyaErr = []

        redLya = []
        redLyaErr = []

        
        for d,i,l,e in zip(difList,cosIncList,lyaWList,lyaErrorList):
            if d >0:
                # blue shifted absorber, but galaxy is REDSHIFTED
                print 'd: ',d
                blue.append(i)
                blueLya.append(l)
                blueLyaErr.append(e)
            else:
                # red shifted absorber, but galaxy is BLUESHIFTED
                red.append(i)
                redLya.append(l)
                redLyaErr.append(e)
        
#         plot1 = hist([red,blue],bins=bins,histtype='bar',color=['Red','blue'],alpha=0.5)
 #        title('Red shifted aborption: Galaxies')
#         xlabel('Inclination (deg)')
#         ylim(0,6)
#         ylabel('Number')
        fig = figure(figsize=(10,2))
        ax = fig.add_subplot(111)
        subplots_adjust(top = 0.85,bottom = 0.25)
        
        hist(red,bins=bins,histtype='bar',color='red',alpha = 0.9)
        print 'average red: ',average(redLya)
        print 'median red: ', median(redLya)
        print 'average(error): ',average(redLyaErr)
        print 'median(error): ',median(redLyaErr)
        print
        title('Red shifted absorption: Galaxies')
        xlabel('Cos(inclination) = b/a')
#         ylim(0,6)
        ylabel('Number')
        savefig('{0}/hist(cos(inclination))_dif_final2_p1.pdf'.format(saveDirectory),format='pdf')


#         fig.add_subplot(212)
#         subplot(212)
        fig = figure(figsize=(10,2))
        ax = fig.add_subplot(111)
        subplots_adjust(top = 0.85,bottom=0.25)

        hist(blue,bins=bins,histtype='bar',color='Blue',alpha = 0.9)
        print 'average blue: ',average(blueLya)
        print 'median blue: ',median(blueLya)
        print 'avereage(error): ',average(blueLyaErr)
        print 'median(error): ',median(blueLyaErr)
        
        title('Blue shifted absorption: Galaxies')
        xlabel('Cos(inclination) = b/a')
#         ylim(0,5)
        ylabel('Number')
        savefig('{0}/hist(cos(inclination))_dif_final2_p2.pdf'.format(saveDirectory),format='pdf')

#         ax = fig.add_subplot(312)
#         subplot(312)    
        fig = figure(figsize=(10,2))
        ax = fig.add_subplot(111)
        subplots_adjust(top = 0.85,bottom=0.25)

        hist(allCosInclinations,bins=bins,histtype='bar',color = 'green',alpha=0.9)
        title('Full Galaxy Sample')
        xlabel('Cos(inclination) = b/a')
        ylabel('Number')
#         tight_layout()

#         savefig('{0}/hist(cos(inclination))_dif_final2_p3.pdf'.format(saveDirectory),format='pdf')
        show()

########################################################################################
########################################################################################

    plotDiameter = False
    if plotDiameter:
        fig = figure()
        ax = fig.add_subplot(111)
        countb = 0
        countr = 0
        count = -1
        labelr = 'Red Shifted Absorber'
        labelb = "Blue Shifted Absorber"
        for d,i,w in zip(difList,majList,lyaWList):
            count +=1
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
                    
            plot1 = scatter(i,w,c=color,s = 50)
        
        xlabel('Major Axis (kpc)')
        ylabel(r'Equivalent Width ($\rm m\AA$)')
        ax.grid(b=None,which='major',axis='both')
        ylim(0,1200)
        ax.legend(scatterpoints=1)
        savefig('{0}/W_vs_diameter.pdf'.format(saveDirectory),format='pdf')
#         show()


#     plotInc = False
#     if plotInc:
#         print 'incLIst: ',incList
#         fig = figure()
#         ax = fig.add_subplot(111)
#         bins = [0,10,20,30,40,50,60,70,80,90]
#         plot1 = hist(incList,bins=bins,histtype='bar')
#         xlabel('Inclination (deg)')
#         ylabel('Number')
#         ylim(0,10)
# 
#         show()
#         
#     plotCosInc = False
#     if plotCosInc:
#         fig = figure()
#         ax = fig.add_subplot(111)
#         bins = [0.0,0.10,0.20,0.30,0.40,0.50,0.60,0.70,0.80,0.90]
#         plot1 = hist(cosIncList,bins=bins,histtype='bar')
#         xlabel('Cos(inclination)')
#         ylabel('Number')
#         ylim(0,7)
# 
#         show()
        
########################################################################################
   
    plotImpact = False
    if plotImpact:
        fig = figure()
        ax = fig.add_subplot(111)
        countb = 0
        countr = 0
        count = -1
        labelr = 'Red Shifted Absorber'
        labelb = "Blue Shifted Absorber"
        for d,i,w in zip(difList,impactList,lyaWList):
            count +=1
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
                    
            plot1 = scatter(i,w,c=color,s = 50)
        
        xlabel('Impact Parameter (kpc)')
        ylabel(r'Equivalent Width ($\rm m\AA$)')
        ax.grid(b=None,which='major',axis='both')
        ylim(0,1200)
        ax.legend(scatterpoints=1)
        savefig('{0}/W_vs_impact_final2.pdf'.format(saveDirectory),format='pdf')
           
########################################################################################       
        
    plotImpact2 = False
    if plotImpact2:
        fig = figure()
        ax = fig.add_subplot(111)
        countb = 0
        countr = 0
        count = -1
        labelr = 'Red Shifted Absorber'
        labelb = "Blue Shifted Absorber"
        for d,i,w,m in zip(difList,impactList,lyaWList,majList):
            count +=1
            if d>0:
                # galaxy is behind absorber, so gas is blue shifted
                color = 'Blue'
                if countb == 0:
                    countb +=1
                    plotb = ax.scatter(i/m,w,c='Blue',s=50,label= labelb)
            if d<0:
                # gas is red shifted compared to galaxy
                color = 'Red'
                if countr == 0:
                    countr +=1
                    plotr = ax.scatter(i/m,w,c='Red',s=50,label= labelr)
                    
            plot1 = scatter(i/m,w,c=color,s = 50)
            
        # make the legend work properly
#         labelr = 'Red Shifted Absorber'
#         labelb = "Blue Shifted Absorber"
#         plotb = scatter(i[countb]/m[countb],w[countb],c='Blue',s=50,label= labelb)
#         plotr = scatter(i[countr]/m[countr],w[countr],c='Red',s=50,label= labelr)
        
        xlabel('Impact Parameter / Diameter')
        ylabel(r'Equivalent Width ($\rm m\AA$)')
        ax.grid(b=None,which='major',axis='both')
        ylim(0,1200)
        xlim(-1,70)

        ax.legend(scatterpoints=1)
        savefig('{0}/W_vs_impact-diam_final2.pdf'.format(saveDirectory),format='pdf')
#         show()
        
########################################################################################

    plotNaV = False
    if plotNaV:
        fig = figure()
        ax = fig.add_subplot(111)
        countb = 0
        countr = 0
        count = -1
        labelr = 'Red Shifted Absorber'
        labelb = "Blue Shifted Absorber"
        for d,i,n in zip(difList,impactList,naList):
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
            
        xlabel('Impact Parameter (kpc)')
        ylabel(r'Na(v) (cm$^{\rm -2}$)')
        legend(scatterpoints=1)
        ax.grid(b=None,which='major',axis='both')
        ylim(-1,5e14)
        xlim(0,500)
        savefig('{0}/NaV_vs_impact_final.pdf'.format(saveDirectory),format='pdf')
#         show()       
        
########################################################################################        

    plotNaV2 = False
    if plotNaV2:
        fig = figure()
        ax = fig.add_subplot(111)
        countb = 0
        countr = 0
        count = -1
        labelr = 'Red Shifted Absorber'
        labelb = "Blue Shifted Absorber"
        for d,i,n,m in zip(difList,impactList,naList,majList):
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
            
        xlabel(r'Impact Parameter / Diameter')
        ylabel(r'Na(v) (cm$^{\rm -2}$)')
        legend(scatterpoints=1)
        ax.grid(b=None,which='major',axis='both')
        ylim(-1,5e14)
        xlim(-1,70)
        savefig('{0}/NaV_vs_impact-diam_final.pdf'.format(saveDirectory),format='pdf')
#         show()
             
#########################################################################################
             
    plotAz = False
    if plotAz:
        fig = figure()
        ax = fig.add_subplot(111)
        countb = 0
        countr = 0
        count = -1
        labelr = 'Red Shifted Absorber'
        labelb = "Blue Shifted Absorber"
        for d,a,w,m in zip(difList,azList,lyaWList,majList):
            if d>0:
                # galaxy is behind absorber, so gas is blue shifted
                color = 'Blue'
                if countb == 0:
                    countb +=1
                    plotb = ax.scatter(a,w,c='Blue',s=50,label= labelb)
            if d<0:
                # gas is red shifted compared to galaxy
                color = 'Red'
                if countr == 0:
                    countr +=1
                    plotr = ax.scatter(a,w,c='Red',s=50,label= labelr)
                
            plot1 = scatter(a,w,c=color,s=50)
            
        xlabel(r'Azimuth (deg)')
        ylabel(r'Equivalent Width ($\rm m\AA$)')
        legend(scatterpoints=1)
        ax.grid(b=None,which='major',axis='both')
        ylim(0,1200)
        xlim(0,90)

        savefig('{0}/W_vs_azimuth_final.pdf'.format(saveDirectory),format='pdf')
#         show()

########################################################################################

    # plot azimuth normalized by galaxy size
    plotAz2 = False
    if plotAz2:
        fig = figure()
        ax = fig.add_subplot(111)
        countb = 0
        countr = 0
        count = -1
        labelr = 'Red Shifted Absorber'
        labelb = "Blue Shifted Absorber"
        
        # give some stats:
        lessThan45 = 0
        for a in azList:
            if a <=45:
                lessThan45 +=1
        print '{0}/{1} have az <= 45 degrees'.format(lessThan45,len(azList))
        print 'average, median azimuth: ',average(azList),', ',median(azList)
        
        for d,a,w,m in zip(difList,azList,lyaWList,majList):
            if d>0:
                # galaxy is behind absorber, so gas is blue shifted
                color = 'Blue'
                if countb == 0:
                    countb +=1
                    plotb = ax.scatter(a/m,w,c='Blue',s=50,label= labelb)
            if d<0:
                # gas is red shifted compared to galaxy
                color = 'Red'
                if countr == 0:
                    countr +=1
                    plotr = ax.scatter(a/m,w,c='Red',s=50,label= labelr)
                
            plot1 = scatter(a/m,w,c=color,s=50)
            
        xlabel(r'Azimuth / Major Axis')
        ylabel(r'Equivalent Width ($\rm \AA$)')
        legend(scatterpoints=1)
        ax.grid(b=None,which='major',axis='both')
#         ylim(0,600)
#         xlim(0,90)

#         savefig('/Users/David/Research Documents/inclination/W_vs_az-diameter_final.pdf',format='pdf')
#         show()

########################################################################################

    plotInc = False
    if plotInc:
        fig = figure()
        ax = fig.add_subplot(111)
        countb = 0
        countr = 0
        count = -1
        labelr = 'Red Shifted Absorber'
        labelb = "Blue Shifted Absorber"
        for d,a,w,m in zip(difList,incList,lyaWList,majList):
            if d>0:
                # galaxy is behind absorber, so gas is blue shifted
                color = 'Blue'
                if countb == 0:
                    countb +=1
                    plotb = ax.scatter(a,w,c='Blue',s=50,label= labelb)
            if d<0:
                # gas is red shifted compared to galaxy
                color = 'Red'
                if countr == 0:
                    countr +=1
                    plotr = ax.scatter(a,w,c='Red',s=50,label= labelr)
                
            plot1 = scatter(a,w,c=color,s=50)
            
        xlabel(r'Inclination (deg)')
        ylabel(r'Equivalent Width ($\rm m\AA$)')
        legend(scatterpoints=1)
        ax.grid(b=None,which='major',axis='both')
        ylim(-1,1200)
        xlim(0,90)
        savefig('{0}/W_vs_inclination_final2.pdf'.format(saveDirectory),format='pdf')
#         show()

########################################################################################

    plotCosInc = False
    if plotCosInc:
        fig = figure()
        ax = fig.add_subplot(111)
        countb = 0
        countr = 0
        count = -1
        labelr = 'Red Shifted Absorber'
        labelb = "Blue Shifted Absorber"
        for d,a,w,m in zip(difList,cosIncList,lyaWList,majList):
            if d>0:
                # galaxy is behind absorber, so gas is blue shifted
                color = 'Blue'
                if countb == 0:
                    countb +=1
                    plotb = ax.scatter(a,w,c='Blue',s=50,label= labelb)
            if d<0:
                # gas is red shifted compared to galaxy
                color = 'Red'
                if countr == 0:
                    countr +=1
                    plotr = ax.scatter(a,w,c='Red',s=50,label= labelr)
                
            plot1 = scatter(a,w,c=color,s=50)
            
        xlabel(r'Cos(inclination) = b/a')
        ylabel(r'Equivalent Width ($\rm m\AA$)')
        legend(scatterpoints=1)
        ax.grid(b=None,which='major',axis='both')
        ylim(-1,1200)
        xlim(0,1)
        savefig('{0}/W_vs_cos(inclination)_final2.pdf'.format(saveDirectory),format='pdf')
#         show()

########################################################################################

    plotCosInc2 = False
    if plotCosInc2:
        # colormap the velocity difference of the absorber
        averaged =[]
#         for i in lyaWList:
#             pass

        blueMap = cm.Blues
        redMap = cm.Reds
        
        fig = figure()
        ax = fig.add_subplot(111)
        countb = 0
        countr = 0
        count = -1
#         labelr = 'Red Shifted Absorber'
#         labelb = "Blue Shifted Absorber"
        rdif = []
        rcosInc = []
        rlyaW = []
        rMaj = []
        bdif = []
        bcosInc = []
        blyaW = []
        bMaj = []
        for d,i,w,m in zip(difList,cosIncList,lyaWList,majList):
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
        
        plotr = ax.scatter(rcosInc, rlyaW, cmap=redMap, c=rdif, s=50, vmin=0, vmax=400)
        norm = matplotlib.colors.Normalize(vmin = 0, vmax = 400)
        cbarRed = plt.colorbar(plotr,cmap=redMap,orientation='horizontal')
        cbarRed.set_label('vcorr - absorber velocity (km/s)')
        
        plotb = ax.scatter(bcosInc, blyaW, cmap=blueMap, c=bdif, s=50, vmin=0, vmax=400)
        cbarBlue = plt.colorbar(plotb,cmap=blueMap,orientation='vertical')
        cbarBlue.set_label('vcorr - absorber velocity (km/s)')
#       cbar.ax.set_yticklabels(ticks)
                    
        xlabel(r'Cos(inclination) = b/a')
        ylabel(r'Equivalent Width ($\rm \AA$)')
#         legend(scatterpoints=1)
        ax.grid(b=None,which='major',axis='both')
#         ylim(0,600)
#         savefig('/Users/David/Research Documents/inclination/W_vs_cos(inclination)_final.pdf',format='pdf')
        show()
        
        
        
        fig = figure()
        ax = fig.add_subplot(211)
        hist(rdif,color='red')
        
        ax = fig.add_subplot(212)
        hist(bdif,color='blue')
        show()
        

########################################################################################

    plotPA = False
    if plotPA:
        fig = figure()
        ax = fig.add_subplot(111)
        countb = 0
        countr = 0
        count = -1
        labelr = 'Red Shifted Absorber'
        labelb = "Blue Shifted Absorber"
        for d,a,w,m in zip(difList,paList,lyaWList,majList):
            if d>0:
                # galaxy is behind absorber, so gas is blue shifted
                color = 'Blue'
                if countb == 0:
                    countb +=1
                    plotb = ax.scatter(a,w,c='Blue',s=50,label= labelb)
            if d<0:
                # gas is red shifted compared to galaxy
                color = 'Red'
                if countr == 0:
                    countr +=1
                    plotr = ax.scatter(a,w,c='Red',s=50,label= labelr)
                
            plot1 = scatter(a,w,c=color,s=50)
            
        xlabel(r'Position Angle (deg)')
        ylabel(r'Equivalent Width ($\rm m\AA$)')
        legend(scatterpoints=1)
        ax.grid(b=None,which='major',axis='both')
        ylim(-1,1200)
        xlim(0,180)
        savefig('{0}/W_vs_positionAngle_final.pdf'.format(saveDirectory),format='pdf')
#         show()



    # the random data
#     x = array(azList)
#     y = array(lyaWList)
    
#     fig = figure()
#     ax = fig.add_subplot(211)    
#     plot1 = hist(x,bins=25,histtype='bar')
#     show()
    
#     print 'x:,',x
#     print 'y: ',y
#     print 'type: ',type(x)
# 
#     nullfmt   = NullFormatter()         # no labels
# 
#     # definitions for the axes 
#     left, width = 0.1, 0.65
#     bottom, height = 0.1, 0.65
#     bottom_h = left_h = left+width+0.02
# 
#     rect_scatter = [left, bottom, width, height]
#     rect_histx = [left, bottom_h, width, 0.2]
#     rect_histy = [left_h, bottom, 0.2, height]
# 
#     # start with a rectangular Figure
#     plt.figure(1, figsize=(8,8))
# 
#     axScatter = plt.axes(rect_scatter)
#     axHistx = plt.axes(rect_histx)
#     axHisty = plt.axes(rect_histy)
# 
#     # no labels
#     axHistx.xaxis.set_major_formatter(nullfmt)
#     axHisty.yaxis.set_major_formatter(nullfmt)
# 
#     # the scatter plot:
#     axScatter.scatter(x, y)
# 
#     # now determine nice limits by hand:
#     binwidth = 0.25
#     xmax = np.max(np.fabs(x))
#     ymax = np.max(np.fabs(y))
#     print 'ymax: ',ymax
#     print 'xmax: ',xmax
#     xymax = np.max( [np.max(np.fabs(x)), np.max(np.fabs(y))] )
# 
# #     print 'xymax',xymax
#     lim = ( int(xymax/binwidth) + 1) * binwidth
# 
#     axScatter.set_xlim( (0, xlim) )
#     axScatter.set_ylim( (0, ylim) )
# 
#     bins = arange(0, lim + binwidth, binwidth)
#     axHistx.hist(x, bins=bins)
#     axHisty.hist(y, bins=bins, orientation='horizontal')
# 
#     axHistx.set_xlim( axScatter.get_xlim() )
#     axHisty.set_ylim( axScatter.get_ylim() )
# 
#     plt.show()
               
    theFile.close()
    
###############################################################################

if __name__=="__main__":
    # do the work
    main()