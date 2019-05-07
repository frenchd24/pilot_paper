#!/usr/bin/env python

'''
By David French (frenchd@astro.wisc.edu)

$Id:  inclination1.py, v 1.0 07/02/2014

Plot some stuff for the 100largest initial results


'''

import sys
# sys.path.append("/usr/users/frenchd/local/lib/python2.6/site-packages/")
import os
import csv
# import string
# import warnings
# import urllib
# import numpy
from pylab import *
import atpy
from math import *
# from astropy.io.votable import parse,tree

# from vo.table import parse
# import vo.tree

# import Queue
# import threading
# import urllib2
# import time

__author__ = "David M. French - frenchd@astro.wisc.edu"
__version__="1.0"


def isNumber(s):
    try:
        float(s)
        return True
    except Exception:
        return False

    
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
    filename = '/Users/David/Research Documents/100largest_identified_yes.csv'
    theFile = open(filename,'rU')
    
    reader = csv.DictReader(theFile)
    
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
    paList = []
    vcorrList = []
    majList = []
    difList = []
    combList = []
    
    for line in reader:
        lyaV = line['Lya_v']
        lyaW = line['Lya_W']
        if isNumber(lyaV):
            print 'lyaV: ',lyaV
            
            lyaVList.append(lyaV)
            
            # remove the pm bit from each measurement
            i = lyaW.find('pm')
            lyaW2 = lyaW[:i]
            lyaError = lyaW[i+2:]
            
            print 'lyaW: ',lyaW2
            print 'lyaError: ',lyaError
            
            # populate lists
            lyaVList.append(float(lyaV))
            lyaWList.append(float(lyaW2))
            lyaErrorList.append(float(lyaError))
            
            impactP = float(line['impactParameter (kpc)'])
            b = float(line['b'].partition('pm')[0])
            bList.append(b)
            impactList.append(impactP)
            
            az = line['azimuth (deg)']
            RC3pa = line['RC3pa (deg)']
            if isNumber(az):
                azList.append(abs(float(az)))
            elif isNumber(RC3pa):
                if str(RC3pa) != '-99':
                    galaxyRA, galaxyDec = eval(line['degreesJ2000RA_DecGalaxy'])
                    AGNra, AGNdec = eval(line['degreesJ2000RA_DecAGN'])
                    galaxyDist = line['distGalaxy (Mpc)']
                    az = calculateAzimuth(galaxyRA,galaxyDec,AGNra,AGNdec,galaxyDist,RC3pa)
                    azList.append(abs(az))
            else:
                azList.append(-99)
                
            inc = line['inclination (deg)']
            if isNumber(inc):
                incList.append(int(inc))
            else:
                incList.append(-99)
            
            pa = line['positionAngle (deg)']
            RC3pa = line['RC3pa (deg)']
            if isNumber(pa):
                finalPA = pa
            elif isNumber(RC3pa):
                if float(RC3pa) != -99:
                    finalPA = RC3pa
                    
            else:
                finalPA = -99
                
            paList.append(finalPA)
            vcorr = float(line['vcorrGalaxy (km/s)'])
            dif = vcorr - float(lyaV)
            difList.append(dif)
            vcorrList.append(vcorr)
                        
            maj = float(line['majorAxis (kpc)'])
            majList.append(maj)
            
            comb = maj*impactP
            combList.append(comb)
            
            na = str(line['Na']).split()[0]
            print 'na: ',na
            na = eval(na.replace('E','e'))
            naList.append(na)
        
    fig = figure()
    ax = fig.add_subplot(111)
    
#     plot1 = ax.scatter(bList,lyaWList,lw=0)
#     xlabel('impact Parameter (kpc)')
#     ylabel('W')

#     plot1 = ax.scatter(azList,lyaWList,lw=0)
#     xlabel('azimuth (deg)')
#     ylabel('W')
#     xlim(-10,180)
#     ylim(0,1000)

#     plot1 = ax.scatter(bList,difList,lw=0)
#     xlabel('impact Parameter (kpc)')
#     ylabel('galaxy velocity - absorber velocity')
#     xlim(-10,180)
#     ylim(0,1000)

#     print 'incList: ',incList
#     plot1 = ax.scatter(incList,lyaWList,lw=0)
#     xlabel('inclination (deg)')
#     ylabel('W')
    
#     plot1 = ax.scatter(impactList,,lw=0)
#     xlabel('comb')
#     ylabel('na')
#     ylim(0,1000)

    plot1 = hist(azList,bins=15)

    show()

            
#         objectInfoList = [bestName,oldName,ra_dec,altnames]
#         row = dict((f,o) for f,o in zip(fieldnames,objectInfoList))
#         writer.writerow(row)
#         print
#         print
               
    
    theFile.close()
#     writerOutFile.close()
    
###############################################################################

if __name__=="__main__":
    # do the work
    main()