#!/usr/bin/env python

'''
By David French (frenchd@astro.wisc.edu)

$Id:  inclination2.py, v 2.0 07/02/2014

Plot some stuff for the 100largest initial results

v2: update for use with new table: 100largest_identified_with_galaxies.csv

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
    filename = '/Users/David/Research Documents/100largest_identified_with_galaxies2.csv'
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
    envList = []
    morphList = []
    galaxyNameList = []
    cosIncList= []
    
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
        RC3pa = line['RC3pa (deg)']
        b = line['b'].partition('pm')[0]

        
        if isNumber(lyaV) and isNumber(b):
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
            elif isNumber(RC3pa):
                if float(RC3pa) != -99:
                    finalPA = float(RC3pa)
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
            
            if galaxyName !='x' and isNumber(env) and isNumber(lyaV) and isNumber(na) and inc !=-99:
                lyaV = float(lyaV)
                lyaVList.append(lyaV)
                lyaWList.append(lyaW2)
                naList.append(float(na))
                bList.append(float(b))
                impactList.append(float(impact))
                azList.append(float(az))
                incList.append(float(inc))
                paList.append(float(finalPA))
                vcorrList.append(float(vcorr))
                majList.append(float(maj))
                difList.append(vcorr - lyaV)
                envList.append(float(env))
                morphList.append(localType)
                galaxyNameList.append(galaxyName)
                cosIncList.append(cos(pi/180 * float(inc)))
                
#             elif galaxyName == 'x' and isNumber(na):
#               lyaV = float(lyaV)
#                 lyaVList.append(lyaV)
#                 lyaWList.append(lyaW2)
#                 naList.append(float(na))
#                 bList.append(float(b))
#                 envList.append(float(env))
#               
        
    fig = figure()
    ax = fig.add_subplot(111)    
    i = 0
    for t in morphList:
        
        if t == 's':
            # spiral
            color = 'Blue'
        elif t == 'e':
            # elliptical
            color = 'Red'
        else:
            # unknown
            color = 'Green'
        
        if envList[i] >=2:
            marker = '^'
        else:
            marker = '.'
            
        newAzList = []
        for a,g in zip(azList,galaxyNameList):
        	if a == -99:
        		a2 = a
        	elif a >=0 and a <=90:
        		a2 =a
        	elif a >90 and a<=180:
        		if a-90 < 180-a:
        			# closer to 90
        				a2 = 180 - a
        		else:
        			#closer to 180 -> reset
        			a2 = 180-a
        	elif a >180 and a <=270:
        		if a-180 < 270-a:
        			# closer to 180
        			a2 = 270 - a
        		else:
        			a2 = 270-a
        	elif a >270 and a<= 360:
        		if a-270 < 360-a:
        			a2 = 360-a
        		else:
        			a2 = 360-a
        	else:
        		print '?: ',a
        		a2 = -99
        	
        	print g,' : ',a, ', a2 = ',a2
        	newAzList.append(a2)
            
#         plot1 = ax.scatter(majList[i],lyaWList[i],c=color,lw=0,marker=marker,s=60)
#         xlabel('major Axis (kpc)')
#         ylabel('Lya W (AA)')
#         xlim(0,70)
# 
#         plot1 = ax.scatter(bList[i],lyaWList[i],lw=0,c=color,marker=marker,s=60)
#         xlabel('b parameter')
#         ylabel('Lya W (AA)')
# 
        plot1 = ax.scatter(azList[i],lyaWList[i],lw=0,c=color,marker=marker,s=60)
        xlabel('azimuth (deg)')
        ylabel('Lya W (AA)')
#         xlim(-1,91)
        ylim(0,500)

#         plot1 = ax.scatter(paList[i],lyaWList[i],lw=0,c=color,marker=marker,s=60)
#         xlabel('PA (deg)')
#         ylabel('Lya W (AA)')
#         xlim(-1,91)
#         ylim(0,500)
# 
#         plot1 = ax.scatter(incList[i],lyaWList[i],lw=0,c=color,marker=marker,s=60)
#         xlabel('inclination (deg)')
#         ylabel('Lya W (AA)')
#         xlim(0,90)

#         plot1 = ax.scatter(cosIncList[i],lyaWList[i],lw=0,c=color,marker=marker,s=60)
#         xlabel('cos(inclination)')
#         ylabel('Lya W (AA)')
#         xlim(0,90)
    
#         plot1 = ax.scatter(impactList[i],lyaWList[i],lw=0,c=color,marker=marker,s=60)
#         xlabel('Impact Paramter (kpc)')
#         ylabel('Lya W (AA)')
    
#         plot1 = ax.scatter(envList[i],lyaWList[i],lw=0,c=color)
#         xlabel('# galaxies in environment')
#         ylabel('Lya W (AA)')

        i +=1


    show()
               
    theFile.close()
    
###############################################################################

if __name__=="__main__":
    # do the work
    main()