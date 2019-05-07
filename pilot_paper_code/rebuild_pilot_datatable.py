#!/usr/bin/env python

'''
By David French (frenchd@astro.wisc.edu)

$Id:  rebuild_pilot_datatable.py, v 1.0 07/08/2015

Take in LG_correlation_combined3.csv and trim it down to just the used sightlines and 
correlations

NOT FINISHED

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

    
def returnLinDiameters(major,minor,distance):
    # input major and minor in arcsec, distance in Mpc
    # outputs major and minor in kpc
    newMajor = math.tan(math.radians(float(major)))*(distance*1000)
    newMinor = math.tan(math.radians(float(minor)))*(distance*1000)
    return (newMajor,newMinor)
    
    

def returnAngDiameters(major,minor,distance):
    # input distances in mpc, major and minor is in kpc
    # outputs angular diameters in arcsec
    newMajor = math.atan((float(major)/1000)/float(distance))*(1/3600)
    newMinor = math.atan((float(minor)/1000)/float(distance))*(1/3600)
    return (newMajor,newMinor)
    
    

###########################################################################

    
def main():
    # assuming 'theFile' contains one name per line, read the file
    
    
    if getpass.getuser() == 'David':
        galaxyFilename = '/Users/David/Research_Documents/NewGalaxyTable5.csv'
        filename = '/Users/David/Research_Documents/inclination/LG_correlation_combined3.csv'

    elif getpass.getuser() == 'frenchd':
        galaxyFilename = '/usr/users/frenchd/gt/NewGalaxyTable5.csv'
        filename = '/usr/users/frenchd/inclination/LG_correlation_combined3.csv'

    else:
        print 'Could not determine username. Exiting.'
        sys.exit()
    
    
    theFile = open(filename,'rU')
    galaxyFile = open(galaxyFilename,'rU')
    pickleFile = open(pickleFilename,'wt')
    
    reader = csv.DictReader(theFile)
    galaxyReader = csv.DictReader(galaxyFile)
    
    
    # create and write output table header and fieldnames
#     fieldnames = ('preferredName','oldName','J2000RA_Dec','alternativeNames')
#     writerOutFile = open(outname,'wt')
#     writer = createCSVTable(writerOutFile,fieldnames)
    
    
    # overall structure: fullDict is a dictionary with all the lines and their data in it
    # separated into 'associated' and 'ambiguous' as the two keys. Associated contains
    # all the lists of data for lines associated with a galaxy. Ambiguous contains all
    # the lists of data for lines not unambiguously associated (could be many galaxies
    # or none)
    
    
    for line in reader:
        #grab all the values
#         print 'line: ',line
        AGNra_dec = eval(line['degreesJ2000RA_DecAGN'])
        galaxyRA_Dec = eval(line['degreesJ2000RA_DecGalaxy'])
        lyaV = line['Lya_v']
        lyaW = line['Lya_W']
        env = line['environment']
        galaxyName = line['galaxyName']
        impact = line['impactParameter (kpc)']
        galaxyDist = line['distGalaxy (Mpc)']
        pa = line['positionAngle (deg)']
        RC3pa = line['RC3pa (deg)']
        morph = line['morphology']
        vcorr = line['vcorrGalaxy (km/s)']
        maj = line['majorAxis (kpc)']
        min = line['minorAxis (kpc)']
        inc = line['inclination (deg)']
        az = line['corrected_az (deg)']
        b = line['b'].partition('pm')[0]
        include = line['include']


        # do some counting
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
        
        # for all lines
        if isNumber(lyaV) and isNumber(b):
        
            # split up column density estimates and errors
            na = str(line['Na']).split()[0]
            na = float(eval(na.replace('E','e')))
            
            # split up equivalent width measurements and errors
            i = lyaW.find('pm')
            lyaW2 = float(lyaW[:i])
            lyaErr = float(str(lyaW.split('pm')[1]))
            
            lyaV = float(lyaV)
            b = float(b)
            
            if isNumber(impact):
                impact = float(impact)
            
            if isNumber(vcorr):
                vcorr = float(vcorr)
                if not isNumber(galaxyDist):
                    galaxyDist = vcorr/71.0
                
                # difference in velocity between the absorption and the galaxy
                # if dif >0, then the absorber is BLUESHIFTED
                # if dif <0, then the absorber is REDSHIFTED
                dif = vcorr - float(lyaV)
            else:
                vcorr = 'x'
                dif = 'x'
        
            # azimuth
            if isNumber(az):
                az = abs(float(az))
            else:
                az = -99
                
            # position angle
            if isNumber(pa):
                pa = float(pa)
            elif isNumber(RC3pa):
                pa = float(RC3pa)
            else:
                pa = -99
            
            # now compute a new one with the working azimuth code
            if isNumber(pa) and isNumber(galaxyDist):
                newAz = calculateAzimuth(galaxyRA_Dec[0], galaxyRA_Dec[1], AGNra_dec[0], AGNra_dec[1], galaxyDist, pa)
            else:
                print 'no good: ', galaxyRA_Dec[0], galaxyRA_Dec[1], AGNra_dec[0], AGNra_dec[1], galaxyDist, pa
                newAz = -99
            
        
            # inclination
            if isNumber(inc):
                inc = float(inc)
                cosInc = cos(pi/180 *inc)
            else:
                inc = -99
                cosInc = -99
            
            # major axis diameter and fancy inclination calculation
            if isNumber(maj):
                maj = float(maj)
                if isNumber(min):
                    # minimum axis ratio: q0 = 0.13
                    q0 = 0.13
                    fancyInc = calculateFancyInclination(maj,min,q0)
                    
                    fancyCosInc = cos(pi/180 * float(fancyInc))
                else:
                    fancyInc = -99
                    fancyCosInc = -99
                    
            else:
                maj = -99
                fancyInc = -99
                fancyCosInc = -99
                
            
                
            # try to determine a morphology
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
            
            # for associated lines
#             if galaxyName !='x' and isNumber(env) and isNumber(lyaV) and isNumber(na) and inc !=-99:
            if include == 'yes':
                lyaVList.append(lyaV)
                lyaWList.append(lyaW2)
                lyaErrorList.append(lyaErr)
                naList.append(na)
                bList.append(b)
                impactList.append(impact)
                azList.append(az)
                newAzList.append(newAz)
                incList.append(inc)
                fancyIncList.append(fancyInc)
                cosIncList.append(cosInc)
                fancyCosIncList.append(fancyCosInc)
                paList.append(pa)
                vcorrList.append(vcorr)
                majList.append(maj)
                difList.append(dif)
                envList.append(float(env))
                morphList.append(localType)
                galaxyNameList.append(galaxyName)

            # for ambiguous lines that have a measured galaxy velocity
            if include == 'no' and isNumber(dif):
                lyaVListAmb.append(lyaV)
                lyaWListAmb.append(lyaW2)
                lyaErrorListAmb.append(lyaErr)
                naListAmb.append(na)
                bListAmb.append(b)
                impactListAmb.append(impact)
                azListAmb.append(az)
                newAzListAmb.append(newAz)
                incListAmb.append(inc)
                fancyIncListAmb.append(fancyInc)
                cosIncListAmb.append(cosInc)
                fancyCosIncListAmb.append(fancyCosInc)
                paListAmb.append(pa)
                vcorrListAmb.append(vcorr)
                majListAmb.append(maj)
                difListAmb.append(dif)
                envListAmb.append(float(env))
                morphListAmb.append(localType)
                galaxyNameListAmb.append(galaxyName)



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


    # populate the dictionary
    fullDict['lyaVList'] = lyaVList
    fullDict['lyaWList'] = lyaWList
    fullDict['lyaErrorList'] = lyaErrorList
    fullDict['naList'] = naList
    fullDict['bList'] = bList
    fullDict['impactList'] = impactList
    fullDict['azList'] = azList
    fullDict['newAzList'] = newAzList
    fullDict['incList'] = incList
    fullDict['fancyIncList'] = fancyIncList
    fullDict['cosIncList'] = cosIncList
    fullDict['fancyCosIncList'] = fancyCosIncList
    fullDict['paList'] = paList
    fullDict['vcorrList'] = vcorrList
    fullDict['majList'] = majList
    fullDict['difList'] = difList
    fullDict['envList'] = envList
    fullDict['morphList'] = morphList
    fullDict['galaxyNameList'] = galaxyNameList
        
    fullDict['lyaVListAmb'] = lyaVListAmb
    fullDict['lyaWListAmb'] = lyaWListAmb
    fullDict['lyaErrorListAmb'] = lyaErrorListAmb
    fullDict['naListAmb'] = naListAmb
    fullDict['bListAmb'] = bListAmb
    fullDict['impactListAmb'] = impactListAmb
    fullDict['azListAmb'] = azListAmb
    fullDict['newAzListAmb'] = newAzListAmb
    fullDict['incListAmb'] = incListAmb
    fullDict['fancyIncListAmb'] = fancyIncListAmb
    fullDict['cosIncListAmb'] = cosIncListAmb
    fullDict['fancyCosIncListAmb'] = fancyCosIncListAmb
    fullDict['paListAmb'] = paListAmb
    fullDict['vcorrListAmb'] = vcorrListAmb
    fullDict['majListAmb'] = majListAmb
    fullDict['difListAmb'] = difListAmb
    fullDict['envListAmb'] = envListAmb
    fullDict['morphListAmb'] = morphListAmb
    fullDict['galaxyNameListAmb'] = galaxyNameListAmb


    # grab all inclinations and position angles in the galaxy dataset
    allInclinations = []
    allCosInclinations = []
    allFancyInclinations = []
    allCosFancyInclinations = []
    allPA = []
    for line in galaxyReader:
        diameters = eval(line['linDiameters (kpc)'])
        inc = line['inclination (deg)']
        pa = line['positionAngle (deg)']
        
        if isNumber(diameters[0]) and isNumber(inc):
            allCosInclinations.append(cos(pi/180 * round(float(inc),1)))
            allInclinations.append(round(float(inc),1))
            
            # computer fancy inclination with minimum: q0 = 0.13
            q0 = 0.13
            if isNumber(diameters[1]):
                fancyInc = calculateFancyInclination(diameters[0], diameters[1],q0)
                cosFancyInc = cos(pi/180 * float(fancyInc))
            else:
                fancyInc = -99
                
            allFancyInclinations.append(fancyInc)
            allCosFancyInclinations.append(cosFancyInc)
            
        if isNumber(pa):
            allPA.append(float(pa))
            
    
    fullDict['allPA'] = allPA
    fullDict['allInclinations'] = allInclinations
    fullDict['allCosInclinations'] = allCosInclinations
    fullDict['allFancyInclinations'] = allFancyInclinations
    fullDict['allCosFancyInclinations'] = allCosFancyInclinations
            
    
    pickle.dump(fullDict,pickleFile)
    pickleFile.close()
    theFile.close()
    galaxyFile.close()
    
    
###############################################################################
###############################################################################
###############################################################################
###############################################################################


if __name__=="__main__":
    # do the work
    main()
    