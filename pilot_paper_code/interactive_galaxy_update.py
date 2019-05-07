#!/usr/bin/env python

'''
By David M. French (frenchd@astro.wisc.edu)

$Id: interactive_galaxy_update.py, v 1.0 09/21/2016

Enter an angular size and PA and this calculates and returns linear size, azimuth,
inclination and fancyInclination


'''

import sys
import os
import csv
# import string
import warnings
import numpy
# import atpy
import getpass
from utilities import *
import math


# from astropy.io.votable import parse,tree

###########################################################################

def main():    
    
    # ask for the info

    quit = False
    while not quit:

        major = raw_input("Major axis (arcsec): ")
        minor = raw_input("Minor axis (arcsec): ")
        dist = raw_input("Dist (Mpc): ")
        pa = raw_input("PA (deg): ")
        
        if major == 'q' or minor =='q' or dist =='q' or pa=='q':
            quit = True
            break
        
        else:
            
    
        # search the table for this name
        for line in tableReader:
            preferredName = line['preferredName'].upper()
            oldName = line['oldName'].upper()
            alternatives = eval(line['alternativeNames'])
            found = False
            if response == preferredName or response == oldName:
#                 typeResponse = 'n'
                typeResponse = 'f'

                while typeResponse != 'b' and typeResponse != 'f':
                    typeResponse = raw_input("Return basic data (b) or full (f)?: ")
                    
                dummy = printOutInfo(line,typeResponse)
                found = True
                break
        
            else:
                for name in alternatives:
                    name = name.upper()
                    if name == response:
                        typeResponse = 'n'
                        while typeResponse != 'b' and typeResponse != 'f':
                            typeResponse = raw_input("Return basic data (b) or full (f)?: ")
                            
                        printOutInfo(line,typeResponse)
                        found = True
                        break

            if found:
                break
                            
        if not found:
            quitAnswer = 'x'
            while quitAnswer != 'n' and quitAnswer != 'y':
                quitAnswer = raw_input('Sorry, this galaxy was not found, search again? (y/n) ')
            if quitAnswer == 'n':
                quit = True
            else:
                quit = False
        
        if found:
            againAnswer = 'x'
            while againAnswer != 'y' and againAnswer != 'n':
                againAnswer = raw_input("Search again? (y/n) ")
            if againAnswer == 'y':
                quit = False
            if againAnswer == 'n':
                quit = True
                break
            
 
    theFile.close()
    print "Done."
    print
    
###############################################################################

if __name__=="__main__":
    # do the work
    main()
