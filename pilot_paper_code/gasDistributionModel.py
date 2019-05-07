 #!/usr/bin/env python

'''
By David French (frenchd@astro.wisc.edu)

$Id:  gasDistributionModel.py, v 1.0 05/06/2015

Toy model for my pilot paper results. Simulates a flow of gas into or out of a galaxy, 
at a specific velocity and angle wrt the galaxy major axis, then gives the resulting 
COS measurement for a sightline at a given azimuth and impact parameter.




'''

import sys
import os
import csv
# import string
# import warnings
# import urllib
# import numpy
from pylab import *
# import atpy
import math
import getpass
import random
from utilities import *

# import scipy.optimize as optimization
# import pickle
# import itertools


# import numpy as np
# import matplotlib.pyplot as plt
# from matplotlib.ticker import NullFormatter



def bfind(obj, st):
    # "better find". Just returns true or false instead of an index
    
    if str(obj).find(st) !=-1:
        return True
    else:
        return False




##########################################################################################
##########################################################################################


def projection(a,b):
    # projection of vector b onto vector a
    
    proj_vector = dotproduct(a,b) * a / length(a)**2
#     proj_mag = length(proj_vector)
    
    return proj_vector


def gasSphereModel(azimuth,impact,velocity,angle):
    # 


def gasVectorModel(azimuth,impact,velocity,angle):
    #
    # azimuth and impact give the position of the QSO sightline
    # velocity gives the velocity of the gas flow wrt the galaxy (can be inflow or 
    # outflow)
    # 
    # outputs the expected velocity of the measured gas
    
    # call the sightline in the x direction.
    qsoX = 1.
    qsoY = 0
    qsoZ = 0
    qsoRho = math.sqrt(qsoX**2 + qsoY**2 + qsoZ**2)
    qsoThe = math.atan(qsoY / qsoX)
#     qsoPhi = math.atan(math.sqrt(qsoX**2 + qsoY**2) / qsoZ)
    qsoPhi = math.pi/2
    
    # sightline in the x direction
    qsoVector_sphere = array([qsoRho,qsoThe,qsoPhi])
    qsoVector = array([qsoX, qsoY, qsoZ])
    
    # galaxy major axis - NOT USED
    majAxisVector_sphere = array([1,math.pi/2, math.pi/2])
    majAxisVector = array([1,0,0])
    
    # gas vector: velocity of the flow is rho, azimuth is 90 - phi, and angle is theta
    gasRho = velocity
    gasThe = math.pi/2 - angle*(math.pi/180.)
    gasPhi = math.pi/2 - azimuth*(math.pi/180.)
    
    gasX = gasRho * math.sin(gasPhi) * math.cos(gasThe)
    gasY = gasRho * math.sin(gasPhi) * math.sin(gasThe)
    gasZ = gasRho * math.cos(gasPhi)
    
    gasVector_sphere = array([gasRho, gasThe, gasPhi])
    gasVector = array([gasX, gasY, gasZ])
    
    # projection of gasVector onto qsoVector
    proj = projection(qsoVector, gasVector)
    
    # return the projection of the gas vector onto the qso sightline vector
    return proj



def gasDistributionModel(azimuth, impact, velocity, opAngle):
    # uses gasVectorModel to build up a cone of gas of opening angle = opAngle
    # i.e. make a bunch of gasVectorModel's in a cone and call it a coherent flow
    
    # return what?

    pass
    

def test_entry():
    # test the gasVectorModel
    
    ans = 'x'
    ans = raw_input("Enter azimuth, velocity, angle: ")
    while ans != 'q':
        ans = raw_input("Enter azimuth, velocity, angle: ")
        
        if ans == 'q':
            sys.exit()
            
        else:
            try:
                az, vel, ang = eval(ans)
            except Exception,e:
                sys.stdout.write(e)
                az = 'x'
                vel = 'x'
                ang = 'x'
            
        if isNumber(az):
            
            sys.stdout.flush()
            sys.stdout.write("Azimuth = {0}, velocity = {1}, angle = {2}\r".format(az,vel,ang))
            print
            
            
            proj = gasVectorModel(az*math.pi/180,1,vel,ang*math.pi/180)
            
            # calculate magnitude with direction in x-axis (- for neg, + for pos x direction)
            proj_mag = sign(proj[0]) * length(proj)
            sys.stdout.write("Projection: {0}\r".format(proj_mag))
            print


def test_distribution():
    # test gasVectorModel by feeding it a cone distribution and plotting the results
    
    # opening angle of cone
    openingAngle = 65
    
    # number of vectors
    density = 20
    
    # angle above the plane of the major axis
    azimuth = 80
    
    # angle of center of cone wrt the x-axis
    angle = 0
    
    # velocity of wind
    velocity = 200
        
#     for ang in np.linspace(angle-openingAngle, angle+openingAngle, num=density):
#         for az in np.linspace(azimuth-openingAngle, azimuth+openingAngle, num = density):
#             
#             # feed this vector into the model
#             projection = gasVectorModel(az,1,velocity,ang)
#             
#             # add the results to the list
#             projDist.append(projection)

    angleDist = []
    projDist = []
    for ang in np.linspace(angle-openingAngle, angle+openingAngle, num=density):
            
        # feed this vector into the model
        proj = gasVectorModel(azimuth,1,velocity,ang)

        # calculate magnitude with direction in x-axis (- for neg, + for pos x direction)
        proj_mag = sign(proj[0]) * length(proj)
        
        # add the results to the list
        projDist.append(proj_mag)
        angleDist.append(ang)
            
    fig = figure()
    ax = fig.add_subplot(111)
#     bins = [0,10,20,30,40,50,60,70,80,90]


#     plot1 = hist(projDist,bins=10,histtype='bar')
    plot1 = scatter(angleDist,projDist)
    title("Distribution of gasVectorModel results")
    xlabel('Angle (centered around 90deg, 10deg opening angle')
    ylabel('Projection')
#     xlim(0,90)
#     ylim(0,10)
    
    show()
    

def main():
    
    # the correlation table filename and path   
    if getpass.getuser() == 'David':
        correlationFilename = '/Users/David/Research_Documents/gt/NGT3-TG6_500Correlation_full_500cutoff.csv'
        galaxyFilename = '/Users/David/Research_Documents/gt/NewGalaxyTable5.csv'
        
    elif getpass.getuser() == 'frenchd':
        correlationFilename = '/usr/users/frenchd/gt/NGT3-TG6_500Correlation_full_500cutoff.csv'
        galaxyFilename = '/usr/users/frenchd/gt/NewGalaxyTable5.csv'
    else:
        print 'Could not determine username. Exiting.'
        sys.exit()
    
    # open the correlation file
#     theFile = open(correlationFilename,'rU')
#     galaxyFile = open(galaxyFilename,'rU')
#     reader = csv.DictReader(theFile)
#     galaxyReader = csv.DictReader(galaxyFile)


#     test_entry()
    test_distribution()

#     a = array([1.,0,0])
#     b = array([-2.,1.,-1.])
#     
#     proj = projection(a,b)
#     print 'proj: ',proj
#     print 'length: ',length(proj)
#     print 'direction: ',sign(proj[0])

#     azimuth = 0
#     velocity = 200
#     ang = -10
#     
#     # feed this vector into the model
#     proj = gasVectorModel(azimuth,1,velocity,ang)
# 
#     # calculate magnitude with direction in x-axis (- for neg, + for pos x direction)
#     proj_mag = sign(proj[0]) * length(proj)
#     
#     print 'proj: ',proj
#     print 'proj_mag: ',proj_mag
    

if __name__=="__main__":
    main()
    
    