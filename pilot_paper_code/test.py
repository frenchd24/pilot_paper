#!/usr/bin/env python
"""
By David French (frenchd@astro.wisc.edu)

$Id: test.py 08/27/14

"""

from pylab import *
import sys
import os
import tempfile
import csv
import string
from math import *
import ast
# from matplotlib.patches import Ellipse
# from matplotlib.collections import EllipseCollection
# import matplotlib.cm as cm


file = open("/Users/David/Research Documents/TARGETLISTupdate_4.csv",'ru')
reader = csv.DictReader(file)
targets = {}
for line in reader:
	name = line['targetName']
	inst = line['instrument']
	type = line['type']
  
	if inst == 'C' and type == 'QSO':
		if targets.has_key(name):
			targets[name]+=1
		else:
			targets[name]=1
			
keys = targets.keys()
for i in sort(keys):
	print i

print 'number of targets: ',len(keys)

file.close()
