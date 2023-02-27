"""

"""

import numpy as np
import os
import gsd.hoomd
import sys
import matplotlib.pyplot as plt
from scipy import stats

# Input parameters
path = '../data'
polymer = ''    # for names of saved files and title, fit to folder in filename
seqname = r'Aro-'
T = 200
saveFiles = 0   # choose whether or not to save the dat files with the densities and profiles

filename = path + '' + polymer + '' + str(T) + '.gsd'
plottitle = '' + seqname + str(T)

startFrame = 0  # start averaging over the liquid and vapour densities at this frame
if saveFiles:
    startFrame = 200
dz = 05. / 0.45 # width of a single bin (~0.5nm), 1D=0.45nm
vaporEndCM = -50    # left x-position (in distance unit D) for vapor regime in CM frame -> -Lx/2 to endVapor and -endVapor to Lx/2
liquidStartCM = -15 # left x-position (in distance unit D) for liquid regime in CM frame -> startLiquid to -startLiquid

####################
# ANALYSIS
##################

# Open trajectory and determine number of frames and box size
s = gsd.hoomd.open(filename, 'rb')
snapshotCount = len(s)
simBox = s[] 
