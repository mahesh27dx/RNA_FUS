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

# filename = path + '' + polymer + '' + str(T) + '.gsd'
# plottitle = '' + seqname + str(T)

startFrame = 0  # start averaging over the liquid and vapour densities at this frame
if saveFiles:
    startFrame = 200
dz = 05. / 0.45 # width of a single bin (~0.5nm), 1D=0.45nm
# vaporEndCM = -50    # left x-position (in distance unit D) for vapor regime in CM frame -> -Lx/2 to endVapor and -endVapor to Lx/2
# liquidStartCM = -15 # left x-position (in distance unit D) for liquid regime in CM frame -> startLiquid to -startLiquid

####################
# ANALYSIS
##################

# Open trajectory and determine number of frames and box size
filename = "../output_files/FUS_dump_300.gsd"
s = gsd.hoomd.open(filename, 'rb')
snapshotCount = len(s)
simBox = s[0].configuration.box
masses = s[0].particles.mass
screenTimer = int((snapshotCount - startFrame) / 20)
if screenTimer == 0:
    screenTimer == 1

# Determine number of sizes and size of bins
histoVolume = simBox[0] * simBox[1] * dz
histoCount = int(simBox[2] / dz) + 1
histoBins = np.linspace(-simBox[2] / 2.0, simBox[2] / 2.0, num=histoCount)

# Prepare histograms
distribMonomerPos = []
distribMonomerCMPos = []
print('number of frames: ', snapshotCount)

for i in range(startFrame, snapshotCount):
    if i % screenTimer == 0:
        print(str(i) + "/" + str(snapshotCount))

    # Compute density distribution of monomers in lab frame
    distribMonomerPos.append(np.histogram(s[i].particles.position[:, 2], bins=histoBins,
                                            weights=masses, density=False))

    # Compute density distribution of monomers in CM frame
    # (first remove CM from the box boundaries, then pos = CM)
    particlesPosShifted = s[i].particles.position[:]
    tempDistrib = np.histogram(particlesPosShifted[:, 2], bins=histoBins,
                                weights=masses, density=False)[0]
    maxDensBin = np.argmax(tempDistrib)
    particlesPosShifted[:, 2] = particlesPosShifted[:, 2] - (maxDensBin * dz - simBox[2] / 2.0)
