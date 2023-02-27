#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('/localscratch/yadavmah/mog2/github/RNA_FUS/analysis/fig_format')
import _LondonUnderground_Colours
import defaults
defaults.set_plot_defaults()
import gsd, gsd.hoomd

with gsd.hoomd.open("../output_files/test_1/FUS_dump.gsd", 'rb') as trajectoryFile:
    print(len(trajectoryFile))
    snap = trajectoryFile[0:10]
    print(len(snap))
logFile = np.loadtxt("../output_files/test_1/energies.log")
time = logFile[:, 0]
temperature = logFile[:, 1]
potential_energy = logFile[:, 2]
kinetic_energy = logFile[:, 3]

# numberOfParticles = snap.particles.N
position = snap.particles.position[:]

x = position[:, 0]
y = position[:, 1]
z = position[:, 2]

# print(f"The shape of position:::{position.shape}")
# print(f"The total number of particles are:::{numberOfParticles}")
def plot_potential_energy():

    # plt.xscale('log', nonpositive='clip')
    # plt.yscale('log', nonpositive='clip')

    plt.plot(x, potential_energy, c="#FF0000", label='')

    plt.xlabel(r"$ x\; [units]$")
    plt.ylabel(r"$Potential \; energy \; [units]$")
    # lgd = r"$$"
    plt.autoscale(enable=True, axis='both', tight=True)
    plt.savefig("/localscratch/yadavmah/mog2/github/RNA_FUS/figures/x_potentialEne.pdf")
    print("#########, \nThe program completes without any error")
    plt.show()

#plot_potential_energy()
