#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import datetime
import sys,os
import math
import numpy as np
import gsd, gsd.hoomd, gsd.pygsd
import hoomd, hoomd.md
#import azplugins
try:
    import azplugins
except:
    from hoomd import azplugins
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import hoomd_util as hu

stat_file = 'input_files/stats_module.dat'

dt = 0.001
T = 300
simulation_steps = int(10e5)
frames_output = 1000
# box_length = 700

# Input parameters for all amino acids
aa_param_dict = hu.aa_stats_from_file(stat_file)
aa_type = list(aa_param_dict.keys())
aa_mass = []
aa_charge = []
aa_sigma = []
aa_lambda = []
for i in aa_type:
    aa_mass.append(aa_param_dict[i][0])
    aa_charge.append(aa_param_dict[i][1])
    aa_sigma.append(aa_param_dict[i][2]/10.) # divide by 10 to consvert angs-> nm
    aa_lambda.append(aa_param_dict[i][3])

hoomd.context.initialize("")

# Read the "starting_config.gsd"
system = hoomd.init.read_gsd('FUS_initial_snapshot.gsd')
# system.replicate(nx=nx, ny=ny, nz=nz)
snapshot = system.take_snapshot()
all_group = hoomd.group.all()
center_group =hoomd.group.rigid_center()
non_rigid_group = hoomd.group.nonrigid()
moving_group = hoomd.group.union('moving_group', center_group, non_rigid_group)

# Harmonic potential between bonded particles
harmonic = hoomd.md.bond.harmonic()
harmonic.bond_coeff.set(['AA_bond', 'rigid_1', 'rigid_2'], k = 8360, r0 = 0.381)

# Neighbourslist and exclusions
# buffer1 = box_length / 5
# cell = hoomd.md.nlist.cell(r_buff=0.4)
cell = hoomd.md.nlist.tree()
cell.reset_exclusions(exclusions=['1-2', 'body'])

# Non-bonded Pairwise interactions for two rigid bodies
nb = azplugins.pair.ashbaugh(r_cut = 0, nlist = cell)
for i in aa_type:
    for j in aa_type:
        nb.pair_coeff.set(i, j, lam = (aa_param_dict[i][3] + aa_param_dict[j][3])/2.,
                        epsilon=0.8368, sigma=(aa_param_dict[i][2] + aa_param_dict[j][2])/10.0/2.0, r_cut = 2.0)
    nb.pair_coeff.set(i, 'R', lam = 0.0, epsilon = 0.0, sigma = 0.0, r_cut = 0.0)
    nb.pair_coeff.set(i, 'Z', lam = 0.0, epsilon = 0.0, sigma = 0.0, r_cut = 0.0)
    nb.pair_coeff.set('R', i, lam = 0.0, epsilon = 0.0, sigma = 0.0, r_cut = 0.0)
    nb.pair_coeff.set('Z', i, lam = 0.0, epsilon = 0.0, sigma = 0.0, r_cut = 0.0)
nb.pair_coeff.set('R', 'R', lam = 0.0, epsilon = 0.0, sigma = 0.0, r_cut = 0.0)
nb.pair_coeff.set('R', 'Z', lam = 0.0, epsilon = 0.0, sigma = 0.0, r_cut = 0.0)
nb.pair_coeff.set('Z', 'Z', lam = 0.0, epsilon = 0.0, sigma = 0.0, r_cut = 0.0)

# Electrostatics
yukawa = hoomd.md.pair.yukawa(r_cut = 0.0, nlist = cell)
for i, atom1 in enumerate(aa_type):
    atom1 = aa_type[i]
    for j, atom2 in enumerate(aa_type):
        atom2 = aa_type[j]
        yukawa.pair_coeff.set(atom1, atom2, epsilon=aa_param_dict[atom1][1]*aa_param_dict[atom2][1]*1.73136, kappa = 1.0, r_cut = 3.5)
    yukawa.pair_coeff.set(atom1, 'R', epsilon = 0.0, kappa = 1.0, r_cut = 0.0)
    yukawa.pair_coeff.set(atom1, 'Z', epsilon = 0.0, kappa = 1.0, r_cut = 0.0)
    yukawa.pair_coeff.set('R', atom1, epsilon = 0.0, kappa = 1.0, r_cut = 0.0)
    yukawa.pair_coeff.set('Z', atom1, epsilon = 0.0, kappa = 1.0, r_cut = 0.0)
yukawa.pair_coeff.set('R', 'R', epsilon = 0.0, kappa = 1.0, r_cut = 0.0)
yukawa.pair_coeff.set('R', 'Z', epsilon = 0.0, kappa = 1.0, r_cut = 0.0)
yukawa.pair_coeff.set('Z', 'Z', epsilon = 0.0, kappa = 1.0, r_cut = 0.0)

hoomd.md.integrate.mode_standard(dt=dt)

# Set up integrator
ld = hoomd.md.integrate.langevin(group=moving_group, kT=T, seed=39999)

# hoomd.update.box_resize(L=hoomd.variant.linear_interp([(0,system.box.Lx),
#                        (resize_steps-500,box_length)]),scale_particles=True)
# for i,name in enumerate(aa_type):
#     ld.set_gamma(name, gamma=aa_mass[i]/1000.0)

# langevin.set_gamma('R', gamma=1)
# langevin.set_gamma('Z', gamma=1)
# langevin.set_gamma('R', gamma=prot_mass_arr/1000.0)
# langevin.set_gamma('Z', gamma=prot_mass_arr/1000.0)

#print(snap.particles.position)
hoomd.dump.gsd('FUS_dump.gsd', period=simulation_steps/frames_output, group=all_group, overwrite=True, truncate=True)
# hoomd.analyze.log(filename="potential_ene.log", quantities=['potential_energy'], period=100, overwrite=False)
hoomd.run(simulation_steps)
#################################################################################################
# ## resize.gsd contains replicated system with 100 chains and box size of 15x15x15 nm
# ----------------------------------------------------------------------------------------------
# ## 3.0 Extend 15nm cube to 15x15x280nm slab
################################################################################################
# def extend(s):
#     boxdim = s.configuration.box[:3]
#     zmin,zmax,dz = -boxdim[2]/2., boxdim[2]/2., boxdim[2]
#     pos1 =  s.particles.position
#     pos = pos1.copy()
#     skip=0
#     ncomp=1
#     for k in range(ncomp):
#         nchain = int(s.particles.N/chain_length)
#         nres = chain_length
#         for i in range(nchain):
#             mol_coord = pos[i*nres+skip:(i+1)*nres+skip,2]
#             for j in range(1,nres):
#                 dist2 = (mol_coord[j] - mol_coord[j-1])**2
#                 if dist2 > 8:
#                     excess = np.sign(mol_coord[j] - mol_coord[j-1])*dz
#                     mol_coord[j] = mol_coord[j] - excess
#                 com = np.mean(mol_coord)
#                 if com < zmin:
#                     mol_coord += dz
#                 elif com > zmax:
#                     mol_coord -= dz
#             pos[i*nres+skip:(i+1)*nres+skip,2] = mol_coord
#         skip += nchain*nres
#     return pos
#
#
# f = gsd.pygsd.GSDFile(open('rigid_FUS_start.gsd','rb'))
# t = gsd.hoomd.HOOMDTrajectory(f)
# s1 = t[0]
# s = gsd.hoomd.Snapshot()
# s.particles.N = s1.particles.N
# s.particles.types = s1.particles.types
# s.particles.typeid = s1.particles.typeid
# s.particles.mass = s1.particles.mass
# s.particles.charge = s1.particles.charge
# s.particles.position = extend(s1)
# s.bonds.N = s1.bonds.N
# s.bonds.types = s1.bonds.types
# s.bonds.typeid = s1.bonds.typeid
# s.bonds.group = s1.bonds.group
# s.configuration.box = s1.configuration.box
# s.configuration.dimensions=3
# s.configuration.box = [s1.configuration.box[0],s1.configuration.box[1],slab_z_length,0,0,0]
# s.configuration.step = 0
# outfile = gsd.hoomd.open('box2slab_extend.gsd','wb')
# outfile.append(s)
# outfile.close()
# #################################################################################################
# # ### Minimized slab formed and saved in minimize.gsd
# #---------------------------------------------------------------------------------------------
# # ## 4.0. Run a production slab simulation using minimize.gsd from previous step
# ################################################################################################
# hoomd.context.initialize()
# system = hoomd.init.read_gsd('box2slab_extend.gsd')
#
# n_steps = production_steps # 1 microseconds
#
# fileroot = 'Production'
# nl = hoomd.md.nlist.cell()
#
# ## Bonds
# # Harmonic potential between bonded particles
# harmonic = hoomd.md.bond.harmonic()
# harmonic.bond_coeff.set(['AA_bond', 'rigid_1', 'rigid_2'], k = 8368, r0 = bond_length)
#
# # Non-bonded Pairwise interactions
# nl.reset_exclusions(exclusions=['1-2', 'body'])
# nb = azplugins.pair.ashbaugh(r_cut = 0, nlist = nl)
# for i in aa_type:
#     for j in aa_type:
#         nb.pair_coeff.set(i, j, lam = (aa_param_dict[i][3] + aa_param_dict[j][3])/2.,
#                         epsilon=0.8368, sigma=(aa_param_dict[i][2] + aa_param_dict[j][2])/10./2., r_cut=2.0)
#     nb.pair_coeff.set(i, 'R', lam=0., epsilon=0., sigma=0., r_cut=0.)
#     nb.pair_coeff.set(i, 'Z', lam=0., epsilon=0, sigma=0, r_cut=0.)
#     nb.pair_coeff.set('R', i, lam=0., epsilon=0., sigma=0., r_cut=0.)
#     nb.pair_coeff.set('Z', i, lam=0., epsilon=0, sigma=0, r_cut=0.)
# nb.pair_coeff.set('R', 'R', lam=0., epsilon=0, sigma=0, r_cut=0)
# nb.pair_coeff.set('R', 'Z', lam=0., epsilon=0, sigma=0, r_cut=0)
# nb.pair_coeff.set('R', 'Z', lam=0., epsilon=0, sigma=0, r_cut=0)
# nb.pair_coeff.set('Z', 'Z', lam=0., epsilon=0, sigma=0, r_cut=0)
#
# # Electrostatics
# yukawa = hoomd.md.pair.yukawa(r_cut=0.0, nlist=nl)
# # yukawa.pair_coeff.set('R','Z', epsilon=1.73136, kappa=1.0, r_cut=3.5)
# for i, atom1 in enumerate(aa_type):
#     atom1 = aa_type[i]
#     for j, atom2 in enumerate(aa_type):
#         atom2 = aa_type[j]
#         yukawa.pair_coeff.set(atom1, atom2, epsilon=aa_param_dict[atom1][1]*aa_param_dict[atom2][1]*1.73136, kappa=1.0, r_cut=3.5)
#     yukawa.pair_coeff.set(atom1, 'R', epsilon=0, kappa=1.0, r_cut=0)
#     yukawa.pair_coeff.set(atom1, 'Z', epsilon=0, kappa=1.0, r_cut=0)
#     yukawa.pair_coeff.set('R', atom1, epsilon=0, kappa=1.0, r_cut=0)
#     yukawa.pair_coeff.set('Z', atom1, epsilon=0, kappa=1.0, r_cut=0)
# yukawa.pair_coeff.set('R', 'R', epsilon=0, kappa=1.0, r_cut=0)
# yukawa.pair_coeff.set('R', 'Z', epsilon=0, kappa=1.0, r_cut=0)
# yukawa.pair_coeff.set('Z', 'Z', epsilon=0, kappa=1.0, r_cut=0)
#
# # # Grouping of the particles
# all_group = hoomd.group.all()
# center_group =hoomd.group.rigid_center()
# non_rigid_group = hoomd.group.nonrigid()
# moving_group = hoomd.group.union('moving_group', center_group, non_rigid_group)
#
# ## Set up integrator
# hoomd.md.integrate.mode_standard(dt=production_dt) # Time units in ps
# temp = production_T*0.00831446
# integrator = hoomd.md.integrate.langevin(group=moving_group, kT=temp, seed=399991) # Temp is kT/0.00831446
# #for cnt,i in enumerate(aakeys):
#     #integrator.set_gamma(i,gamma=aamass[cnt]/1000.0)
# ## Outputs
# hoomd.analyze.log(filename=fileroot+'.log', quantities=['potential_energy', 'pressure_xx', 'pressure_yy', 'pressure_zz', 'temperature','lx','ly','lz'], period=100000, overwrite=False, header_prefix='#')
# hoomd.analyze.log(filename='stress.log', quantities=['pressure_xy', 'pressure_xz', 'pressure_yz'], period=100000, overwrite=False, header_prefix='#') # Output stress tensor
# hoomd.dump.gsd('restart_tmp1.gsd', period=1000000, group=all_group, truncate=True)
# hoomd.dump.gsd('restart_tmp2.gsd', period=1000000, group=all_group, truncate=True, phase=500000)
# hoomd.dump.dcd(fileroot+'_dump.dcd', period=100000, group=all_group, overwrite=False)
#
# ## Run simulation
# hoomd.run_upto(1000, limit_hours=48)
########################################################################################################
