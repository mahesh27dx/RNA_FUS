#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import datetime
import sys,os
import math
import numpy as np
import gsd, gsd.hoomd, gsd.pygsd
import hoomd, hoomd.md
try:
    import azplugins
except:
    from hoomd import azplugins
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import hoomd_util as hu
import pathlib
import argparse

valid_input_formats = ['.gsd']
valid_output_formats = ['.gsd']

parser = argparse.ArgumentParser(
    description=
    "Initialize the system from input file and perform MD runs")
# parser.add_argument(
#     'ifile',
#     type= str,
#     help = f'Input file (allowed formats: {valid_input_formats})'
# )
# parser.add_argument(
#     'ofile',
#     type = str
#     help = f'Output file (allowed formats: {valid_output_formats})'
# )
# ifile = pathlib.Path(args.ifile)
# ofile = pathlin.Path(args.ofile)

parser.add_argument('-dt', '--dt', type=float, help='time step size')
parser.add_argument('-time' ,'--time', type=int, help='simulation run time')
parser.add_argument('-temp' ,'--temp', type=int, help='temperature for the current simulation run')
parser.add_argument('-period', '--period', type=int, help='number of time steps between file dumps.')
args = parser.parse_args()

stat_file = 'input_files/stats_module_RNA.dat'
filein_FUS = 'input_files/calpha_FUS.pdb'
stat_rna = 'input_files/rna_stats.dat'

dt = args.dt
simulation_steps = args.time
# equilibration_steps = 1000
# T = args.temp * 0.00831446
T = args.temp
temp = T * 0.00831446
period = args.period

k_B = 1.38064852e-23                 # boltzmann constant [J/K]
N_A = 6.02214076e23                  # Avogadro constant [1/mol]
epsilon_in_kcal_mol = 0.2            # energy unit [kcal/mol]: 1epsilon = 0.2 kcal/mol
rna_length = int(250)

## Input parameters for all RNA bases
rna_param_dict = hu.rna_stats_from_file(stat_rna)
rna_type = list(rna_param_dict.keys())
rna_mass = []
rna_charge = []
rna_sigma = []
rna_lambda = []
for i in rna_type:
    rna_mass.append(rna_param_dict[i][0])
    rna_charge.append(rna_param_dict[i][1])
    rna_sigma.append(rna_param_dict[i][2]/10.) # divide by 10 to consvert angs-> nm
    rna_lambda.append(rna_param_dict[i][3])

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

print(rna_param_dict)
## Now we can translate the entire sequence into a sequence of numbers and
# assign to each Amino acid of the sequence its stats
prot_id, prot_mass, prot_charge, prot_sigma, prot_position = hu.aa_stats_sequence(filein_FUS, aa_param_dict)
# prot_id = np.array(prot_id)
prot_position_array = np.array(prot_position) / 10.
# prot_position_array = prot_position_array + 8
prot_length = len(prot_id)
prot_total_mass = np.sum(prot_mass)
prot_mass_arr = np.array(prot_mass)

# Add RNA chains
for i in range(rna_length):
    # prot_id.append(len(rna_type)-1)
    prot_id.append(len(aa_type)-1)
    prot_mass.append(329.2)
    prot_charge.append(-1)
# print(f"The rna_type:::{len(rna_type)-1}")
# print(f"The aa_type:::{aa_type}")
# print(len(prot_id))

# Relative positions of the constituent particles of the two rigid bodies
### First rigid body relative position and moment of inertia
rigid_1_coord_1 = np.sum(prot_position_array[284:371, 0] * prot_mass[284:371]) / np.sum(prot_mass[284:371])
rigid_1_coord_2 = np.sum(prot_position_array[284:371, 1] * prot_mass[284:371]) / np.sum(prot_mass[284:371])
rigid_1_coord_3 = np.sum(prot_position_array[284:371, 2] * prot_mass[284:371]) / np.sum(prot_mass[284:371])
relative_pos_rigid_1 = prot_position_array[284:371] - np.array([[rigid_1_coord_1, rigid_1_coord_2, rigid_1_coord_3]])

I_general_rigid_1 = hu.protein_moment_inertia(relative_pos_rigid_1, prot_mass[284:371]) # Moment of inertia
I_diag_rigid_1, E_vec_rigid_1 = np.linalg.eig(I_general_rigid_1)
I_diag_rigid_1 = I_diag_rigid_1.reshape(-1, 3)

## Second rigid body relative position and moment of inertia
rigid_2_coord_1 = np.sum(prot_position_array[421:453, 0] * prot_mass[421:453]) / np.sum(prot_mass[421:453])
rigid_2_coord_2 = np.sum(prot_position_array[421:453, 1] * prot_mass[421:453]) / np.sum(prot_mass[421:453])
rigid_2_coord_3 = np.sum(prot_position_array[421:453, 2] * prot_mass[421:453]) / np.sum(prot_mass[421:453])
relative_pos_rigid_2 = prot_position_array[421:453] - np.array([[rigid_2_coord_1, rigid_2_coord_2, rigid_2_coord_3]])

I_general_rigid_2 = hu.protein_moment_inertia(relative_pos_rigid_2, prot_mass[421:453])
I_diag_rigid_2, E_vec_rigid_2 = np.linalg.eig(I_general_rigid_2)
I_diag_rigid_2 = I_diag_rigid_2.reshape(-1, 3)

# Context initialize
hoomd.context.initialize("")

## Types for the rigid bodies
type_rigid_1 = [aa_type[prot_id[i]] for i in range(284, 371)]
type_rigid_2 = [aa_type[prot_id[i]] for i in range(421, 453)]

# Read the "starting_config.gsd"
system = hoomd.init.read_gsd('output_files/init_RNA_snap.gsd')
# system = hoomd.init.read_gsd('output_files/starting_config_RNA.gsd')

rigid = hoomd.md.constrain.rigid()

## First rigid body
rigid.set_param('R',
                types = type_rigid_1,
                positions = relative_pos_rigid_1)

## Second rigid body
rigid.set_param('Z',
                types = type_rigid_2,
                positions = relative_pos_rigid_2)

rigid.create_bodies(create=False)
# rigid.create_bodies()

# system.replicate(nx=1, ny=1, nz=1)
snapshot = system.take_snapshot()

all_group = hoomd.group.all()
center_group =hoomd.group.rigid_center()
non_rigid_group = hoomd.group.nonrigid()
moving_group = hoomd.group.union('moving_group', center_group, non_rigid_group)

# Harmonic potential between bonded particles
harmonic = hoomd.md.bond.harmonic()
harmonic.bond_coeff.set(['AA_bond', 'rigid_1', 'rigid_2'], k = 8360, r0 = 0.381)
harmonic.bond_coeff.set('NT_bond', k=8360, r0=0.5)
# Neighbourslist and exclusions
cell = hoomd.md.nlist.cell()
cell.set_params(r_buff=0.5)
cell.reset_exclusions(exclusions=['1-2', 'body'])

# aa_type = np.append(aa_type, rna_type)
print(rna_type[0])
# exit()
# Hydrophobic interactions
nb = azplugins.pair.ashbaugh(r_cut = 0, nlist = cell)
for i in aa_type:
    for j in aa_type:
        nb.pair_coeff.set(i, j, lam = (aa_param_dict[i][3] + aa_param_dict[j][3])/2.,
                        epsilon=0.8368, sigma=(aa_param_dict[i][2] + aa_param_dict[j][2])/10./2.,
                        r_cut = 2.0)
    nb.pair_coeff.set(i, 'R', lam = 0.0, epsilon = 0.0, sigma = 0.0, r_cut = 0.0)
    nb.pair_coeff.set(i, 'Z', lam = 0.0, epsilon = 0.0, sigma = 0.0, r_cut = 0.0)
    # nb.pair_coeff.set('R', i, lam = 0.0, epsilon = 0.0, sigma = 0.0, r_cut = 0.0)
    # nb.pair_coeff.set('Z', i, lam = 0.0, epsilon = 0.0, sigma = 0.0, r_cut = 0.0)
    nb.pair_coeff.set(i, rna_type[0], lam=0.0, epsilon=0.0, sigma=0.0, r_cut=0.0)
    # nb.pair_coeff.set(rna_type[0], i, lam=0.0, epsilon=0.0, sigma=0.0, r_cut=0.0)
nb.pair_coeff.set('R', 'R', lam = 0.0, epsilon = 0.0, sigma = 0.0, r_cut = 0.0)
nb.pair_coeff.set('R', 'Z', lam = 0.0, epsilon = 0.0, sigma = 0.0, r_cut = 0.0)
nb.pair_coeff.set('Z', 'Z', lam = 0.0, epsilon = 0.0, sigma = 0.0, r_cut = 0.0)
# nb.pair_coeff.set()

# Electrostatics
yukawa = hoomd.md.pair.yukawa(r_cut = 0.0, nlist = cell)
for i, atom1 in enumerate(aa_type):
    atom1 = aa_type[i]
    for j, atom2 in enumerate(aa_type):
        atom2 = aa_type[j]
        yukawa.pair_coeff.set(atom1, atom2,
                epsilon=aa_param_dict[atom1][1]*aa_param_dict[atom2][1]*1.73136,
                kappa = 1.0, r_cut = 3.5)
    yukawa.pair_coeff.set(atom1, 'R', epsilon = 0.0, kappa = 1.0, r_cut = 0.0)
    yukawa.pair_coeff.set(atom1, 'Z', epsilon = 0.0, kappa = 1.0, r_cut = 0.0)
    yukawa.pair_coeff.set(atom1, rna_type[0], epsilon = 0.0, kappa = 1.0, r_cut = 0.0)
    # yukawa.pair_coeff.set('R', atom1, epsilon = 0.0, kappa = 1.0, r_cut = 0.0)
    # yukawa.pair_coeff.set('Z', atom1, epsilon = 0.0, kappa = 1.0, r_cut = 0.0)
yukawa.pair_coeff.set('R', 'R', epsilon = 0.0, kappa = 1.0, r_cut = 0.0)
yukawa.pair_coeff.set('R', 'Z', epsilon = 0.0, kappa = 1.0, r_cut = 0.0)
yukawa.pair_coeff.set('Z', 'Z', epsilon = 0.0, kappa = 1.0, r_cut = 0.0)

hoomd.md.integrate.mode_standard(dt=dt)

# Set up integrator
ld = hoomd.md.integrate.langevin(group=moving_group, kT=temp, seed=399999)

# hoomd.update.box_resize(L=hoomd.variant.linear_interp([(0,system.box.Lx),
#                        (resize_steps-500,box_length)]),scale_particles=True)
for i,name in enumerate(aa_type):
    ld.set_gamma(name, gamma=aa_mass[i]/1000.0)
ld.set_gamma('R', gamma=aa_mass[i]/1000.0)
ld.set_gamma('Z', gamma=aa_mass[i]/1000.0)


hoomd.dump.gsd('output_files/RNA_FUS_dump' + '_' + str(T) + '_' +str(simulation_steps) +'.gsd', period=period, group=all_group, overwrite=True)
## Look at the unwrap_rigid variable for the .dcd dump
# hoomd.dump.dcd('output_files/FUS_dump' + '_' + str(T) + '.dcd', period=period, group=all_group, overwrite=True)

hoomd.analyze.log(filename='output_files/energies_rna' + '_' + str(T) + '.log',
                  quantities=['temperature', 'potential_energy', 'kinetic_energy'],
                  period=period, overwrite=True)
hoomd.analyze.log(filename='output_files/pressure_rna' + '_' + str(T) + '.log',
                  quantities=['pressure', 'pressure_xx', 'pressure_yy', 'pressure_zz',
                  'pressure_xy', 'pressure_xz', 'pressure_yz'], period=period, overwrite=True)

hoomd.run(simulation_steps)
