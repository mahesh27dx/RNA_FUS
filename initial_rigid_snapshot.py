# -*- coding: utf-8 -*-
#!/usr/bin/python3

# this script starts an MD simulation of multiple chains of a chosen sequence of the protein hnRNPA1
# short-range pairwise potential: HPS model (accounts for both protein-protein and protein-splvent interactions)

import sys,os
import numpy as np
import gsd, gsd.hoomd
import hoomd, hoomd.md
from hoomd import azplugins
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import hoomd_util as hu

# UNITS: distance -> nm   (!!!positions and sigma in files are in agstrom!!!)
#        mass -> amu
#        energy -> kJ/mol
# ### MACROs
# production_dt=0.01 # Time step for production run in picoseconds

stat_file = 'input_files/stats_module.dat'
filein_FUS = 'input_files/calpha_FUS.pdb'

if __name__=='__main__':
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
        aa_sigma.append(aa_param_dict[i][2]/10.)
        aa_lambda.append(aa_param_dict[i][3])

    # Now we can translate the entire sequence into a sequence of numbers and
    # assign to each Amino acid of the sequence its stats
    prot_id, prot_mass, prot_charge, prot_sigma, prot_position = hu.aa_stats_sequence(filein_FUS, aa_param_dict)
    prot_position_array = np.array(prot_position) # converts the position list into an array
    prot_length = len(prot_id)
    prot_total_mass = np.sum(prot_mass)
    prot_mass_arr = np.array(prot_mass)
    # poly1 = :284 ## find a clever way to assign the index
    # Relative positions of the constituent particles of the two rigid bodies
    ### First rigid body relative position and moment of inertia
    rigid_1_coord_1 = np.sum(prot_position_array[284:371, 0] * prot_mass[284:371]) / np.sum(prot_mass[284:371])
    rigid_1_coord_2 = np.sum(prot_position_array[284:371, 1] * prot_mass[284:371]) / np.sum(prot_mass[284:371])
    rigid_1_coord_3 = np.sum(prot_position_array[284:371, 2] * prot_mass[284:371]) / np.sum(prot_mass[284:371])
    relative_pos_rigid_1 = prot_position_array[284:371] - np.array([[rigid_1_coord_1, rigid_1_coord_2, rigid_1_coord_3]])
    # print(relative_pos_rigid_1.shape)
    I_general_rigid_1 = hu.protein_moment_inertia(relative_pos_rigid_1, prot_mass[284:371]) # Moment of inertia
    I_diag_rigid_1, E_vec_rigid_1 = np.linalg.eig(I_general_rigid_1)
    I_diag_rigid_1 = I_diag_rigid_1.reshape(-1, 3)

    ### Second rigid body relative position and moment of inertia
    rigid_2_coord_1 = np.sum(prot_position_array[421:453, 0] * prot_mass[421:453]) / np.sum(prot_mass[421:453])
    rigid_2_coord_2 = np.sum(prot_position_array[421:453, 1] * prot_mass[421:453]) / np.sum(prot_mass[421:453])
    rigid_2_coord_3 = np.sum(prot_position_array[421:453, 2] * prot_mass[421:453]) / np.sum(prot_mass[421:453])
    relative_pos_rigid_2 = prot_position_array[421:453] - np.array([[rigid_2_coord_1, rigid_2_coord_2, rigid_2_coord_3]])
    # print(relative_pos_rigid_2.shape)

    I_general_rigid_2 = hu.protein_moment_inertia(relative_pos_rigid_2, prot_mass[421:453]) # Moment of inertia
    I_diag_rigid_2, E_vec_rigid_2 = np.linalg.eig(I_general_rigid_2)
    I_diag_rigid_2 = I_diag_rigid_2.reshape(-1, 3)

    # Position vectors of two rigid bodies
    relative_pos_rigid = np.append([[rigid_1_coord_1], [rigid_1_coord_2], [rigid_1_coord_3]], [[rigid_2_coord_1], [rigid_2_coord_2], [rigid_2_coord_3]]).reshape(-1, 3)
    # print(relative_pos_rigid.shape)

    # Simulation parameters
    nx = 5
    nx = 5
    nx = 4
    boxsize = 15.0
    slab_z_length = 280.0
    resize_dt = 0.01 # Time step in picoseconds for box resizing
    resize_T = 300 # Temperature for resizing run in Kelvin
    production_dt = 0.01
    temp = resize_T * 8.3144598/1000.

    # Build a HOOMD snapshot of the system
    snap = gsd.hoomd.Snapshot()
    snap.particles.types = aa_type + ['R'] + ['Z']  # types of particles, ['R'] and ['Z'] are two rigid bodies

    pos_chain_1 = prot_position_array[:284].reshape(-1, 3)
    pos_chain_2 = prot_position_array[371:421].reshape(-1, 3)
    pos_chain_3 = prot_position_array[453:].reshape(-1, 3)
    rigid_1_coord = np.array([[rigid_1_coord_1, rigid_1_coord_2, rigid_1_coord_3]])
    rigid_2_coord = np.array([[rigid_2_coord_1, rigid_2_coord_2, rigid_2_coord_3]])
    position = np.concatenate((pos_chain_1, rigid_1_coord, pos_chain_2, rigid_2_coord,
                                pos_chain_3), axis=0)
    print(f"The shape of position is::{position.shape} ")

    snap.particles.position = position # positions of the free particles and rigid body constituent particles in the system
    snap.particles.N = len(position)

    mass_chain_1 = prot_mass[:284]
    mass_chain_2 = prot_mass[371:421]
    mass_chain_3 = prot_mass[453:]
    mass_rigid_1 = np.sum(prot_mass[284:371])
    mass_rigid_2 = np.sum(prot_mass[421:453])
    mass = np.concatenate((mass_chain_1, mass_rigid_1, mass_chain_2,
                            mass_rigid_2, mass_chain_3), axis=None)
    snap.particles.mass = mass

    #charge = np.empty(0, dtype=float)
    charge_chain_1 = prot_charge[:284]
    charge_chain_2 = prot_charge[371:421]
    charge_chain_3 = prot_charge[453:]
    charge_rigid_1 = np.sum(prot_charge[284:371])
    charge_rigid_2 = np.sum(prot_charge[421:453])
    charge = np.concatenate((charge_chain_1, charge_rigid_1, charge_chain_2,
                            charge_rigid_2, charge_chain_3), axis=None)
    # print(charge.shape)
    snap.particles.charge = charge

    ## Moment of inertia of the system
    MOI_chain_1 = np.zeros((len(prot_charge[:284]), 3), dtype=float)
    MOI_chain_2 = np.zeros((len(prot_charge[371:421]), 3), dtype=float)
    MOI_chain_3 = np.zeros((len(prot_charge[453:]), 3), dtype=float)
    # print(MOI_chain_1)
    moment_inertia = np.concatenate((MOI_chain_1, I_diag_rigid_1, MOI_chain_2,
                                    I_diag_rigid_2, MOI_chain_3), axis=0)
    # print(I_diag_rigid_2)
    snap.particles.moment_inertia = moment_inertia

    ## Orientation of the system
    orien_chain_1 = np.zeros((len(prot_charge[:284]), 4), dtype=float)
    orien_chain_2 = np.zeros((len(prot_charge[371:421]), 4), dtype=float)
    orien_chain_3 = np.zeros((len(prot_charge[453:]), 4), dtype=float)
    orien_rigid_1 = np.array([[0, 1, 1, 1]])
    orien_rigid_2 = np.array([[0, 1, 1, 1]])
    orientation = np.concatenate((orien_chain_1, orien_rigid_1, orien_chain_2,
                                    orien_rigid_2, orien_chain_3), axis=0)
    #print(orientation)
    snap.particles.orientation = orientation

    ## Composite body associated with each particle. Value [-1] indicates no body
    body = np.empty(0, dtype=int)
    body = np.append(body, [-1] * len(charge_chain_1) + [284] + [-1] * len(charge_chain_2)
                        + [334] + [-1] * len(charge_chain_3))
    snap.particles.body = body

    # The number of bonds between the residues
    bonds = np.empty((0, 2), dtype=int)
    for i in range(len(position) - 1):
        bonds = np.append(bonds, [[i, i+1]], axis=0)
    snap.bonds.group = bonds
    snap.bonds.N = len(bonds)
    snap.bonds.types = ['AA_bond', 'rigid_1', 'rigid_2']

    # Type ID's of the polymer chain
    id_chain_1 = prot_id[:284]
    id_chain_2 = prot_id[371:421]
    id_chain_3 = prot_id[453:]
    type_id = np.arange(0, len(position))
    prot_id = np.array(prot_id)
    type_id[:284] = prot_id[:284]
    type_id[284] = 20
    type_id[285:285+50] = prot_id[371:421]
    type_id[335] = 21
    type_id[336:] = prot_id[453:]
    print(type_id.shape,position.shape,moment_inertia.shape,orientation.shape)
    snap.particles.typeid = type_id

    bond_length = 0.38
    box_length = bond_length * prot_length + 10

    ## Dimensions of the box
    snap.configuration.dimensions = 3
    snap.configuration.box = [box_length, box_length, box_length, 0, 0, 0]
    snap.configuration.step = 0

    # Write the snapshot to the file
    with gsd.hoomd.open(name='starting_config.gsd', mode='wb') as f:
        f.append(snap)
        f.close()

    hoomd.context.initialize("--mode=cpu")

    # Read the "starting_config.gsd"
    system = hoomd.init.read_gsd('starting_config.gsd')

    snapshot = system.take_snapshot()

    # Define the rigid bodies
    ## Types for the rigid bodies
    type_rigid_1 = [aa_type[prot_id[i]] for i in range(284, 371)]
    type_rigid_2 = [aa_type[prot_id[i]] for i in range(421, 453)]

    rigid = hoomd.md.constrain.rigid()

    ### First rigid body
    rigid.set_param('R',
                    types = type_rigid_1,
                    positions = relative_pos_rigid_1[:])

    ### Second rigid body
    rigid.set_param('Z',
                    types = type_rigid_2,
                    positions = relative_pos_rigid_2[:])

    rigid.create_bodies()

    ## Add bonds
    #system.bonds.add('rigid_1', 285, 286)
    # for i in range(len(relative_pos_rigid_1)):
    #     system.bonds.add('rigid_1', 284, 528 + i - len(relative_pos_rigid_1)
    #                                             - len(relative_pos_rigid_2))
    #
    # for i in range(len(relative_pos_rigid_2) - 1):
    #     system.bonds.add('rigid_2', 335, 528 + i - len(relative_pos_rigid_2))
    #
    # # Grouping of the particles
    all_group = hoomd.group.all()
    center_group =hoomd.group.rigid_center()
    non_rigid_group = hoomd.group.nonrigid()
    moving_group = hoomd.group.union('moving_group', center_group, non_rigid_group)

    # Harmonic potential between bonded particles
    harmonic = hoomd.md.bond.harmonic()
    harmonic.bond_coeff.set(['AA_bond', 'rigid_1', 'rigid_2'], k = 8368, r0 = bond_length)

    # Neighbourslist and exclusions
    nl = hoomd.md.nlist.cell()
    nl.reset_exclusions(exclusions=['1-2', 'body'])

    # Non-bonded Pairwise interactions
    nb = azplugins.pair.ashbaugh(r_cut = 0, nlist = nl)
    for i in aa_type:
        for j in aa_type:
            nb.pair_coeff.set(i, j, lam = (aa_param_dict[i][3] + aa_param_dict[j][3])/2.,
                            epsilon=0.8368, sigma=(aa_param_dict[i][2] + aa_param_dict[j][2])/10./2., r_cut=2.0)
        nb.pair_coeff.set(i, 'R', lam=0., epsilon=0., sigma=0., r_cut=0.)
        nb.pair_coeff.set(i, 'Z', lam=0., epsilon=0, sigma=0, r_cut=0.)
        nb.pair_coeff.set('R', i, lam=0., epsilon=0., sigma=0., r_cut=0.)
        nb.pair_coeff.set('Z', i, lam=0., epsilon=0, sigma=0, r_cut=0.)
    nb.pair_coeff.set('R', 'R', lam=0., epsilon=0, sigma=0, r_cut=0)
    nb.pair_coeff.set('R', 'Z', lam=0., epsilon=0, sigma=0, r_cut=0)
    nb.pair_coeff.set('R', 'Z', lam=0., epsilon=0, sigma=0, r_cut=0)
    nb.pair_coeff.set('Z', 'Z', lam=0., epsilon=0, sigma=0, r_cut=0)

    # Electrostatics
    yukawa = hoomd.md.pair.yukawa(r_cut=0.0, nlist=nl)
    # yukawa.pair_coeff.set('R','Z', epsilon=1.73136, kappa=1.0, r_cut=3.5)
    for i, atom1 in enumerate(aa_type):
        atom1 = aa_type[i]
        for j, atom2 in enumerate(aa_type):
            atom2 = aa_type[j]
            yukawa.pair_coeff.set(atom1, atom2, epsilon=aa_param_dict[atom1][1]*aa_param_dict[atom2][1]*1.73136, kappa=1.0, r_cut=3.5)
        yukawa.pair_coeff.set(atom1, 'R', epsilon=0, kappa=1.0, r_cut=0)
        yukawa.pair_coeff.set(atom1, 'Z', epsilon=0, kappa=1.0, r_cut=0)
        yukawa.pair_coeff.set('R', atom1, epsilon=0, kappa=1.0, r_cut=0)
        yukawa.pair_coeff.set('Z', atom1, epsilon=0, kappa=1.0, r_cut=0)
    yukawa.pair_coeff.set('R', 'R', epsilon=0, kappa=1.0, r_cut=0)
    yukawa.pair_coeff.set('R', 'Z', epsilon=0, kappa=1.0, r_cut=0)
    yukawa.pair_coeff.set('Z', 'Z', epsilon=0, kappa=1.0, r_cut=0)

    hoomd.md.integrate.mode_standard(dt=production_dt)

    # Set up integrator
    langevin = hoomd.md.integrate.langevin(group=moving_group, kT=temp, seed=399991)
    for i,name in enumerate(aa_type):
        langevin.set_gamma(name, gamma=aa_mass[i]/1000.0)
    # langevin.set_gamma('R', gamma=1)
    # langevin.set_gamma('Z', gamma=1)
    # langevin.set_gamma('R', gamma=prot_mass_arr/1000.0)
    # langevin.set_gamma('Z', gamma=prot_mass_arr/1000.0)

    hoomd.dump.gsd('rigid_FUS_start.gsd', period=1, group=all_group, truncate=True)
    hoomd.run(1000)