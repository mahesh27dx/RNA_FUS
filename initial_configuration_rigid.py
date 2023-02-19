# -*- coding: utf-8 -*-

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
box_length = 1000

stat_file = 'input_files/stats_module.dat'
filein_FUS = 'input_files/calpha_FUS.pdb'

if __name__=='__main__':
    # Input parameters for all the amino acids
    aa_param_dict = hu.aa_stats_from_file(stat_file)
    aa_type = list(aa_param_dict.keys())
    aa_mass = []
    aa_charge = []
    aa_sigma = []
    aa_lambda = []
    for k in aa_type:
        aa_mass.append(aa_param_dict[k][0])
        aa_charge.append(aa_param_dict[k][1])
        aa_sigma.append(aa_param_dict[k][2]/10.)
        aa_lambda.append(aa_param_dict[k][3])

    # Now we can translate the entire sequence into a sequence of numbers and
    # assign to each a.a. of the sequence its stats
    FUS_id, FUS_mass, FUS_charge, FUS_sigma, FUS_pos = hu.aa_stats_sequence(filein_FUS, aa_param_dict)
    FUS_pos_arr = np.array(FUS_pos)/10.
    FUS_length = len(FUS_id)
    FUS_tot_mass = np.sum(FUS_mass)

    FUS_cof_pos_rigid1 = np.array([[np.sum(FUS_pos_arr[284:371,0] * FUS_mass[284:371])/np.sum(FUS_mass[284:371]),
                   np.sum(FUS_pos_arr[284:371,1] * FUS_mass[284:371])/np.sum(FUS_mass[284:371]),
                   np.sum(FUS_pos_arr[284:371,2] * FUS_mass[284:371])/np.sum(FUS_mass[284:371])]])

    FUS_cof_pos_rigid2 = np.array([[np.sum(FUS_pos_arr[421:453,0] * FUS_mass[421:453])/np.sum(FUS_mass[421:453]),
                   np.sum(FUS_pos_arr[421:453,1] * FUS_mass[421:453])/np.sum(FUS_mass[421:453]),
                   np.sum(FUS_pos_arr[421:453,2] * FUS_mass[421:453])/np.sum(FUS_mass[421:453])]])
    FUS_rel_pos_rigid1 = (FUS_pos_arr[284:371] - FUS_cof_pos_rigid1)
    FUS_rel_pos_rigid2 = (FUS_pos_arr[421:453] - FUS_cof_pos_rigid2)
    FUS_cof_pos_rigid = np.append(FUS_cof_pos_rigid1, FUS_cof_pos_rigid2, axis=0)

    # Rigid body I moment of inertia
    I_rigid1 = hu.protein_moment_inertia(FUS_rel_pos_rigid1, FUS_mass[284:371])
    I_diag_rigid1, E_vec_rigid1 = np.linalg.eig(I_rigid1)
    I_diag_rigid1 = I_diag_rigid1.reshape(-1,3)

    # Rigid body II moment of inertia
    I_rigid2 = hu.protein_moment_inertia(FUS_rel_pos_rigid2, FUS_mass[421:453])
    I_diag_rigid2, E_vec_rigid2 = np.linalg.eig(I_rigid2)
    I_diag_rigid2 = I_diag_rigid2.reshape(-1,3)

    # Initialize bond
    nbonds_FUS_poly1 = 283
    bond_pairs_poly1 = np.zeros((nbonds_FUS_poly1, 2),dtype=int)
    for i in range(0, nbonds_FUS_poly1):
        bond_pairs_poly1[i,:] = np.array([i, i+1])
    # print(bond_pairs_poly1)

    nbonds_FUS_poly2 = 49
    bond_pairs_poly2 = np.zeros((nbonds_FUS_poly2, 2),dtype=int)
    for i in range(0, nbonds_FUS_poly2):
        bond_pairs_poly2[i,:] = np.array([nbonds_FUS_poly1 + i + 2, nbonds_FUS_poly1 + i + 3])
    print("################")
    # print(bond_pairs_poly2)

    nbonds_FUS_poly3 = 70
    bond_pairs_poly3 = np.zeros((nbonds_FUS_poly3, 2),dtype=int)
    for i in range(0, nbonds_FUS_poly3):
        bond_pairs_poly3[i,:] = np.array([nbonds_FUS_poly1 + nbonds_FUS_poly2 + i + 4, nbonds_FUS_poly1 + nbonds_FUS_poly2 + i + 5])
    # print(bond_pairs_poly3)

    # Now we can build HOOMD data structure for one single frame
    #number_of_rigid_particles = 121

    s=gsd.hoomd.Snapshot()

    ## Types of the rigid bodies and the amino acids
    s.particles.types = aa_type  +['R'] + ['Z']    # 21 types

    position = np.concatenate((FUS_pos_arr[:284].reshape(-1,3), FUS_cof_pos_rigid1,
                               FUS_pos_arr[371:421].reshape(-1,3), FUS_cof_pos_rigid2,
                               FUS_pos_arr[453:].reshape(-1,3)), axis=0)
    # position = np.append([FUS_cof_pos_rigid], FUS_pos_arr_poly_re)
    s.particles.position = position

    ## Total number of particles in the system
    s.particles.N = len(position)
    print(s.particles.N)

    ## Total mass of the residues in the system
    mass=np.concatenate((FUS_mass[:284],np.sum(FUS_mass[284:371]),
                         FUS_mass[371:421], np.sum(FUS_mass[421:453]),
                         FUS_mass[453:]),axis=None)
    #print("mass",mass.shape)
    s.particles.mass = mass

    ## Total charge on the residues in the system
    charge = np.empty(0, dtype=float)
    charge = np.concatenate((FUS_charge[:284], np.sum(FUS_charge[284:371]),
                             FUS_charge[371:421], np.sum(FUS_charge[421:453]),
                             FUS_charge[453:]),axis=None)
    s.particles.charge = charge

    ## Moment of inertia of the two rigid bodies
    moment_inertia = np.concatenate((np.zeros((len(FUS_charge[:284]), 3), dtype=float),
                                     I_diag_rigid1, np.zeros((len(FUS_charge[371:421]), 3), dtype=float),
                                     I_diag_rigid2, np.zeros((len(FUS_charge[453:]), 3), dtype=float)), axis=0)
    s.particles.moment_inertia = moment_inertia

    ## Orientation of each particle (size:: (N, 4))
    orientation = np.concatenate((np.zeros((len(FUS_charge[:284]), 4), dtype=float),
                                  np.array([[0,1,1,1]]), np.zeros((len(FUS_charge[371:421]), 4), dtype=float),
                                  np.array([[0,1,1,1]]), np.zeros((len(FUS_charge[453:]), 4), dtype=float)), axis=0)
    s.particles.orientation = orientation

    ## Stores the composite body associated with each particle. Value [-1] indicates no body.
    body = np.empty(0, dtype=int)
    body = np.append(body, [-1] * len(FUS_charge[:284]) + [284] + [-1] * len(FUS_charge[371:421])
                                                        + [334] + [-1] * len(FUS_charge[453:]))
    s.particles.body = body

    ## TYPE ID's of the polymer chains
    # FUS_id_poly1 = FUS_id[:284]
    # FUS_id_poly2 = FUS_id[371:421]
    # FUS_id_poly3 = FUS_id[455:]
    # print(len(FUS_id))
    # type_id = np.empty(0, dtype=int)
    # type_id = np.append(type_id, FUS_id_poly1 + FUS_id_poly2+FUS_id_poly3 + [20] + [21])
    # print(type_id.shape)
    # # type_id = np.append((type_id,rigid_1_id, rigid_2_id, FUS_id_poly1, FUS_id_poly2, FUS_id_poly3))
    # s.particles.typeid = type_id
    # #print(type(s.particles.typeid))

    ## The number of bonds between the residues.
    bonds = np.empty((0,2),dtype=int)
    for i in range(len(position) - 1):
        bonds = np.append(bonds, [[i, i+1]], axis=0)
    s.bonds.group = bonds
    s.bonds.N = len(bonds)
    s.bonds.types = ['AA_bond','rigid1', 'rigid2']
    typeid = np.arange(0, len(position))
    FUS_id = np.array(FUS_id)
    # print(typeid.shape, FUS_id.shape)
    # print(np.sort(FUS_id))
    # print(len(typeid))
    typeid[:284]=FUS_id[:284]
    typeid[284]=20
    typeid[285:285+50] = FUS_id[371:421]
    typeid[335]=21
    typeid[336:]=FUS_id[453:]
    # print(typeid.shape,position.shape,moment_inertia.shape,orientation.shape)
    # print(typeid,typeid.shape)
    s.particles.typeid=typeid
    # s.bonds.N = nbonds_FUS_poly1 + nbonds_FUS_poly2 + nbonds_FUS_poly3
    # s.bonds.types = ['AA_bond']
    # s.bonds.typeid = [0] * (nbonds_FUS_poly1 + nbonds_FUS_poly2 + nbonds_FUS_poly3)
    # s.bonds.group = np.concatenate((bond_pairs_poly1, bond_pairs_poly2, bond_pairs_poly3), axis=0)

    s.configuration.dimensions = 3
    s.configuration.box = [box_length,box_length,box_length,0,0,0]
    s.configuration.step = 0

    rel_charge1 = [aa_type[FUS_id[i]] for i in range(284, 371)]
    rel_charge2 = [aa_type[FUS_id[i]] for i in range(421, 453)]
    # print(len(rel_charge1),len(FUS_rel_pos_rigid1[:]))
    # print(len(rel_charge2),len(FUS_rel_pos_rigid2[:]))
    # exit()
    # print(len(FUS_rel_pos_rigid1),len(FUS_rel_pos_rigid2))
    with gsd.hoomd.open(name='FUS_start.gsd', mode='wb') as f:
        f.append(s)
        f.close()

    # exit()

    # Defining rigid body
    hoomd.context.initialize("--mode=cpu")
    system = hoomd.init.read_gsd('FUS_start.gsd')
    #all_p = hoomd.group.all()

    snapshot = system.take_snapshot()

    rigid = hoomd.md.constrain.rigid()
    rigid.set_param('R',
                    types= rel_charge1,
                    positions=FUS_rel_pos_rigid1[:])
    rigid.set_param('Z',
                    types=rel_charge2,
                    positions=FUS_rel_pos_rigid2[:])
  #  print(rigid.create_bodies(False))
    rigid.create_bodies()
    for i in range(len(FUS_rel_pos_rigid1)):
        # print(284, 528+i-len(FUS_rel_pos_rigid1)-len(FUS_rel_pos_rigid2))
        system.bonds.add('rigid1', 284, 528+i-len(FUS_rel_pos_rigid1)-len(FUS_rel_pos_rigid2))
    for i in range(len(FUS_rel_pos_rigid2)-1):
        system.bonds.add('rigid2', 335, 528+i-len(FUS_rel_pos_rigid2))
    #system.bonds.add('AA_bond', 285, 286)

    all_group = hoomd.group.all()
    # Group particles that belong to rigid bodies
    center_group = hoomd.group.rigid_center()
    non_rigid_group = hoomd.group.nonrigid()
    moving_group = hoomd.group.union('moving_group', center_group, non_rigid_group)

    # Harmonic potential between bonded particles
    bond_length = 0.38
    harmonic = hoomd.md.bond.harmonic()
    harmonic.bond_coeff.set(['AA_bond', 'rigid1', 'rigid2'], k=8368, r0=bond_length)

    #Neighborlist and exclusions
    nl = hoomd.md.nlist.cell()
    nl.reset_exclusions(exclusions=['1-2', 'body'])

    # Pairwise interactions
    nb = azplugins.pair.ashbaugh(r_cut=0, nlist=nl)
    for i in aa_type:
        for j in aa_type:
            nb.pair_coeff.set(i, j, lam=(aa_param_dict[i][3] + aa_param_dict[j][3])/2., epsilon=0.8368, sigma=(aa_param_dict[i][2] + aa_param_dict[j][2])/10./2., r_cut=3.0)
        nb.pair_coeff.set(i, 'R', lam=0., epsilon=0, sigma=0, r_cut=0)
        nb.pair_coeff.set(i, 'Z', lam=0., epsilon=0, sigma=0, r_cut=0)
    nb.pair_coeff.set('R', 'R', lam=0., epsilon=0, sigma=0, r_cut=0)
    nb.pair_coeff.set('R', 'Z', lam=0., epsilon=0, sigma=0, r_cut=0)
    nb.pair_coeff.set('R', 'Z', lam=0., epsilon=0, sigma=0, r_cut=0)
    nb.pair_coeff.set('Z', 'Z', lam=0., epsilon=0, sigma=0, r_cut=0)
    # Electrostatics
    yukawa = hoomd.md.pair.yukawa(r_cut=0.0, nlist=nl)
    # yukawa.pair_coeff.set('R','Z', epsilon=1.73136, kappa=1.0, r_cut=3.5)
    for i, atom1 in enumerate(aa_type):
        atom1 = aa_type[i]
        yukawa.pair_coeff.set(atom1, 'R', epsilon=1.73136, kappa=1.0, r_cut=3.5)
        for j, atom2 in enumerate(aa_type):
            atom2 = aa_type[j]
            yukawa.pair_coeff.set(atom1, atom2, epsilon=aa_param_dict[atom1][1]*aa_param_dict[atom2][1]*1.73136, kappa=1.0, r_cut=3.5)
        yukawa.pair_coeff.set(atom1, 'R', epsilon=0, kappa=1.0, r_cut=0)
        yukawa.pair_coeff.set(atom1, 'Z', epsilon=0, kappa=1.0, r_cut=0)
    yukawa.pair_coeff.set('R', 'R', epsilon=0, kappa=1.0, r_cut=0)
    yukawa.pair_coeff.set('R', 'Z', epsilon=0, kappa=1.0, r_cut=0)
    yukawa.pair_coeff.set('Z', 'Z', epsilon=0, kappa=1.0, r_cut=0)

    # Set up integrator
    resize_dt=0.01 # Time step in picoseconds for box resizing
    resize_T=300 # Temperature for resizing run in Kelvin
    production_dt = 10000
    temp = resize_T * 8.3144598/1000.
    hoomd.md.integrate.mode_standard(dt=production_dt)
    # exit()
    langevin = hoomd.md.integrate.langevin(group=moving_group, kT=temp, seed=399991)
    for i,name in enumerate(aa_type):
        langevin.set_gamma(name, gamma=aa_mass[i]/1000.0)
    langevin.set_gamma('R', gamma=0.0001)
    langevin.set_gamma('Z', gamma=0.0001)
    # langevin.set_gamma('R', gamma=FUS_mass_arr[:]/1000.0)
    # langevin.set_gamma('Z', gamma=FUS_mass_arr[:]/1000.0)

    hoomd.dump.gsd('rigid_FUS_start.gsd', period=1, group=all_group, truncate=True)
    hoomd.run(1)
