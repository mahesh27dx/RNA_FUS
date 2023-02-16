# -*- coding: utf-8 -*-

import sys,os
import numpy as np
import gsd, gsd.hoomd
import hoomd, hoomd.md
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import hoomd_util as hu

# UNITS: distance -> nm   (!!!positions and sigma in files are in agstrom!!!)
#        mass -> amu
#        energy -> kJ/mol
# ### MACROs
production_dt=0.01 # Time step for production run in picoseconds
box_length=50

stat_file = 'input_stats/stats_module.dat'
filein_FUS = 'input_stats/calpha_FUS.pdb'

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
    FUS_rel_pos_rigid1 = FUS_pos_arr[284:371] - FUS_cof_pos_rigid1
    FUS_rel_pos_rigid2 = FUS_pos_arr[421:453] - FUS_cof_pos_rigid2
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
        bond_pairs_poly1[i,:] = np.array([i+2, i+3])
    # print(bond_pairs_poly1)
    
    nbonds_FUS_poly2 = 49
    bond_pairs_poly2 = np.zeros((nbonds_FUS_poly2, 2),dtype=int)
    for i in range(0, nbonds_FUS_poly2):
        bond_pairs_poly2[i,:] = np.array([nbonds_FUS_poly1 + i + 3, nbonds_FUS_poly1 + i + 4])
    print("################")
    # print(bond_pairs_poly2)
    
    nbonds_FUS_poly3 = 70
    bond_pairs_poly3 = np.zeros((nbonds_FUS_poly3, 2),dtype=int)
    for i in range(0, nbonds_FUS_poly3):
        bond_pairs_poly3[i,:] = np.array([nbonds_FUS_poly1 + nbonds_FUS_poly2 + i + 4, nbonds_FUS_poly1 + nbonds_FUS_poly2 + i + 5])
    #print(bond_pairs_poly3)
    
    # Now we can build HOOMD data structure for one single frame
    #number_of_rigid_particles = 121

    s=gsd.hoomd.Snapshot()
    
    ## Types of the rigid bodies and the amino acids
    s.particles.types = aa_type  +['R'] + ['Z']    # 21 types
    
    # exit()
    ## Postions of the residues in the system
    # FUS_pos_arr_poly = np.concatenate((FUS_pos_arr[:284].reshape(-1,3), FUS_pos_arr[371:421].reshape(-1,3), FUS_pos_arr[455:].reshape(-1,3)), axis=0) 
    # FUS_pos_arr_poly_re = FUS_pos_arr_poly
    # position = np.empty((0, 3), dtype=float)
    position = np.concatenate((FUS_pos_arr[:284].reshape(-1,3), FUS_cof_pos_rigid1, 
                               FUS_pos_arr[371:421].reshape(-1,3), FUS_cof_pos_rigid2, 
                               FUS_pos_arr[455:].reshape(-1,3)), axis=0)
    # position = np.append([FUS_cof_pos_rigid], FUS_pos_arr_poly_re)
    s.particles.position = position
    
    ## Total number of particles in the system
    s.particles.N = len(position) 
    print(s.particles.N)
    ## Total mass of the residues in the system
    mass=np.concatenate((FUS_mass[:284],np.sum(FUS_mass[284:371]), 
                         FUS_mass[371:421], np.sum(FUS_mass[421:454]),
                         FUS_mass[455:]),axis=None)
    #print("mass",mass.shape)
    s.particles.mass = mass

    ## Total charge on the residues in the system
    charge = np.empty(0, dtype=float)
    charge = np.concatenate((FUS_charge[:284], np.sum(FUS_charge[284:371]), 
                             FUS_charge[371:421], np.sum(FUS_charge[421:454]), 
                             FUS_charge[455:]),axis=None)
    s.particles.charge = charge

    ## Moment of inertia of the two rigid bodies
    # moment_inertia_rigid1 = np.zeros((len(FUS_pos_arr[:284]), 3), dtype=float)
    moment_inertia = np.concatenate((np.zeros((len(FUS_charge[:284]), 3), dtype=float), 
                                     I_diag_rigid1, np.zeros((len(FUS_charge[371:421]), 3), dtype=float),
                                     I_diag_rigid2, np.zeros((len(FUS_charge[455:]), 3), dtype=float)), axis=0)
    s.particles.moment_inertia = moment_inertia

    ## Orientation of each particle (size:: (N, 4))
    orientation = np.concatenate((np.zeros((len(FUS_charge[:284]), 4), dtype=float), 
                                  np.array([[0,1,1,1]]), np.zeros((len(FUS_charge[371:421]), 4), dtype=float), 
                                  np.array([[0,1,1,1]]), np.zeros((len(FUS_charge[455:]), 4), dtype=float)), axis=0)
    s.particles.orientation = orientation

    ## Stores the composite body associated with each particle. Value [-1] indicates no body.
    body = np.empty(0, dtype=int)
    body = np.append(body, [-1] * len(FUS_charge[:284]) + [284] + [-1] * len(FUS_charge[371:421]) 
                                                        + [334] + [-1] * len(FUS_charge[455:]))
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
    s.bonds.types = ['1typebon']
    typeid = np.arange(0, len(position))
    FUS_id = np.array(FUS_id)
    print(typeid.shape, FUS_id.shape)
    print(np.sort(FUS_id))
    print(len(typeid))
    typeid[:284]=FUS_id[:284]
    typeid[284]=20
    typeid[285:285+50] = FUS_id[371:421]
    typeid[335]=21
    typeid[336:]=FUS_id[455:]
    print(typeid,typeid.shape)
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

    with gsd.hoomd.open(name='FUS_start.gsd', mode='wb') as f:
        f.append(s)
        f.close()

    # exit()    
    
    # Defining rigid body
    hoomd.context.initialize()
    system = hoomd.init.read_gsd('FUS_start.gsd')
    all_p = hoomd.group.all()

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

    #system.bonds.add('AA_bond', 284, 441)
   # system.bonds.add('AA_bond', 285, 526)

    hoomd.dump.gsd('rigid_FUS_start.gsd', period=1, group=all_p, truncate=True)
    hoomd.run_upto(1, limit_hours=24)