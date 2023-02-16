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
    FUS_cof_pos_rigid = np.append(FUS_cof_pos_rigid1, FUS_cof_pos_rigid2,axis=0)

    # Rigid body I moment of inertia
    I_rigid1 = hu.protein_moment_inertia(FUS_rel_pos_rigid1, FUS_mass[284:371])
    I_diag_rigid1, E_vec_rigid1 = np.linalg.eig(I_rigid1)
    I_diag_rigid1=I_diag_rigid1.reshape(-1,3)

    # Rigid body II moment of inertia
    I_rigid2 = hu.protein_moment_inertia(FUS_rel_pos_rigid2, FUS_mass[421:453])
    I_diag_rigid2, E_vec_rigid2 = np.linalg.eig(I_rigid2)
    I_diag_rigid2 =I_diag_rigid2.reshape(-1,3)
    
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
    print(bond_pairs_poly3)
    

    # Now we can build HOOMD data structure for one single frame
    number_of_rigid_particles = 121

    s=gsd.hoomd.Snapshot()
    
    ## Types of the rigid bodies and the amino acids
    s.particles.types = ['R'] + ['Z'] + aa_type    # 21 types
    
    ## TYPE ID's of the polymer chains
    FUS_id_poly1 = FUS_id[:284]
    FUS_id_poly2 = FUS_id[371:421]
    FUS_id_poly3 = FUS_id[455:]
    type_id = np.empty(0, dtype=float)
    type_id = np.append([len(aa_type)] + [len(aa_type)],  FUS_id_poly1 + FUS_id_poly2 + FUS_id_poly3)
    s.particles.typeid = type_id
   
    ## Postions of the residues in the system
    # FUS_pos_arr_poly = np.concatenate((FUS_pos_arr[:284].reshape(-1,3), FUS_pos_arr[371:421].reshape(-1,3), FUS_pos_arr[455:].reshape(-1,3)), axis=0) 
    # FUS_pos_arr_poly_re = FUS_pos_arr_poly
    # position = np.empty((0, 3), dtype=float)
    position = np.concatenate((FUS_pos_arr[:284].reshape(-1,3), FUS_cof_pos_rigid1, 
                               FUS_pos_arr[371:421].reshape(-1,3), FUS_cof_pos_rigid2, 
                               FUS_pos_arr[455:].reshape(-1,3)), axis=0)
    # position = np.append([FUS_cof_pos_rigid], FUS_pos_arr_poly_re)
    s.particles.position = np.append( FUS_cof_pos_rigid, FUS_pos_arr_poly_re,axis=0)
    
    ## Total number of particles in the system
    s.particles.N = len(s.particles.position) 

    ## Total mass of the residues in the system
    mass = np.empty(0, dtype=float)
    mass = np.append(mass, [np.sum(FUS_mass[284:371]) , np.sum(FUS_mass[421:454])] + FUS_mass[:284] + FUS_mass[371:421] + FUS_mass[455:])
    s.particles.mass = mass

    charge = np.empty(0, dtype=float)
    charge = np.append(charge, [np.sum(FUS_charge[284:371]), np.sum(FUS_charge[421:454])] + FUS_charge[:284] + FUS_charge[371:421] + FUS_charge[455:])
    s.particles.charge = charge

    print(f"FUS_pos:::{type(FUS_pos_arr[:284])}")
    print(f"I_diag_rigid1:::{I_diag_rigid1.ndim}")
    moment_inertia_rigid1 = np.zeros((len(FUS_pos_arr[:284]), 3), dtype=float)
    print(f"The shape moment_inertia rigid 1::::{moment_inertia_rigid1.ndim}")
    # print(f"The shape moment_inertia rigid 2::::{np.shape(moment_inertia_rigid2)}")
    
    # MOI_rigid_1 = np.append(moment_inertia_rigid1, I_diag_rigid1, axis=1)
    moment_inertia = np.concatenate((np.zeros((len(FUS_charge[:284]),3),dtype=float),I_diag_rigid1,np.zeros((len(FUS_charge[371:421]),3),dtype=float),I_diag_rigid2,np.zeros((len(FUS_charge[455:]),3),dtype=float)),axis=0)
    # moment_inertia_rigid2 = np.zeros((len(FUS_pos_arr[371:421]),3),dtype=float)
    # moment_inertia_rigid3 = np.zeros((len(FUS_pos_arr[455:]),3),dtype=float)
    # print(f"The shape moment_inertia rigid 1::::{moment_inertia_rigid1}")
    # print(f"The shape moment_inertia rigid 2::::{np.shape(moment_inertia_rigid2)}")
    # print(f"The type moment_inertia::::{type(moment_inertia)}")
    
    s.particles.moment_inertia = moment_inertia
    # print(f"The len moment_inertia::::{s.particles.moment_inertia}")
    # s.particles.moment_inertia = [I_diag_rigid1[0], I_diag_rigid1[1], I_diag_rigid1[2], I_diag_rigid2[0], I_diag_rigid2[1], I_diag_rigid2[2]] + [0,0,0]*(FUS_length - number_of_rigid_particles + 1)

    orientation = np.empty((0, 4), dtype=float)
    orientation = np.append(orientation, [(1, 0, 0, 0)]* (FUS_length - number_of_rigid_particles + 2), axis=0)
    s.particles.orientation = orientation

    body = np.empty(0, dtype=float)
    body = np.append(body, [0] + [-1]*(FUS_length - number_of_rigid_particles + 1))
    s.particles.body = body

    s.bonds.N = nbonds_FUS_poly1 + nbonds_FUS_poly2 + nbonds_FUS_poly3
    s.bonds.types = ['AA_bond']
    s.bonds.typeid = [0]*(nbonds_FUS_poly1 + nbonds_FUS_poly2 + nbonds_FUS_poly3)
    s.bonds.group = np.concatenate((bond_pairs_poly1, bond_pairs_poly2, bond_pairs_poly3), axis=0)
    # print(f"s.bonds.N:::{s.bonds.N}")
    # print(f"FUS_rel_pos_rigid1:::{len(FUS_rel_pos_rigid2)}")
    # print(f"s.bonds.typeid:::{len(s.bonds.typeid)}")
    # print(f"FUS_pos_arr[421:453]:::{len(FUS_pos_arr[421:454])}") 
    #  print(f"s.particles.MOI:::{len(s.particles.moment_inertia)}")
    # print(f"s.particles.typeid:::{len(s.particles.typeid)}")
    # print(f"s.bonds.group:::{s.bonds.group.shape}")
    # print(f"s.particles.mass:::{s.particles.mass.shape}")
    
    # exit()
    s.configuration.dimensions = 3
    s.configuration.box = [box_length,box_length,box_length,0,0,0]
    s.configuration.step = 0

    # print(f"FUS_rel_pos_rigid1:::{type(FUS_rel_pos_rigid1)}")
    # print(f"FUS_rel_pos_rigid1:::{len(FUS_rel_pos_rigid2)}")
    # print(f"FUS_pos_arr[284:371]:::{len(FUS_pos_arr[284:371])}")
    # print(f"FUS_pos_arr[421:453]:::{len(FUS_pos_arr[421:454])}") 
    #  print(f"s.particles.MOI:::{len(s.particles.moment_inertia)}")
    # print(f"s.particles.typeid:::{len(s.particles.typeid)}")
    # print(f"FUS ID:::{type(FUS_id[284:371])}")
    # print(f"s.particles.mass:::{s.particles.mass.shape}")
    
    rel_charge1 = [aa_type[FUS_id[i]] for i in range(284, 371)]
    rel_charge2 = [aa_type[FUS_id[i]] for i in range(421, 453)]
    print(len(rel_charge1),len(FUS_rel_pos_rigid1[:]))
    print(len(rel_charge2),len(FUS_rel_pos_rigid2[:]))
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

    # system.bonds.add('AA_bond', 284, 440)
    # system.bonds.add('AA_bond', 285, 526)

    hoomd.dump.gsd('rigid_FUS_start.gsd', period=1, group=all_p, truncate=True)
    hoomd.run_upto(1, limit_hours=24)