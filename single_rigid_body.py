# -*- coding: utf-8 -*-

import sys,os
import numpy as np
import gsd, gsd.hoomd
import hoomd, hoomd.md

import hoomd_util as hu

# UNITS: distance -> nm   (!!!positions and sigma in files are in agstrom!!!)
#        mass -> amu
#        energy -> kJ/mol
# ### MACROs
production_dt=0.01 # Time step for production run in picoseconds
box_length=50

# rigid_body1_length = 284
stat_file = 'input_files/stats_module.dat'
# filein_ck1d = 'input_stats/CA_ck1delta.pdb'
filein_FUS = 'input_files/calpha_FUS.pdb'

# ############################
# ### ONE RIGID BODY #################
# ###########################
# if __name__=='__main__':
#     # Input parameters for all the amino acids 
#     aa_param_dict = hu.aa_stats_from_file(stat_file)
#     aa_type = list(aa_param_dict.keys())    
#     aa_mass = []
#     aa_charge = []
#     aa_sigma = []
#     aa_lambda = []
#     for k in aa_type:
#         aa_mass.append(aa_param_dict[k][0])
#         aa_charge.append(aa_param_dict[k][1])
#         aa_sigma.append(aa_param_dict[k][2]/10.)
#         aa_lambda.append(aa_param_dict[k][3])
    
#     # Now we can translate the entire sequence into a sequence of numbers and 
#     # assign to each a.a. of the sequence its stats
#     FUS_id, FUS_mass, FUS_charge, FUS_sigma, FUS_pos = hu.aa_stats_sequence(filein_FUS, aa_param_dict)
#     FUS_pos_arr = np.array(FUS_pos)/10.
#     # FUS_pos_arr = FUS_pos_arr + 7
#     print(f"The FUS_pos_arr:::{FUS_pos_arr.shape}")
#     FUS_length = len(FUS_id)       
#     print(f"the length of FUS::{FUS_length}")
#     FUS_tot_mass = np.sum(FUS_mass)   
#     FUS_cof_pos = (np.sum(FUS_pos_arr[:,0] * FUS_mass)/FUS_tot_mass, 
#                    np.sum(FUS_pos_arr[:,1] * FUS_mass)/FUS_tot_mass, 
#                    np.sum(FUS_pos_arr[:,2] * FUS_mass)/FUS_tot_mass)
    
#     print(f"FUS_cof_pos::{FUS_cof_pos}")
#     FUS_rel_pos = FUS_pos_arr - FUS_cof_pos
#     print(f"The FUS_rel_pos:::{FUS_rel_pos.shape}")
#     # Rigid body I moment of inertia
#     I_rigid1 = hu.protein_moment_inertia(FUS_rel_pos, FUS_mass)
#     I_diag_rigid1, E_vec_rigid1 = np.linalg.eig(I_rigid1)
#     FUS_diag_pos_rigid1 = np.dot( E_vec_rigid1.T, FUS_rel_pos.T).T

#     # Initialize bond
#     nbonds_FUS_poly1 = 300
#     bond_pairs_poly1 = np.zeros((nbonds_FUS_poly1, 2),dtype=int)
#     for i in range(0, nbonds_FUS_poly1):
#         bond_pairs_poly1[i,:] = np.array([i+2, i+3])
#     # 
#     nbonds_FUS_poly2 = 126
#     bond_pairs_poly2 = np.zeros((nbonds_FUS_poly2, 2),dtype=int)
#     for i in range(0, nbonds_FUS_poly2):
#         bond_pairs_poly2[i,:] = np.array([nbonds_FUS_poly1 + i + 4, nbonds_FUS_poly1 + i + 5])
#     # print(f"Bond pairs of polymer chain2::: {bond_pairs_poly2}")

#     #print(f"Bond pairs of polymer chain3::: {bond_pairs_poly1+bond_pairs_poly2+bond_pairs_poly3}")
#     # print(nbonds_FUS_poly1 + nbonds_FUS_poly2 + nbonds_FUS_poly3 + number_of_rigid_particles)
    
    
#     # Now we can build HOOMD data structure for one single frame
#     number_of_rigid_particles = 100
    
#     s=gsd.hoomd.Snapshot()
    
#     # FUS_length = FUS_length -1
#     s.particles.N = FUS_length - number_of_rigid_particles # (FUS_length-86)
#     # s.particles.N = FUS_length  # (FUS_length-86)
#     print(f"s.particles.N::{s.particles.N}")

#     s.particles.types = aa_type + ['R']   # 21 types

#     FUS_id_poly1 = FUS_id[:300]
#     FUS_id_poly2 = FUS_id[401:]
#     # s.particles.typeid = [len(aa_type)] + [len(aa_type)] +  FUS_id_poly1 + FUS_id_poly2 + FUS_id_poly3   
#     # print(f"s.particles.typeid::{len(s.particles.typeid)}")

#     FUS_pos_arr_poly = np.concatenate((FUS_pos_arr[:300], FUS_pos_arr[401:]), axis=None) 
#     print(f"FUS_cof_pos:::{len(FUS_cof_pos)}")
#     # print(f"FUS_pos_arr_poly_re:::{FUS_pos_arr_poly_re.shape}")
#     s.particles.position = np.append( [FUS_cof_pos], FUS_pos_arr_poly)
#     # s.particles.position =  FUS_pos_arr_poly_re
#     print(f"s.particles.position:::{s.particles.position.shape}")


#     s.particles.mass = [FUS_tot_mass] + FUS_mass[:300] + FUS_id[401:]
#     s.particles.charge = [0] + FUS_charge
#     s.particles.moment_inertia = [I_diag_rigid1[0], I_diag_rigid1[1], I_diag_rigid1[2]] + [0,0,0]*(FUS_length - number_of_rigid_particles) 
#     s.particles.orientation = [(1, 0, 0, 0)] * (FUS_length - number_of_rigid_particles + 1)
#     s.particles.body = [0] + [-1]*(FUS_length - number_of_rigid_particles)
    
#     s.bonds.N = nbonds_FUS_poly1 + nbonds_FUS_poly2
#     s.bonds.types = ['AA_bond']
#     s.bonds.typeid = [0]*(nbonds_FUS_poly1 + nbonds_FUS_poly2)
#     s.bonds.group = np.concatenate((bond_pairs_poly1, bond_pairs_poly2), axis=None)
#     # print(s.bonds.group.shape)
    
#     s.configuration.dimensions = 3
#     s.configuration.box = [box_length,box_length,box_length,0,0,0] 
#     s.configuration.step = 0
    
#     with gsd.hoomd.open(name='FUS_start.gsd', mode='wb') as f:
#         f.append(s)
#         f.close()
        
#     # Defining rigid body
#     hoomd.context.initialize()
#     system = hoomd.init.read_gsd('FUS_start.gsd')
#     all_p = hoomd.group.all()
    
#     snapshot = system.take_snapshot()
#     #print(snapshot.particles.body)

#     rigid = hoomd.md.constrain.rigid()
#     rigid.set_param('R', types=[aa_type[FUS_id[i]] for i in range(285,371)],
#                     positions=FUS_rel_pos[285:371])
    
#   #  print(rigid.create_bodies(False))
#     rigid.create_bodies()

#     system.bonds.add('AA_bond', 284, 440)
    
#     hoomd.dump.gsd('rigid_FUS_start.gsd', period=1, group=all_p, truncate=True)
#     hoomd.run_upto(1, limit_hours=24)
