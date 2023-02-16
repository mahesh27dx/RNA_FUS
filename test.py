# -*- coding: utf-8 -*-

import argparse
import gsd, gsd.hoomd
import hoomd, hoomd.md
import numpy as np
import sys,os
import warnings

import hoomd_util as hu

warnings.simplefilter(action='ignore', category=FutureWarning)

parser = argparse.ArgumentParser(description='Create two rigid bodies of the protein')
parser.add_argument('-f', '--foo', help='Enter the size of the simulation box(L)', required=True)
parser.add_argument('-b', '--bar', hep='Enter the path of the stats module file', required=True)
#box_length=50

# # rigid_body1_length = 284
# stat_file = 'input_stats/stats_module.dat'
# # filein_ck1d = 'input_stats/CA_ck1delta.pdb'
# filein_FUS = 'input_stats/calpha_FUS.pdb'

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
        
    
#     N = len(pos) 

#     hoomd.context.initialize("")

#     s = hoomd.data.make_snapshot(N=N, box=hoomd.dat.boxdim(Lx=Lx, Ly=Ly, Lz=Lz),
#                                 particle_types=['R', 'X'], bond_types=['polymer'],
#                                 angle_types=['angles180', 'angles90'])

#     s.particles.position[:] = pos[:]
    
#     system = hoomd.init.read_snapshot(snapshot)
    # rigid = hoomd.md.constrain.rigid()
    # rigid.set_param('R', types=[aa_type[FUS_id[i]] for i in range(285,371)],
                    # positions=FUS_rel_pos[285:371])
    
    # rigid.set_param('Z', types=[aa_type[FUS_id[i]] for i in range(422, 453)],
                    # positions=FUS_rel_pos[422:453])
                    
  #  print(rigid.create_bodies(False))
    # rigid.create_bodies()
    
