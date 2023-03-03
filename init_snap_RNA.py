# -*- coding: utf-8 -*-
#!/usr/bin/python3

"""
this scri

"""

import datetime
import os
import sys
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


## Creating subfolders for results and setting date and time
DATETIME = datetime.datetime.now()
current_datetime = ( str(DATETIME.day) + '.' + str(DATETIME.month) + '.' + str(DATETIME.year))

## create directory in folder 'results' to save the results from single simulation.
# output = os.makedirs("output_files/" + current_datetime + "/")
# with open('output/testing.txt', 'w') as f:
#     f.close()

## Input files
stat_file = 'input_files/stats_module.dat'
rna_stat_file = 'input_files/rna_stats.dat'
filein_FUS = 'input_files/calpha_FUS.pdb'
polyAlength = int(50)
##  SIMULATION PARAMETERS
dt = 0.001

if __name__=='__main__':
    ## Input parameters for RNA
    # rna_param_dict = hu.rna_stats_from_file(rna_stat_file)
    # rna_type = list(rna_param_dict.keys())
    # rna_mass = []
    # rna_charge = []
    # rna_sigma = []
    # rna_lambda = []
    # for i in rna_type:
    #     rna_mass.append(rna_param_dict[i][0])
    #     rna_charge.append(rna_param_dict[i][1])
    #     rna_sigma.append(rna_param_dict[i][2]/10.) # divide by 10 to consvert angs-> nm
    #     rna_lambda.append(rna_param_dict[i][3])

    ## Input parameters for all amino acids
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

    ## Now we can translate the entire sequence into a sequence of numbers and
    # assign to each Amino acid of the sequence its stats
    prot_id, prot_mass, prot_charge, prot_sigma, prot_position = hu.aa_stats_sequence(filein_FUS, aa_param_dict)
    # Add RNA chains

    for i in range(polyAlength):
        # prot_id.append(len(rna_type)-1)
        prot_id.append(len(aa_type)-1)
        prot_mass.append(329.2)
        prot_charge.append(-1)
    # print(f"The rna_type:::{len(rna_type)-1}")
    # print(f"The aa_type:::{aa_type}")
    # exit()
    prot_position_array = np.array(prot_position) / 10.
    # prot_position_array = prot_position_array + 8
    prot_length = len(prot_id)
    prot_total_mass = np.sum(prot_mass)
    prot_mass_arr = np.array(prot_mass)

    # Relative positions of the constituent particles of the two rigid bodies
    ### First rigid body relative position and moment of inertia
    rigid_1_coord_1 = np.sum(prot_position_array[284:371, 0] * prot_mass[284:371]) / np.sum(prot_mass[284:371])
    rigid_1_coord_2 = np.sum(prot_position_array[284:371, 1] * prot_mass[284:371]) / np.sum(prot_mass[284:371])
    rigid_1_coord_3 = np.sum(prot_position_array[284:371, 2] * prot_mass[284:371]) / np.sum(prot_mass[284:371])
    relative_pos_rigid_1 = prot_position_array[284:371] - np.array([[rigid_1_coord_1, rigid_1_coord_2, rigid_1_coord_3]])
    relative_pos_rigid_1 = relative_pos_rigid_1

    I_general_rigid_1 = hu.protein_moment_inertia(relative_pos_rigid_1, prot_mass[284:371]) # Moment of inertia
    I_diag_rigid_1, E_vec_rigid_1 = np.linalg.eig(I_general_rigid_1)
    I_diag_rigid_1 = I_diag_rigid_1.reshape(-1, 3)

    ## Second rigid body relative position and moment of inertia
    rigid_2_coord_1 = np.sum(prot_position_array[421:453, 0] * prot_mass[421:453]) / np.sum(prot_mass[421:453])
    rigid_2_coord_2 = np.sum(prot_position_array[421:453, 1] * prot_mass[421:453]) / np.sum(prot_mass[421:453])
    rigid_2_coord_3 = np.sum(prot_position_array[421:453, 2] * prot_mass[421:453]) / np.sum(prot_mass[421:453])
    relative_pos_rigid_2 = prot_position_array[421:453] - np.array([[rigid_2_coord_1, rigid_2_coord_2, rigid_2_coord_3]])
    relative_pos_rigid_2 = relative_pos_rigid_2

    I_general_rigid_2 = hu.protein_moment_inertia(relative_pos_rigid_2, prot_mass[421:453])
    I_diag_rigid_2, E_vec_rigid_2 = np.linalg.eig(I_general_rigid_2)
    I_diag_rigid_2 = I_diag_rigid_2.reshape(-1, 3)

    ## Position vectors of two rigid bodies
    relative_pos_rigid = np.append([[rigid_1_coord_1], [rigid_1_coord_2],
                                    [rigid_1_coord_3]], [[rigid_2_coord_1],
                                    [rigid_2_coord_2], [rigid_2_coord_3]]).reshape(-1, 3)

    ## Build a HOOMD snapshot of the system
    snap = gsd.hoomd.Snapshot()

    # Build initial position of RNA as a linear chain
    rna_bond_length = 0.5
    rna_pos = []
    for i in range(polyAlength):
        rna_pos.append((10, 0, (i - int(polyAlength/2))*rna_bond_length))
        # print((i-int(polyAlength/2))*rna_bond_length)
    rna_pos = np.array(rna_pos)
    print(f"The shape of the RNA_pos::{rna_pos.shape}")
    pos_chain_1 = prot_position_array[:284].reshape(-1, 3)
    pos_chain_2 = prot_position_array[371:421].reshape(-1, 3)
    pos_chain_3 = prot_position_array[453:].reshape(-1, 3)

    rigid_1_coord = np.array([[rigid_1_coord_1, rigid_1_coord_2, rigid_1_coord_3]])
    rigid_2_coord = np.array([[rigid_2_coord_1, rigid_2_coord_2, rigid_2_coord_3]])
    position = np.concatenate((pos_chain_1, rigid_1_coord, pos_chain_2, rigid_2_coord,
                                pos_chain_3), axis=0)
    snap.particles.N = len(position) + polyAlength
    # print(f"The total number of beads:::{len(position) + polyAlength}")


    new_pos = np.append(position, rna_pos).reshape(-1,3)
    snap.particles.position = new_pos

    # types of particles, ['R'] and ['Z'] are two rigid bodies
    # snap.particles.types = aa_type + ['R'] + ['Z']  + rna_type
    snap.particles.types = aa_type + ['R'] + ['Z']
    # print(f"The length of the prot_mass::{len(prot_mass)}")

    mass_chain_1 = prot_mass[:284]
    mass_chain_2 = prot_mass[371:421]
    mass_chain_3 = prot_mass[453:526]
    rna_mass = prot_mass[526:]
    mass_rigid_1 = np.sum(prot_mass[284:371])
    mass_rigid_2 = np.sum(prot_mass[421:453])
    mass = np.concatenate((mass_chain_1, mass_rigid_1, mass_chain_2,
                            mass_rigid_2, mass_chain_3, rna_mass), axis=None)
    # print(f"The shape of the mass:::{mass.shape}")
    snap.particles.mass = mass

    #charge = np.empty(0, dtype=float)
    charge_chain_1 = prot_charge[:284]
    charge_chain_2 = prot_charge[371:421]
    charge_chain_3 = prot_charge[453:526]
    rna_charge = prot_charge[526:]
    charge_rigid_1 = np.sum(prot_charge[284:371])
    charge_rigid_2 = np.sum(prot_charge[421:453])
    charge = np.concatenate((charge_chain_1, charge_rigid_1, charge_chain_2,
                            charge_rigid_2, charge_chain_3, rna_charge), axis=None)
    # print(f"The shape of the charge:::{charge.shape}")
    snap.particles.charge = charge

    ## Moment of inertia of the system
    MOI_chain_1 = np.zeros((len(prot_charge[:284]), 3), dtype=float)
    MOI_chain_2 = np.zeros((len(prot_charge[371:421]), 3), dtype=float)
    MOI_chain_3 = np.zeros((len(prot_charge[453:526]), 3), dtype=float)
    MOI_rna = np.zeros((len(prot_charge[526:]), 3), dtype=float)
    moment_inertia = np.concatenate((MOI_chain_1, I_diag_rigid_1, MOI_chain_2,
                                    I_diag_rigid_2, MOI_chain_3, MOI_rna), axis=0)
    snap.particles.moment_inertia = moment_inertia

    ## Orientation of the system
    orien_chain_1 = np.zeros((len(prot_charge[:284]), 4), dtype=float)
    orien_chain_2 = np.zeros((len(prot_charge[371:421]), 4), dtype=float)
    orien_chain_3 = np.zeros((len(prot_charge[453:526]), 4), dtype=float)
    orien_rna = np.zeros((len(prot_charge[526:]), 4), dtype=float)
    # orien_rigid_2 = np.array([[0, 1, 1, 1]])
    orien_rigid_1 = np.array([[1, 0, 0, 0]])
    orien_rigid_2 = np.array([[1, 0, 0, 0]])
    orientation = np.concatenate((orien_chain_1, orien_rigid_1, orien_chain_2,
                                    orien_rigid_2, orien_chain_3, orien_rna), axis=0)
    snap.particles.orientation = orientation

    ## Composite body associated with each particle. Value [-1] indicates no body
    body = np.empty(0, dtype=int)
    body = np.append(body, [-1] * len(charge_chain_1) + [284] + [-1] * len(charge_chain_2)
                        + [334] + [-1] * len(charge_chain_3) + [-1] * len(rna_charge))
    # print(f"The shape of the body:::{body.shape}")
    snap.particles.body = body

    ## The number of bonds between the residues
    bond_length = 0.381

    bonds = np.empty((0, 2), dtype=int)
    for i in range(len(position) - 1):
        bonds = np.append(bonds, [[i, i+1]], axis=0)


    rna_bonds = np.empty((0, 2), dtype=int)
    for i in range(polyAlength - 1):
        rna_bonds = np.append(rna_bonds, [[i, i+1]], axis=0)

    new_bonds = np.append(bonds, rna_bonds).reshape(-1,2)
    print(f"The shape of the bonds:::{bonds.shape}")
    print(f"The shape of the rna_bonds:::{rna_bonds.shape}")
    print(f"The shape of the new_bonds:::{new_bonds.shape}")

    snap.bonds.N = len(new_bonds)

    bond_pairs = np.zeros((len(new_bonds), 2), dtype=int)
    for i in range(0, len(position) - 1):
        print('%s-%s-A' % (i, i+1))
        bond_pairs[i, :] = np.array([i, i+1])

    for cnt, i in enumerate(range(len(position),  len(new_bonds)+1)):
        print('%s-%s-B' % (i, i+1))
        bond_pairs[cnt+len(position) - 1, :] = np.array([i, i+1])
    # exit()

    snap.bonds.group = bond_pairs

    # exit()
    snap.bonds.types = ['AA_bond', 'rigid_1', 'rigid_2', 'NT_bonds']

    ## Type ID's of the polymer chain
    poly_chain_one = prot_id[:284]
    # id_chain_1 = prot_id[:284]
    id_chain_1 = poly_chain_one
    # id_chain_2 = prot_id[371:421]
    # id_chain_3 = prot_id[453:]
    # print(f"The len(position):::{len(position)}")
    # print(f"The len(new_pos):::{len(new_pos)}")
    type_id = np.arange(0, len(new_pos))
    # print(f"The type_id before assigning prot_id:::{len(type_id)}")
    prot_id = np.array(prot_id)

    type_id[:284] = prot_id[:284]
    type_id[284] = 20
    type_id[285:285+50] = prot_id[371:421]
    type_id[335] = 21
    # print(len(prot_id))
    # print(len(type_id))
    type_id[336:409] = prot_id[453:526]
    # type_id[336:409] = prot_id[453:526]
    # type_id[405:413] = prot_id[526:530]
    # print(len(prot_id))
    # print(len(type_id))

    # print(len(prot_id[526:]))
    # print(len(type_id[336:]))
    # exit()
    # rna_id = np.arange(0, len(rna_pos))
    # print(len(rna_id))
    # rna_id = np.array(prot_id[526:])
    type_id[409:] = prot_id[526:]

    snap.particles.typeid = type_id

    #box_length = bond_length * prot_length + 10
    Lx = 100
    Ly = 100
    Lz = 150

    ## Dimensions of the box
    snap.configuration.dimensions = 3
    snap.configuration.box = [Lx, Ly, Lz, 0, 0, 0]
    snap.configuration.step = 0

    ## Write the snapshot to the file
    with gsd.hoomd.open(name='output_files/starting_config_RNA.gsd' , mode='wb') as fout:
        fout.append(snap)
        fout.close()

    hoomd.context.initialize("")

    ## Read the "starting_config.gsd"
    system = hoomd.init.read_gsd('output_files/starting_config_RNA.gsd')
    # no of chains nx=3, ny=3, nz=3=27
    # system.replicate(nx=3, ny=3, nz=3)
    snapshot = system.take_snapshot()

    ## Types for the rigid bodies
    type_rigid_1 = [aa_type[prot_id[i]] for i in range(284, 371)]
    type_rigid_2 = [aa_type[prot_id[i]] for i in range(421, 453)]

    rigid = hoomd.md.constrain.rigid()

    ## First rigid body
    rigid.set_param('R',
                    types = type_rigid_1,
                    positions = relative_pos_rigid_1)

    ## Second rigid body
    rigid.set_param('Z',
                    types = type_rigid_2,
                    positions = relative_pos_rigid_2)

    rigid.create_bodies()

    ## Grouping of the particles
    all_group = hoomd.group.all()

    hoomd.dump.gsd('output_files/init_RNA_snap.gsd', period=None, group=all_group, overwrite=True)

    hoomd.run(0)