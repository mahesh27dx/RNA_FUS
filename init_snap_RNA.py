# -*- coding: utf-8 -*-
#!/usr/bin/python3

"""
this scri
"""
import datetime
import itertools
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
stat_rna = 'input_files/rna_stats.dat'
filein_FUS = 'input_files/calpha_FUS.pdb'
rna_length = int(250)
#box_length = bond_length * prot_length + 10
"""
These length dimensions are working
Lx = 75
Ly = 75
Lz = 150
"""
Lx = 15
Ly = 15
Lz = 280
box_length = Lx * Ly * Lz
# dt = 0.001
# T = 300

if __name__=='__main__':
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
    print(rna_type)
    # exit()
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
    prot_position_array = np.array(prot_position) / 10.
    prot_length = len(prot_id)
    prot_total_mass = np.sum(prot_mass)
    prot_mass_arr = np.array(prot_mass)

    # Add RNA chains
    rna_id = []
    for i in range(rna_length):
        rna_id.append(len(rna_type))
        prot_mass.append(329.2)
        prot_charge.append(-1)
    # print(f"The rna_type:::{len(rna_type)-1}")
    # print(f"The aa_type:::{aa_type}")

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

     # Build initial position of RNA as a linear chain
    rna_bond_length = 0.5
    rna_pos = []
    for i in range(rna_length):
        rna_pos.append((5, 0, (i - int(rna_length/2))*rna_bond_length))
        # print((i-int(polyAlength/2))*rna_bond_length)
    rna_pos = np.array(rna_pos)
    """
    rna_pos_2 = []
    for i in range(rna_length):
        rna_pos_2.append((10, 0, (i - int(rna_length/2))*rna_bond_length))
        # print((i-int(polyAlength/2))*rna_bond_length)
    rna_pos_2 = np.array(rna_pos_2)
    # print(f"The shape of the RNA_pos::{rna_pos.shape}")
    rna_pos = np.append(rna_pos_1, rna_pos_2).reshape(-1, 3)
    """
    pos_chain_1 = prot_position_array[:284].reshape(-1, 3)
    pos_chain_2 = prot_position_array[371:421].reshape(-1, 3)
    pos_chain_3 = prot_position_array[453:].reshape(-1, 3)
    rigid_1_coord = np.array([[rigid_1_coord_1, rigid_1_coord_2, rigid_1_coord_3]])
    rigid_2_coord = np.array([[rigid_2_coord_1, rigid_2_coord_2, rigid_2_coord_3]])
    position = np.concatenate((pos_chain_1, rigid_1_coord, pos_chain_2, rigid_2_coord,
                                pos_chain_3), axis=0)
    new_pos = np.append(position, rna_pos).reshape(-1,3) # append order can be changed to (rna_pos, position)
    # print(f"The shape of the new_position:::{new_pos[1]}")

    ## The number of bonds between the residues
    bond_length = 0.381
    poly_bonds = np.empty((0, 2), dtype=int)
    for i in range(len(position) - 1):
        poly_bonds = np.append(poly_bonds, [[i, i+1]], axis=0)
    # print(f"The length of the poly_bonds::{poly_bonds}")
    print("#######################")

    rna_bonds = np.empty((0, 2), dtype=int)
    for i in range(len(poly_bonds) + rna_length):
        if i >= 409:
            rna_bonds = np.append(rna_bonds, [[i, i+1]], axis=0)
    print(f"The length of rna_bonds:::{len(rna_bonds)}")

    new_bonds = np.append(poly_bonds, rna_bonds).reshape(-1,2)

    nbonds = len(new_bonds)

    bond_pairs = np.zeros((len(new_bonds), 2), dtype=int)
    for i in range(0, len(position)-1):
        # print('%s-%s-A' % (i, i+1))
        bond_pairs[i, :] = np.array([i, i+1])

    for i in range(len(position), len(new_bonds)+1):
        # print('%s-%s-B' % (i, i+1))
        bond_pairs[i-1,:] = np.array([i,i+1])


    # Type ID's of the polymer chain
    type_id = np.arange(0, len(new_pos))
    prot_id = np.array(prot_id)
    # print(f"The shape of the prot_id::{prot_id.shape}")
    # print(f"The shape of the type_id::{type_id.shape}")
    print(len(rna_id[:]))
    print(prot_id[:])
    # exit()
    type_id[:284] = prot_id[:284]
    type_id[284] = 20
    type_id[285:285+50] = prot_id[371:421]
    type_id[335] = 21
    type_id[336:409] = prot_id[453:526]
    type_id[409:] = rna_id[:]

    # Mass of the system
    mass_chain_1 = prot_mass[:284]
    mass_chain_2 = prot_mass[371:421]
    mass_chain_3 = prot_mass[453:526]
    rna_mass = prot_mass[526:]
    mass_rigid_1 = np.sum(prot_mass[284:371])
    mass_rigid_2 = np.sum(prot_mass[421:453])
    mass = np.concatenate((mass_chain_1, mass_rigid_1, mass_chain_2,
                            mass_rigid_2, mass_chain_3, rna_mass), axis=None)
    # print(f"The shape of the mass:::{mass.shape}")

    # Charge of the system
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

    ## Moment of inertia of the system
    MOI_chain_1 = np.zeros((len(prot_charge[:284]), 3), dtype=float)
    MOI_chain_2 = np.zeros((len(prot_charge[371:421]), 3), dtype=float)
    MOI_chain_3 = np.zeros((len(prot_charge[453:526]), 3), dtype=float)
    MOI_rna = np.zeros((len(prot_charge[526:]), 3), dtype=float)
    moment_inertia = np.concatenate((MOI_chain_1, I_diag_rigid_1, MOI_chain_2,
                                    I_diag_rigid_2, MOI_chain_3, MOI_rna), axis=0)

    ## Orientation of the system
    orien_chain_1 = np.zeros((len(prot_charge[:284]), 4), dtype=float)
    orien_chain_2 = np.zeros((len(prot_charge[371:421]), 4), dtype=float)
    orien_chain_3 = np.zeros((len(prot_charge[453:526]), 4), dtype=float)
    orien_rna = np.zeros((len(prot_charge[526:]), 4), dtype=float)
    orien_rigid_1 = np.array([[1, 0, 0, 0]])
    orien_rigid_2 = np.array([[1, 0, 0, 0]])
    orientation = np.concatenate((orien_chain_1, orien_rigid_1, orien_chain_2,
                                    orien_rigid_2, orien_chain_3, orien_rna), axis=0)

    ## Composite body associated with each particle. Value [-1] indicates no body
    body = np.empty(0, dtype=int)
    body = np.append(body, [-1] * len(charge_chain_1) + [284] + [-1] * len(charge_chain_2)
                        + [334] + [-1] * len(charge_chain_3) + [-1] * len(rna_charge))
    # print(f"The shape of the body:::{body.shape}")
    # exit()

    ## Build a HOOMD snapshot of the system
    snap = gsd.hoomd.Snapshot()
    snap.particles.N = len(position) + rna_length
    snap.particles.types = aa_type + ['R'] + ['Z']
    snap.particles.position = new_pos
    snap.particles.mass = mass
    snap.particles.charge = charge
    snap.particles.moment_inertia = moment_inertia
    snap.particles.orientation = orientation
    snap.particles.body = body
    snap.bonds.N = len(bond_pairs)
    snap.bonds.group = bond_pairs
    snap.bonds.types = ['AA_bond', 'rigid_1', 'rigid_2', 'NT_bond']
    snap.particles.typeid = type_id
    snap.configuration.dimensions = 3
    snap.configuration.box = [Lx, Ly, Lz, 0, 0, 0]
    snap.configuration.step = 0

    # Write the snapshot to the file
    with gsd.hoomd.open(name='output_files/starting_config_RNA_test.gsd' , mode='wb') as fout:
        fout.append(snap)
        fout.close()

    hoomd.context.initialize("")

    ## Read the "starting_config.gsd"
    system = hoomd.init.read_gsd('output_files/starting_config_RNA.gsd')
    snapshot = system.take_snapshot(bonds=True, pairs=True)

    # Types for the rigid bodies
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
    print(len(position)+rna_length+ 86+ 32)
    system.bonds.add('NT_bond', 283, len(position)+rna_length)
    # system.bonds.add('AA_bond', 500, 285)
    system.bonds.add('NT_bond', len(position)+rna_length+ 86, 285)
    system.bonds.remove(283)
    system.bonds.remove(284)
    # system.bonds.add('AA_bond', 334, 501)
    system.bonds.add('NT_bond', 334, len(position)+rna_length+ 87)
    system.bonds.remove(334)
    system.bonds.remove(335)
    # system.bonds.add('AA_bond', 336, 532)
    system.bonds.add('NT_bond', 336, len(position)+rna_length+ 86 + 32)
    # system.replicate(nx=3, ny=3, nz=3)
    ## Grouping of the particles
    all_group = hoomd.group.all()

    hoomd.dump.gsd('output_files/init_RNA_snap_test.gsd', period=1, group=all_group, overwrite=True)

    hoomd.run(1)
