"""
Scripts for initializing multiple RNA chains around the protein molecule
"""

import datetime
import itertools
import os
import sys
import numpy as np
import gsd, gsd.hoomd
import hoomd, hoomd.md
try:
    import azplugins
except:
    from hoomd import azplugins
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import hoomd_util as hu

# Input files
stat_file = 'input_files/stats_module.dat'
stat_rna = 'input_files/rna_stats.dat'
filein_FUS = 'input_files/calpha_FUS.pdb'

rna_length = int(250)
n_rna_molecules = int(2)
# Slab dimensions
Lx, Ly, Lz = 15, 15, 280
box_length = Lx * Ly * Lz

if __name__ == '__main__':
    # Input parameters for all RNA bases
    rna_param_dict = hu.rna_stats_from_file(stat_rna)
    rna_type = list(rna_param_dict.keys())
    print(f"The rna_types are:::{rna_type}")
    rna_mass = []
    rna_charge = []
    rna_sigma = []
    rna_lambda = []
    for i in rna_type:
        rna_mass.append(rna_param_dict[i][0])
        rna_charge.append(rna_param_dict[i][1])
        rna_sigma.append(rna_param_dict[i][2]/10.)
        rna_lambda.append(rna_param_dict[i][3])
    # print(f"The mass of the RNA :::{rna_mass}")
    # print(f"The charge of the RNA :::{rna_charge}")
    # print(f"The sigma of the RNA :::{rna_sigma}")
    # print(f"The lambda of the RNA :::{rna_lambda}")

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
        aa_sigma.append(aa_param_dict[i][2]/10.) #divide by 10 to consvert angs-> nm
        aa_lambda.append(aa_param_dict[i][3])

    # Now we can assign to amino acids to the sequence read from the .pdb file
    prot_id, prot_mass, prot_charge, prot_sigma, prot_position =            hu.aa_stats_sequence(filein_FUS, aa_param_dict)
    prot_position_array = np.array(prot_position) / 10.
    prot_length = len(prot_id)
    prot_total_mass = np.sum(prot_mass)
    prot_mass_arr = np.array(prot_mass)
    # print(prot_id)

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

    # Add RNA chain
    rna_id = []
    for j in range(n_rna_molecules):
        for i in range(rna_length):
            rna_id.append(len(rna_type)-1)
    # print(f"The rna_id is :::{rna_id}")

    # Build initial position of RNA as a linear polymer chain
    rna_bond_length = 0.5
    rna_positions = []
    for j in range(n_rna_molecules):
        for i in range(rna_length):
            rna_positions.append((5, 0, (i - int(rna_length/2)) * rna_bond_length))
        # print((i - int(rna_length/2)) * rna_bond_length)
    rna_positions = np.array(rna_positions)
    print(f"The shape of the RNA pos:::{rna_positions.shape}")

    pos_chain_1 = prot_position_array[:284].reshape(-1, 3)
    pos_chain_2 = prot_position_array[371:421].reshape(-1, 3)
    pos_chain_3 = prot_position_array[453:].reshape(-1, 3)
    rigid_1_coord = np.array([[rigid_1_coord_1, rigid_1_coord_2, rigid_1_coord_3]])
    rigid_2_coord = np.array([[rigid_2_coord_1, rigid_2_coord_2, rigid_2_coord_3]])
    protein_positions = np.concatenate((pos_chain_1, rigid_1_coord, pos_chain_2,
                                        rigid_2_coord,pos_chain_3), axis=0)
    system_positions = np.append(protein_positions, rna_positions).reshape(-1,3) # append order can be changed to (rna_pos, position)
    print(f"The shape of the Protein positions:::{protein_positions.shape}")
    print(f"The shape of the system positions:::{system_positions.shape}")

    # The number of bonds between the amino acids
    protein_bonds = np.empty((0, 2), dtype=int)
    for i in range(len(protein_positions) - 1):
        protein_bonds = np.append(protein_bonds, [[i, i+1]], axis=0)
    print(f"The shape of the protein bonds::{protein_bonds.shape}")

    rna_bonds = np.empty((0, 2), dtype=int)
    for j in range(n_rna_molecules):
        for i in range(len(protein_bonds) + rna_length):
            if i >= 409:
            # print("The bonds are not correct!!!")
        # else:
                rna_bonds = np.append(rna_bonds, [[i, i+1]], axis=0)
    print(f"The shape of the RNA bonds::{rna_bonds.shape}")

    system_bonds = np.append(protein_bonds, rna_bonds).reshape(-1, 3)
    print(f"The shape of the system bonds:::{system_bonds.shape}")
    nbonds = len(system_bonds)

    bond_pairs = np.zeros((len(system_bonds), 2), dtype=int)
    for j in range(n_rna_molecules):
        for i in range(0, len(protein_positions) - 1):
        #print('%s-%s-A' % (i, i+1))
            bond_pairs[i, :] = np.array([i, i+1])
    #print(f"The shape of the protein bond pairs::{bond_pairs}")
        for cnt, i in enumerate(range(len(protein_positions), len(system_bonds) + 1)):
            #print('%s-%s-B' % (i, i+1))
            bond_pairs[i-1, :] = np.array([i, i+1])
    #print(f"The shape of the protein bond pairs::{bond_pairs}")

    # Type Id's
    type_id = np.arange(0, len(system_positions))
    print(f"The shape of the type id::{type_id.shape}")
    #protein chain
    protein_id = np.array(prot_id)
    print(f"The shape of the protein id::{protein_id.shape}")
    type_id[:284] = protein_id[:284]
    type_id[284] = 20
    type_id[285:285+50] = protein_id[371:421]
    type_id[335] = 21
    type_id[336:409] = protein_id[453:526]
    # RNA molecule
    rna_id = np.array(rna_id)
    print(f"The shape of the RNA id::{rna_id.shape}")
    type_id[409:659] = rna_id[:250]
    # print(type_id[409:659] == rna_id[:250])
    type_id[659:] = rna_id[250:]
    # print(type_id[659:] == rna_id[250:])

    # Mass of the system
    mass_chain_1 = prot_mass[:284]
    mass_chain_2 = prot_mass[371:421]
    mass_chain_3 = prot_mass[453:526]
    rna_mass = rna_mass[526:]
    mass_rigid_1 = np.sum(prot_mass[284:371])
    mass_rigid_2 = np.sum(prot_mass[421:453])
    mass = np.concatenate((mass_chain_1, mass_rigid_1, mass_chain_2,
                            mass_rigid_2, mass_chain_3, rna_mass), axis=None)
    print(f"The shape of the mass:::{mass.shape}")

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
                        + [334] + [-1] * len(charge_chain_3) + [-2] * len(rna_charge))
