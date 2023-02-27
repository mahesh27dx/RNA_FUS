# -*- coding: utf-8 -*-
#!/usr/bin/python3

# this script starts an MD simulation of multiple chains of a chosen sequence of the protein hnRNPA1
# short-range pairwise potential: HPS model (accounts for both protein-protein and protein-splvent interactions)

import sys,os
import numpy as np
import gsd, gsd.hoomd, gsd.pygsd
import hoomd, hoomd.md
# import azplugins
try:
    from hoomd import azplugins
except ImportError:
    import azplugins
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import hoomd_util as hu

# UNITS: distance -> nm   (!!!positions and sigma in files are in agstrom!!!)
#        mass -> amu
#        energy -> kJ/mol

# Input files
stat_file = 'input_files/stats_module.dat'
# filein_FUS = 'input_files/calpha_6G99.pdb'
filein_FUS = 'input_files/FUS_LC.dat'

dt = 0.001    # time step (in HOOMD units; equals dt = 1e-14)
T = 300

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
        aa_sigma.append(aa_param_dict[i][2]/10.) # divide by 10 to consvert angs-> nm
        aa_lambda.append(aa_param_dict[i][3])

    ff_para = stat_file
    aalist = {}
    with open(ff_para, 'r') as fid:
        for i in fid:
            if i[0] != '#':
                tmp = i.rsplit()
                aalist[tmp[0]] = np.loadtxt(tmp[1:], dtype=float)
    aakeys = list(aalist.keys())
    # This translates each amino acid into a number, which will be used in HOOMD
    #  For example, GLY is with an ID of 10
    aamass = []
    aacharge = []
    aaradius = []
    aahps = []
    print('aakeys')
    print(aakeys)
    print('aalist[i][1]')
    print(aalist[aakeys[1][1]])
    for i in aakeys:
        aamass.append(aalist[i][0])
        aacharge.append(aalist[i][1])
        aaradius.append(aalist[i][2])
        aahps.append(aalist[i][3])
    ## Now we can translate the entire sequence into a number code according to the order in 'aakeys'
    chain_id = []
    chain_mass = []
    chain_charge = []
    
    prot_id, prot_mass, prot_charge, prot_sigma, prot_position = hu.aa_stats_sequence(filein_FUS, aa_param_dict)
    prot_position_array = np.array(prot_position)/10. # converts the position list into an array
    print(prot_position_array.shape)
    # exit()
    prot_length = len(prot_id)
    prot_total_mass = np.sum(prot_mass)
    prot_mass_arr = np.array(prot_mass)

    # Build a HOOMD snapshot of the system
    snap = gsd.hoomd.Snapshot()

    snap.particles.N = len(prot_position_array)
    print(f"The length of chain::{len(prot_position_array)} ")

    snap.particles.types = aa_type  # types of particles

    # Type ID's of the polymer chain
    id_chain = prot_id[:]
    type_id = np.arange(0, len(prot_position_array))
    prot_id = np.array(prot_id)
    snap.particles.typeid = type_id
    print(f"The shape of the type id:::{type_id.shape}")

    mass_chain = prot_mass[:]
    mass = np.empty(0, dtype=float)
    mass = np.concatenate((mass_chain), axis=None)
    snap.particles.mass = mass
    print(f"The shape of mass:::{mass.shape}")

    charge_chain = prot_charge[:]
    charge = np.empty(0, dtype=float)
    charge = np.concatenate((charge_chain), axis=None)
    snap.particles.charge = charge
    print(f"The shape of charge:::{charge.shape}")

    # pos_chain = prot_position_array[:].reshape(-1, 3)
    print(f"The shape of position is::{prot_position_array.shape} ")
    snap.particles.position = prot_position_array

    # The number of bonds between the residues
    # bonds = np.empty((0, 2), dtype=int)
    # for i in range(len(prot_position_array) - 1):
        # bonds = np.append(bonds, [[i, i+1]], axis=0)
    bonds = prot_length -1
    snap.bonds.N = bonds
    snap.bonds.types = ['AA_bond']
    snap.bonds.typeid = [0]*(bonds)
    # print(f"The number of bonds:::{len(bonds)}")
    bond_pairs = np.zeros((bonds, 2), dtype=int)
    for i in range(bonds):
        bond_pairs[i, :] = np.array([i, i+1])
    snap.bonds.group = bond_pairs

    bond_length = 0.381
    # box_length = bond_length * prot_length + 10
    box_length = 30
    print(box_length)
    print(prot_length)
    # exit()
    ## Dimensions of the box
    snap.configuration.dimensions = 3
    snap.configuration.box = [box_length, box_length, box_length, 0, 0, 0]
    snap.configuration.step = 0

    ## Moment of inertia of the system
    # MOI_chain = np.zeros((len(prot_charge[:]), 3), dtype=float)
    # moment_inertia = np.concatenate((MOI_chain), axis=None).reshape(-1,3)
    # snap.particles.moment_inertia = moment_inertia
    # print(f"The shape of moment of inertia:::{moment_inertia.shape}")

    ## Orientation of the system
    # orien_chain = np.zeros((len(prot_charge[:]), 4), dtype=float)
    # orientation = np.concatenate((orien_chain), axis=0).reshape(-1, 4)
    # snap.particles.orientation = orientation
    # print(f"The shape of orientation:::{orientation.shape}")

    ## Composite body associated with each particle. Value [-1] indicates no body
    # body = np.empty(0, dtype=int)
    # body = np.append(body, [-1] * len(charge_chain))
    # snap.particles.body = body
    # print(f"The shape of body:::{body.shape}")
    # Write the snapshot to the file
    with gsd.hoomd.open(name='starting_config.gsd', mode='wb') as f:
        f.append(snap)
        f.close()

    # hoomd.context.initialize()
    #
    # system = hoomd.init.read_gsd('starting_config.gsd')
    # # # system.replicate(nx=nx, ny=ny, nz=nz)
    # # snapshot = system.take_snapshot()
    #
    # # Harmonic potential between bonded particles
    # harmonic = hoomd.md.bond.harmonic()
    # harmonic.bond_coeff.set('AA_bond', k = 8368, r0 = bond_length)
    # #
    # # # Neighbourslist and exclusions
    # nl = hoomd.md.nlist.cell()
    # nl.reset_exclusions(exclusions=['1-2', 'body'])
    # #
    # # # Non-bonded Pairwise interactions for two rigid bodies
    # nb = azplugins.pair.ashbaugh(r_cut = 0, nlist = nl)
    # for i in aa_type:
    #     for j in aa_type:
    #         nb.pair_coeff.set(i, j, lam = (aa_param_dict[i][3] + aa_param_dict[j][3])/2.,
    #                         epsilon=0.8368, sigma=(aa_param_dict[i][2] + aa_param_dict[j][2])/10./2., r_cut=2.0)
    #
    # # # # Electrostatics
    # yukawa = hoomd.md.pair.yukawa(r_cut=0.0, nlist=nl)
    # for i, atom1 in enumerate(aa_type):
    #     for j, atom2 in enumerate(aa_type):
    #         # yukawa.pair_coeff.set(atom1, atom2, epsilon=1.73136, kappa=1.0, r_cut=3.5)
    #         yukawa.pair_coeff.set(atom1, atom2, epsilon=aa_param_dict[atom1][1]*aa_param_dict[atom2][1]*1.73136, kappa=1.0, r_cut=3.5)
    #
    # # # # Grouping of the particles
    # all_group = hoomd.group.all()
    #
    # hoomd.md.integrate.mode_standard(dt=dt)
    # #
    # # # Set up integrator
    # langevin = hoomd.md.integrate.langevin(group=all_group, kT=300, seed=63535)
    #
    # #hoomd.update.box_resize(L=hoomd.variant.linear_interp([(0,system.box.Lx),
    # #                        (resize_steps-500,box_length)]),scale_particles=True)
    # # for i,name in enumerate(aa_type):
    # #     langevin.set_gamma(name, gamma=aa_mass[i]/1000.0)
    # # langevin.set_gamma('R', gamma=1)
    #
    # hoomd.dump.gsd('rigid_FUS_start.gsd', period=100, group=all_group, truncate=True)
    # #hoomd.analyze.log(filename="potential_ene.log", quantities=['potential_energy'], period=100, overwrite=False)
    # hoomd.run(10000)
