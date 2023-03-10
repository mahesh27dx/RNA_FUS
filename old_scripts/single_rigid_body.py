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
filein_FUS = 'input_files/calpha_FUS.pdb'

# SYSTEM SIZE
nx=5 # Replicate system in x direction
ny=5 # Replicate system in y direction
nz=4 # Replicate system in z direction --> Total system size would be nx*ny*nz
# Please keep magnitudes of nx, ny and nz as close to each other
# as possible to avoid large replicated system box sizes in a specific direction

# ### SIMULATION PARAMETERS
saveTrajectories = 1
numberGSDframes = 1000 # the number of snapshots the sim writes into the .gsd file

dt = 0.005    # time step (in HOOMD units; equals dt = 1e-14)
steps = int(5000)
T_in_K = 300

Lx_in_nm = 15                    # box size x-direction in nm (1D = 4,5A --> 15 nm = 33.333 D)
Ly_in_nm = 15                    # box size y-direction in nm (1D = 4,5A --> 15 nm =  33.333 D)
Lz_in_nm = 150                    # box size z-direction in nm (1D = 4,5A --> 75 nm =  166.667 D)

# Cutoff radius of electrostatic interactions in nm (Coulombic term; HOOMD: Yukawa potential)
rcut_YU_in_nm = 4.0
# Cutoff radius of short range interactions in nm (Ashbaugh Hatch functional form)
rcut_AH_in_nm = 2.0

# SLAB DIMENSIONS
boxsize=15.0 # The x and y dimensions of the slab configuration in nanometers
slab_z_length = 280.0 # The z dimension of the slab configuration in nanometers

# box_length = 900
# box_length = Lx*Ly*Lz
# HOOMD units
# Absolute energy scale of the short-ranged interactions between all pairs of
# amino acids (1 epsilon = 0.2 kcal/mol)
epsilon = 1.0
r0 = 0.38 / 0.45    # bond equilibrium position (see harmonic potential, r0=3.8A, 1D= 4.5 A)
# spring constant (see harmonic potential, k=20J/m^2), high=> stiff: bond length constrained
# to be almost exactly sigma.
k = 2025.0
D = 80.0    # dielectric constant of the solvent (see Yukawa potential)
kappa = 1/2.222 # inverse (due to def of yukawa) of Debye screening length (~1nm)

# Convert temperature and box size in HOOMD units:
k_B = 1.38064852e-23    # Boltzmann constant [J/K]
N_A = 6.02214076e23     # Avogadro number [mol^{-1}]
epsilon_in_kcal_mol = 0.2   # energy unit [kcal/mol]: 1 eps = 0.2 kcal/mol
distUnit_in_A = 4.5     # distance unit [A]: 1D = 4.5 angs

# T = round(T_in_K * k_B / (epsilon_in_kcal_mol * 4184 / N_A), 3) # 1 kcal=4184 J
T = 300
Lx = round(10.0 * Lx_in_nm / distUnit_in_A, 3)
Ly = round(10.0 * Ly_in_nm / distUnit_in_A, 3)
Lz = round(10.0 * Lz_in_nm / distUnit_in_A, 3)

rcut_YU = round(10.0 * rcut_YU_in_nm / distUnit_in_A, 3)
rcut_AH = round(10.0 * rcut_AH_in_nm / distUnit_in_A, 3)

resize_dt = 0.01 # Time step in picoseconds for box resizing
resize_T = 300 # Temperature for resizing run in Kelvin
resize_steps=1000

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

    # Now we can translate the entire sequence into a sequence of numbers and
    # assign to each Amino acid of the sequence its stats
    prot_id, prot_mass, prot_charge, prot_sigma, prot_position = hu.aa_stats_sequence(filein_FUS, aa_param_dict)
    prot_position_array = np.array(prot_position)/10. # converts the position list into an array
    prot_length = len(prot_id)
    prot_total_mass = np.sum(prot_mass)
    prot_mass_arr = np.array(prot_mass)
    # poly1 = :284 ## find a clever way to assign the index
    # Relative positions of the constituent particles of the two rigid bodies
    ### First rigid body relative position and moment of inertia
    rigid_coord_1 = np.sum(prot_position_array[284:371, 0] * prot_mass[284:371]) / np.sum(prot_mass[284:371])
    rigid_coord_2 = np.sum(prot_position_array[284:371, 1] * prot_mass[284:371]) / np.sum(prot_mass[284:371])
    rigid_coord_3 = np.sum(prot_position_array[284:371, 2] * prot_mass[284:371]) / np.sum(prot_mass[284:371])
    relative_pos_rigid = prot_position_array[284:371] - np.array([[rigid_coord_1, rigid_coord_2, rigid_coord_3]])
    # print(relative_pos_rigid_1.shape)
    I_general_rigid = hu.protein_moment_inertia(relative_pos_rigid, prot_mass[284:371]) # Moment of inertia
    I_diag_rigid, E_vec_rigid = np.linalg.eig(I_general_rigid)
    I_diag_rigid = I_diag_rigid.reshape(-1, 3)

    # Build a HOOMD snapshot of the system
    snap = gsd.hoomd.Snapshot()
    snap.particles.types = aa_type + ['R']  # types of particles, ['R'] and ['Z'] are two rigid bodies

    pos_chain_1 = prot_position_array[:284].reshape(-1, 3)
    pos_chain_2 = prot_position_array[371:].reshape(-1, 3)
    rigid_coord = np.array([[rigid_coord_1, rigid_coord_2, rigid_coord_3]])
    position = np.concatenate((pos_chain_1, rigid_coord, pos_chain_2), axis=0)
    print(f"The shape of position is::{position.shape} ")
    snap.particles.position = position # positions of the free particles and rigid body constituent particles in the system
    snap.particles.N = len(position)


    mass_chain_1 = prot_mass[:284]
    mass_chain_2 = prot_mass[371:]
    mass_rigid = np.sum(prot_mass[284:371])
    mass = np.concatenate((mass_chain_1, mass_rigid, mass_chain_2,), axis=None)
    snap.particles.mass = mass
    print(f"The shape of mass:::{mass.shape}")

    charge_chain_1 = prot_charge[:284]
    charge_chain_2 = prot_charge[371:]
    charge_rigid = np.sum(prot_charge[284:371])
    charge = np.concatenate((charge_chain_1, charge_rigid, charge_chain_2), axis=None)
    # print(charge.shape)
    snap.particles.charge = charge
    print(f"The shape of charge:::{charge.shape}")

    ## Moment of inertia of the system
    MOI_chain_1 = np.zeros((len(prot_charge[:284]), 3), dtype=float)
    MOI_chain_2 = np.zeros((len(prot_charge[371:]), 3), dtype=float)
    # print(MOI_chain_1)
    moment_inertia = np.concatenate((MOI_chain_1, I_diag_rigid, MOI_chain_2), axis=0)
    # print(I_diag_rigid_2)
    snap.particles.moment_inertia = moment_inertia
    print(f"The shape of moment of inertia:::{moment_inertia.shape}")

    ## Orientation of the system
    orien_chain_1 = np.zeros((len(prot_charge[:284]), 4), dtype=float)
    orien_chain_2 = np.zeros((len(prot_charge[371:]), 4), dtype=float)
    orien_rigid = np.array([[0, 1, 1, 1]])
    orientation = np.concatenate((orien_chain_1, orien_rigid, orien_chain_2), axis=0)
    snap.particles.orientation = orientation
    # print(f"The shape of orientation:::{orientation.shape}")

    ## Composite body associated with each particle. Value [-1] indicates no body
    body = np.empty(0, dtype=int)
    body = np.append(body, [-1] * len(charge_chain_1) + [284] + [-1] * len(charge_chain_2))
    snap.particles.body = body
    # print(f"The shape of body:::{body.shape}")

    # The number of bonds between the residues
    bonds = np.empty((0, 2), dtype=int)
    for i in range(len(position) - 1):
        bonds = np.append(bonds, [[i, i+1]], axis=0)
    snap.bonds.group = bonds
    snap.bonds.N = len(bonds)
    snap.bonds.types = ['AA_bond', 'rigid_1']
    # print(f"The number of bonds:::{len(bonds)}")

    # Type ID's of the polymer chain
    id_chain_1 = prot_id[:284]
    id_chain_2 = prot_id[371:]
    type_id = np.arange(0, len(position))
    prot_id = np.array(prot_id)
    type_id[:284] = prot_id[:284]
    # print(len(type_id[:284]))
    type_id[284] = 20
    type_id[285:] = prot_id[371:]
    # print(len(type_id[284:]))
    # print(len(prot_id[284:]))

    # type_id[285:285+50] = prot_id[371:421]
    # type_id[335] = 21
    # type_id[336:] = prot_id[453:]
    snap.particles.typeid = type_id
    # print(f"The shape of the type id:::{type_id.shape}")

    bond_length = 0.381
    box_length = bond_length * prot_length + 300
    print(box_length)
    print(prot_length)
    # exit()
    ## Dimensions of the box
    snap.configuration.dimensions = 3
    snap.configuration.box = [box_length, box_length, box_length, 0, 0, 0]
    snap.configuration.step = 1

    # Write the snapshot to the file
    with gsd.hoomd.open(name='starting_config.gsd', mode='wb') as f:
        f.append(snap)
        f.close()

    hoomd.context.initialize("--mode=cpu")

    # Read the "starting_config.gsd"
    system = hoomd.init.read_gsd('starting_config.gsd')
    # system.replicate(nx=nx, ny=ny, nz=nz)
    snapshot = system.take_snapshot()

    # Define the rigid bodies
    type_rigid = [aa_type[prot_id[i]] for i in range(284, 371)]

    rigid = hoomd.md.constrain.rigid()

    ### First rigid body
    rigid.set_param('R',
                    types = type_rigid,
                    positions = relative_pos_rigid[:])

    rigid.create_bodies()

    # # Grouping of the particles
    all_group = hoomd.group.all()
    center_group =hoomd.group.rigid_center()
    non_rigid_group = hoomd.group.nonrigid()
    moving_group = hoomd.group.union('moving_group', center_group, non_rigid_group)


    # Harmonic potential between bonded particles
    harmonic = hoomd.md.bond.harmonic()
    harmonic.bond_coeff.set(['AA_bond', 'rigid_1'], k = 8368, r0 = bond_length)

    # Neighbourslist and exclusions
    nl = hoomd.md.nlist.cell()
    nl.reset_exclusions(exclusions=['1-2', 'body'])

    # Non-bonded Pairwise interactions for two rigid bodies
    nb = azplugins.pair.ashbaugh(r_cut = 0, nlist = nl)
    for i in aa_type:
        for j in aa_type:
            nb.pair_coeff.set(i, j, lam = (aa_param_dict[i][3] + aa_param_dict[j][3])/2.,
                            epsilon=0.8368, sigma=(aa_param_dict[i][2] + aa_param_dict[j][2])/10./2., r_cut=2.0)
        nb.pair_coeff.set(i, 'R', lam=0., epsilon=0., sigma=0., r_cut=0.)
        nb.pair_coeff.set('R', i, lam=0., epsilon=0., sigma=0., r_cut=0.)
    nb.pair_coeff.set('R', 'R', lam=0., epsilon=0, sigma=0, r_cut=0)


    # Electrostatics
    yukawa = hoomd.md.pair.yukawa(r_cut=0.0, nlist=nl)
    for i, atom1 in enumerate(aa_type):
        atom1 = aa_type[i]
        for j, atom2 in enumerate(aa_type):
            atom2 = aa_type[j]
            yukawa.pair_coeff.set(atom1, atom2, epsilon=aa_param_dict[atom1][1]*aa_param_dict[atom2][1]*1.73136, kappa=1.0, r_cut=3.5)
        yukawa.pair_coeff.set(atom1, 'R', epsilon=0, kappa=1.0, r_cut=0)
        # yukawa.pair_coeff.set('R', atom1, epsilon=0, kappa=1.0, r_cut=0)
    yukawa.pair_coeff.set('R', 'R', epsilon=0, kappa=1.0, r_cut=0)


    hoomd.md.integrate.mode_standard(dt=dt)

    # Set up integrator
    langevin = hoomd.md.integrate.langevin(group=moving_group, kT=T, seed=3991)

    #hoomd.update.box_resize(L=hoomd.variant.linear_interp([(0,system.box.Lx),
    #                        (resize_steps-500,box_length)]),scale_particles=True)
    # for i,name in enumerate(aa_type):
        # langevin.set_gamma(name, gamma=aa_mass[i]/1000.0)
    # langevin.set_gamma('R', gamma=1)

    hoomd.dump.gsd('rigid_FUS_start.gsd', period=1, group=all_group, truncate=True)
    #hoomd.analyze.log(filename="potential_ene.log", quantities=['potential_energy'], period=100, overwrite=False)
    hoomd.run(1000)
