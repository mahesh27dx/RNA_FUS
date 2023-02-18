# -*- coding: utf-8 -*-
#!/usr/bin/python3

# this script starts an MD simulation of multiple chains of a chosen sequence of the protein hnRNPA1
# short-range pairwise potential: HPS model (accounts for both protein-protein and protein-splvent interactions)

import math
import numpy as np
import hoomd
import hoomd.md
import sys
import warnings

try:
    from hoomd import azplugins
except ImportError:
    import azplugins

def initSingleChain(Lx, Ly, Lz, monomerSizesInSeq, snapshot):
    """
    This function initializes the positions and bonds of a single straight chain
    in a box of size [Lx, Ly, Lz], with Lz > Lx=Ly

    monomerSizesInSeq contains the sizes of the monomers of the chain
    -> (initial distance between monomer i and i+1) = (size_i, size_i+1)/2
     positions are saved in snapshot (xpos = ypos = 0; zpos centered)
    """

    chainLen = len(monomerSizesInSeq)
    totalLen = sum(monomerSizesInSeq)
    print(f"initialize chain of total length {totalLen}, in Box with Lz = {Lz}")
    if totalLen >= Lz:
        sys.exit('total length of the chain is bigger than Lz')
    if totalLen >= Lx:
        warnings.warn('total length of the chain is bigger than Lx, so the chain might interact with itself')

    zpos = np.zeros(chainLen)
    zpos[0] = -Lz / 2.0 + (Lz - totalLen) / 2.0 + monomerSizesInSeq[0] / 2.0
    for i in range(chainLen - 1):
        zpos[i+1] = zpos[i] + (monomerSizesInSeq[i] + monomerSizesInSeq[i+1]) / 2.0
    snapshot.particles.position[:, 2] = zpos
    bond1 = np.linspace(0, chainLen-2, chainLen-1).astype(int)
    snapshot.bonds.group[:] = np.array([bond1, bond1 + 1]).transpose()

def get_HPSparameters():
    """
    In HOOMD units!
    aminoAcids = [ALA, ARG, ASN, ASP, CYS, GLN, GLU, GLY, HIS, ILE, LEU, LYS,
                    MET, PHE, PRO, SER, THR, TRP, TYR, VAL]
    """
    aminoAcids = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M',
                    'F', 'P', 'S', 'T', 'W', 'Y', 'V']
    lambdaAcids = np.array(
        [0.72973, 0.0, 0.432432, 0.378378, 0.594595, 0.513514, 0.459459, 0.648649, 0.513514, 0.972973,
         0.972973, 0.513514, 0.837838, 1.0, 1.0, 0.594595, 0.675676, 0.945946, 0.864865, 0.891892])
    massAcids = (1 / 57.05) * np.array(
        [71.08, 156.20, 114.10, 115.10, 103.10, 128.10, 129.10, 57.05, 137.10, 113.20, 113.20,
         128.20, 131.20, 147.20, 97.12, 87.08, 101.10, 186.20, 163.20, 99.07])
    sigmaAcids = (1 / 4.5) * np.array(
        [5.04, 6.56, 5.68, 5.58, 5.48, 6.02, 5.92, 4.50, 6.08, 6.18, 6.18, 6.36, 6.18, 6.36, 5.56,
         5.18, 5.62, 6.78, 6.46, 5.86])
    chargeAcids = np.array([0, 1, 0, -1, 0, 0, -1, 0, 0.5, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0])  # in units: e

    return aminoAcids, lambdaAcids, massAcids, sigmaAcids, chargeAcids

#################################################################
######## SET PARAMETERS
#############################################################

protein = 'FUS_LC'
Npoly = 350
proteinStr = np.genfromtxt('input_files/' + protein + '.dat', dtype=str, delimiter=1)
chainLen = len(proteinStr)  # take the whole sequence
print(f"The total length of the cahin is {chainLen}")

model = 'HPS'

saveTrajectories = 1    # 1: save trajectories in gsd file, 0: don't save trajectories
saveThermodynamics = 1  # 1: save temperature and pot/kin energy during the simulation steps
thermoPeriod = 500  # defines the period of the measurements of the thermodynamic quantites
numberGSDframes = 1000  # the number of sim writes into the gsd file
initFromSim = 0 # 1: load snapshot from another sim for initial state; 0; init snapshot from scratch
continueFromSim = 0 # 1: load snapshot from another sim for initial state; 0: init snapshot from scratch

dt = 2.69e-3    # time step (in HOOMD units; equals dt = 1e-14 s)
steps = int(500e6)  # number of simulation steps
T_in_K = float(sys.argv[1]) # temperature - held constant by a thermostat; input in kelvin

Nmonom = chainLen * Npoly   # total number of monomers
Nbonds = Npoly * (chainLen - 1) # total number of bonds

Lx_in_nm = 15   # box size x-direction in nm (1D = 4.5 Angs -> 15nm = 33.333 D)
Ly_in_nm = 15   # box size y-direction in nm (1D = 4.5 Angs -> 15nm = 33.333 D)
Lz_in_nm = 150   # box size z-direction in nm (1D = 4.5 Angs -> 75nm = 166.667 D)

rcut_AH_in_nm = 2.0 # cutoff radius of short-range interactions in nm (Ashbough Hatch functional form)
rcut_YU_in_nm = 4.0 # cutoff radius of electrostatic interactions in nm (Coulombic term; HOOMD: Yukawa potential)

# in HOOMD units
epsilon = 1.0   # absolute energy scale of the short-ranged interactions between all pairs of amino acids (see lj_AH, get value from PLOS Articles) (1 epsilon = 0.2kcal/mol)
r0 = 0.38 / 0.45    # bond equilibrium position (see harmonic potential, r0 = 3.8A, 1D = 4.5A)
k = 2025.0  # spring constant (see harmonic potential, k=20J/m^2), high=> stiff: bondlength constrained to be almost exactly sigma
D = 80.0    # dielectric constant of the solvent (water) (see Yukawa potential)
kappa = 1/2.222 # inverse (due to def of Yukawa) of Debye screening length (~1nm)

seed1 = 1

sout = protein + '_M' + str(chainLen) + '_NP' + str(Npoly) + '_T' + str(int(T_in_K)) + '_Box' + str(Lx_in_nm) +\
                '_' + str(Ly_in_nm) + '_' + str(Lz_in_nm) + '_s' + str(steps)

fromFile = 'MDhnRNPA1seq' + sout + '.gsd'

# convert temperature and boxsize into HOOMD-units:
k_B = 1.38064852e-23                 # boltzmann constant [J/K]
N_A = 6.02214076e23                  # Avogadro constant [1/mol]
epsilon_in_kcal_mol = 0.2            # energy unit [kcal/mol]: 1epsilon = 0.2 kcal/mol
distUnit_in_A = 4.5                  # distance unit [A]: 1D = 4.5A

T = round(T_in_K * k_B / (epsilon_in_kcal_mol * 4184 / N_A), 3)          # 1kcal = 4184J
Lx = round(10.0 * Lx_in_nm / distUnit_in_A, 3)
Ly = round(10.0 * Ly_in_nm / distUnit_in_A, 3)
Lz = round(10.0 * Lz_in_nm / distUnit_in_A, 3)

rcut_YU = round(10.0 * rcut_YU_in_nm / distUnit_in_A, 3)
rcut_AH = round(10.0 * rcut_AH_in_nm / distUnit_in_A, 3)

##################################################################
################# Define the protein
#######################################################
if model == 'HPS':
    aminoAcids, lambdaAcids, massAcids, sigmaAcids, chargeAcids = get_HPSparameters()
else:
    sys.exit('chosen model (parameters) not valid')

# calculate interaction parameter between all pairs of amino acids
numAcids = len(aminoAcids)
lambdaHPS = np.zeros([numAcids, numAcids])
sigmaHPS = np.zeros([numAcids, numAcids])
for i in range(numAcids):
    for j in range(numAcids):
        lambdaHPS[i, j] = ([lambdaAcids[i] + lambdaAcids[j]) / 2    # (mean) hydrophobicity values for each amino acid interaction pair
        sigmaHPS[i, j] = (sigmaAcids[i] + sigmaAcids[j]) / 2    # (mean) sizes for each amino acid interaction pair

proteinStrID = np.zeros(chainLen, int)
massSeq = np.zeros(chainLen)
chargeSeq = np.zeros(chainLen)
sigmaSeq = np.zeros(chainLen)
for i in range(chainLen):
    # try:
        # index = aminoAcids.index(proteinStr[i])
    # except:
        # index = np.where(aminoAcids == proteinStr[i])[0][0]
    index = aminoAcids.index(proteinStr[i])
    proteinStrID[i] = index # find the index of each acid in the given sequence
    massSeq[i] = massAcids[index]   # save masses of the residues according to given sequence
    chargeSeq[i] = chargeAcids[index]   # save charges of residues according to given sequence
    sigmaSeq[i] = sigmaAcids[index] # save sizes of residues according to given sequence

#####################################
######## HOOMD INITIALIZATION
##################################

context = hoomd.context.initialize("")

# Initialize polymers
if initFromSim:         # INITIALIZE SYSTEM:positions and bonds from snapshot (different simulation)
    snapshot = hoomd.data.gsd_snapshot(filename=fromFile, frame=-1)
    system = hoomd.init.read_snapshot(snapshot)
elif continueFromSim:   # INITIALIZE SYSTEM:positions and bonds from snapshot (different simulation)
    system = hoom.init.read_gsd(filename=fromFile, restart=fromFile, frame=-1)
    # stepsForGSDdump = steps
    steps = int(steps + continueGSDdumpPeriod)
else:   # INITIALIZE SYSTEM: positions and bonds from snapshot (snapshot initialised individually)
    snapshot = hoomd.data.make_snapshot(N=Nmonom, box=hoomd.data.boxdim(Lx=Lx, Ly=Ly, Lz=Lz), particle_types=aminoAcids,
                                        bond_types=['polymer'])
    snapshot.bonds.resize(Nbonds)

    # initialise positions and bonds
    if Npoly > 1:
        if protein == 'FUS':
            initManyChains_FUS(Lx, Ly, Lz, Npoly, sigmaSeq, snapshot, r0)
        else:
            initManyChains(Lx, Ly, Lz, Npoly, sigmaSeq, snapshot)
    else:
        initSingleChain(Lx, Ly, Lz, sigmaSeq, snapshot)

    # define the sequence of amino acids (all chains)
    snapshot.particle.typeid[:] = np.tile(proteinStrID, Npoly)
    snapshot.particles.mass[:] = np.tile(massSeq, Npoly)
    snapshot.particles.charge[:] = np.tile(chargeSeq, Npoly)
    snapshot.particles.diameter[:] = np.tile(sigmaSeq, Npoly)
    system = hoomd.init.read_snapshot(snapshot)


# electrostatic interactions (Yukawa) and short-range pairwise interactions (HPS model)
nl = hoomd.md.nlist.cell()  # use cell list for force calculattion
nl.reset_exclusions(exclusions=['1-2']) # np pair interactions between bonded particles (interact only via bond interactions)

np = azplugins.pair.ashbaugh(r_cut = 0, nlist=nl)   # pairwise interactions
yukawa = hoomd.md.pair.yukawa(r_cut=rcut, nlist=nl) # Yukawa potential for electrostatic interactions
for a1 in range(numAcids):
    for a2 in range(numAcids):
        nb.pair_coeff.set(aminoAcids[a1], aminoAcids[a2], lam=lambdaHPS[a1, a2], epsilon=epsilon, sigma=sigmaHPS[a1,a2], r_cut=rcut_AH)
        if chargeAcids[a1] == 0 or chargeAcids[a2] == 0:    # skip pair, because potential is zero (saves time)
            yukawa.pair_coeff.set(aminoAcids[a1], aminoAcids[a2], epsilon=chargeAcids[a1]*chargeAcids[a2] * 4.61008,
            kappa=kappa, r_cut=False) # get prefactor in epsilon from
        else:
            yukawa.pair_coeff.set(aminoAcids[a1], aminoAcids[a2], epsilon=chargeAcids[a1] * chargeAcids[a2] * 4.61008,
            kappa=kappa)

# harmonic potential between bonded particles
if chainLen > 1:
    harmonic = hoomd.md.bond.harmonic()
    harmonic.bond_coeff.set('polymer', k=k, r0=r0)

# integrate at constant temperature
all = hoomd.group.all()
hoomd.md.integrate.mode_standard(dt=dt)
integrator = hoomd.md.integrate.langevin(group=all, kT=T, seed=seed1) # use Langevin dynamics
