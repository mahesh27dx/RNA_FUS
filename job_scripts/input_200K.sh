#!/bin/bash

#SBATCH -J FUS_200K      # Job name
#SBATCH -o fus_200K_%j.out       # Specify stdout output file (%j expands to jobId)
#SBATCH -p m2_gpu           # Partition name
#SBATCH -c 16
#SBATCH -n 1                # Total number of tasks
#SBATCH -t 118:00:00         # Run time (hh:mm:ss)
#SBATCH -A m2_komet331hpc   # Specify allocation to charge against
#SBATCH --mem=20G           # Memory allocation

timestart=$(date +"%s")
#### LOAD MODULES ########
# setpkgs -a phys/HOOMD/2.9.6-fosscuda-2020b-double
#. /home/yadavmah/binaries/env/bin/activate
module load system/CUDA/11.1.1-GCC-10.2.0

srun python3 simulation_run.py -dt="0.001" -time="10000000" -temp="200"  #input/FUS_LC.dat $i


timeend=$(date +"%s")

iduration=$(awk "BEGIN {print ($timeend - $timestart) / 3600}")

echo "Simulation took $duration hours to complete"

#### Following commands can be used #####

# chmod ugo+x a.out
#mv a.out ~/binf77/ising_2d_simulation

# chmod ugo+x a.out
#mv a.out ~/binf77/ising_2d_simulation
