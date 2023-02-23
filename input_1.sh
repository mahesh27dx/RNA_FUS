#!/bin/bash

#SBATCH -J FUS      # Job name
#SBATCH -o fus_%j.out       # Specify stdout output file (%j expands to jobId)
#SBATCH -p m2_gpu           # Partition name
#SBATCH -c 2
#SBATCH -n 1                # Total number of tasks
#SBATCH -t 120:00:00         # Run time (hh:mm:ss)

#SBATCH -A m2_komet331hpc   # Specify allocation to charge against
#SBATCH --mem=8G           # Memory allocation

#mkdir output_10ns

timestart=$(date +"%s")
#### LOAD MODULES ########
# setpkgs -a phys/HOOMD/2.9.6-fosscuda-2020b-double
. /home/yadavmah/binaries/env/bin/activate
module load system/CUDA/11.1.1-GCC-10.2.0
#module load phys/HOOMD/2.9.6-fosscuda-2020b-single


# Execution

#for i in 190 200
#do
srun python3 initial_rigid_snapshot.py  #input/FUS_LC.dat $i
#done

timeend=$(date +"%s")

iduration=$(awk "BEGIN {print ($timeend - $timestart) / 3600}")

echo "Simulation took $duration hours to complete"

#### Following commands can be used #####

# chmod ugo+x a.out
#mv a.out ~/binf77/ising_2d_simulation
