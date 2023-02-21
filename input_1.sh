#!/bin/bash

#SBATCH -J FUS_1_1ms      # Job name
#SBATCH -o fus_%j.out       # Specify stdout output file (%j expands to jobId)
#SBATCH -p m2_gpu           # Partition name
#SBATCH -c 2
#SBATCH -n 1                # Total number of tasks
#SBATCH -t 120:00:00         # Run time (hh:mm:ss)
        # Reserve 1 GPUs
#SBATCH -A m2_komet331hpc   # Specify allocation to charge against
#SBATCH --mem=2G           # Memory allocation

#mkdir output_10ns

timestart=$(date +"%s")
#### LOAD MODULES ########
# setpkgs -a phys/HOOMD/2.9.6-fosscuda-2020b-double
module load system/CUDA/11.1.1-GCC-10.2.0
#module load phys/HOOMD/2.9.6-fosscuda-2020b-single
. /home/yadavmah/env/bin/activate

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

# lang/Python/3.10.8-GCCcore-12.2.0-bare
#  phys/HOOMD/2.9.6-fosscuda-2020b-single
# /cluster/easybuild/broadwell/software/HOOMD/2.9.6-fosscuda-2020b-single/hoomd
# phys/HOOMD/2.9.6-fosscuda-2020b-double
# phys/HOOMD/2.9.6-fosscuda-2019b-single
# phys/HOOMD/2.9.1-fosscuda-2019b-double
   # phys/HOOMD/2.9.1-fosscuda-2019b-single
