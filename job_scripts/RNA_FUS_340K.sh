#!/bin/bash

#SBATCH -J RNA_340K_FUS      # Job name
#SBATCH -o rna_340K_fus_%j.out       # Specify stdout output file (%j expands to jobId)
#SBATCH -p m2_gpu           # Partition name
#SBATCH -n 1               # Total number of tasks
#SBATCH -c 4
#SBATCH -t 120:00:00         # Run time (hh:mm:ss)
#SBATCH --gres=gpu:1
#SBATCH -A m2_komet331hpc   # Specify allocation to charge against
#SBATCH --mem 10GB           # Memory allocation

timestart=$(date +"%s")
#### LOAD MODULES ########

module load phys/HOOMD/2.9.6-fosscuda-2020b-single
# module load system/CUDA/11.1.1-GCC-10.2.0
#. /home/yadavmah/binaries/env/bin/activate

srun python3 rna_sim_run.py -dt="0.01" -time="5000000" -temp="340" -period="10000"

timeend=$(date +"%s")

iduration=$(awk "BEGIN {print ($timeend - $timestart) / 3600}")

echo "Simulation took $duration hours to complete"

#### Following commands can be used #####

# chmod ugo+x a.out
#mv a.out ~/binf77/ising_2d_simulation

# chmod ugo+x a.out
#mv a.out ~/binf77/ising_2d_simulation
