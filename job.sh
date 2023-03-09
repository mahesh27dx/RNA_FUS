#!/bin/bash

timestart=$(date +"%s")

#python3 simulation_run.py -dt="0.001" -time="10000" -temp="300" -period="100"
python3 rna_sim_run.py -dt="0.001" -time="100000" -temp="310" -period="100"


timeend=$(date +"%s")

iduration=$(awk "BEGIN {print ($timeend - $timestart) / 3600}")

echo "Simulation took $duration hours to complete"

#### Following commands can be used #####

# chmod ugo+x a.out
#mv a.out ~/binf77/ising_2d_simulation

# chmod ugo+x a.out
#mv a.out ~/binf77/ising_2d_simulation
