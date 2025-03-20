#!/bin/bash -l
# The -l above is required to get the full environment with modules

# Set the allocation to be charged for this job
# not required if you have set a default allocation
#SBATCH -A 2020-1-30 


# The name of the script is myjob
#SBATCH -J Si320.sh

# Only 1 hour wall-clock time will be given to this job
#SBATCH -t 24:00:00

# Number of nodes
#SBATCH -N 1
# Number of MPI processes per node (the following is actually the default)
#SBATCH --ntasks-per-node=32
# Number of MPI processes.
#SBATCH -n 32

#SBATCH -e sb_error.e
#SBATCH -o sb_output.o

# Run the executable named myexe 
# and write the output into my_output_file
srun -n 32 ./inho > out01.o 2>&1
