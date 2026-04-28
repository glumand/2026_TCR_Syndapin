#!/bin/bash

# Job name
#SBATCH --job-name=2704
# Define format of output, deactivated
#SBATCH --output=2704_%j-%a.out
# Define format of errorfile, deactivated
#SBATCH --error=2704_%j-%a.err
# Define partition
#SBATCH --partition=cpu-single
# Define number of nodes per task
#SBATCH --nodes=1
# Define number of cores per node
#SBATCH --ntasks-per-node=1
# Define walltime
#SBATCH --time=00:30:00
# Define of repetition
#SBATCH -a 0-999
# memory per CPU core
#SBATCH --mem-per-cpu=2gb


# Load R modules
module load math/R
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

# Run R script
Rscript 2704.R