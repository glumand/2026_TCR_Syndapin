#!/bin/bash

# Job name
#SBATCH --job-name=Model_14-04-2026_fit1
# Define format of output, deactivated
#SBATCH --output=Model_14-04-2026_fit1_%j-%a.out
# Define format of errorfile, deactivated
#SBATCH --error=Model_14-04-2026_fit1_%j-%a.err
# Define partition
#SBATCH --partition=cpu-single
# Define number of nodes per task
#SBATCH --nodes=1
# Define number of cores per node
#SBATCH --ntasks-per-node=4
# Define walltime
#SBATCH --time=1:15:00
# Define of repetition
#SBATCH -a 0-0
# memory per CPU core
#SBATCH --mem-per-cpu=2gb


# Load R modules
module load math/R
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

# Run R script
Rscript Model_14-04-2026_fit1.R