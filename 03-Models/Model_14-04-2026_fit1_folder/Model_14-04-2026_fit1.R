#!/usr/bin/env Rscript

# Load packages
try(library(dplyrExtras))
try(library(dtplyr))
try(library(data.table))
try(library(lubridate))
try(library(forcats))
try(library(stringr))
try(library(purrr))
try(library(readr))
try(library(tidyr))
try(library(tibble))
try(library(tidyverse))
try(library(rootSolve))
try(library(deSolve))
try(library(grid))
try(library(ggplot2))
try(library(dMod))
try(library(cOde))
try(library(dplyr))
try(library(stats))
try(library(graphics))
try(library(grDevices))
try(library(utils))
try(library(datasets))
try(library(methods))
try(library(base))
try(library(tidyverse))

# Load environment
load('Model_14-04-2026_fit1_workspace.RData')

# remove random seeds
rm(.Random.seed)
set.seed(as.numeric(Sys.getenv('SLURM_JOB_ID')))

# load shared object if precompiled


files <- list.files(pattern = '.so$')
for (f in files) dyn.load(f)

# List of variablevalues


# Define variable values per run


# Fixed parameters
node_ID = Sys.getenv('SLURM_ARRAY_TASK_ID')
job_ID = Sys.getenv('SLURM_JOB_ID')
jobname = 'Model_14-04-2026_fit1'



# Paste function call
cluster_result <- try({
    mstrust(objfun = obj, center = ini, studyname = jobname, nodeID = node_ID, rinit = 0.1, rmax = 8, iterlim = 10000, sd = 3, parupper = upperB, parlower = lowerB, fits = 1, cores = 1, fixed = NULL, cautiousMode = T)
})
save(cluster_result, file = paste0(jobname,'_', node_ID, '_result.RData'))