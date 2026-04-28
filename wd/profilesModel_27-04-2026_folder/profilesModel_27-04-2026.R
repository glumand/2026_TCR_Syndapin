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
load('profilesModel_27-04-2026_workspace.RData')

# remove random seeds
rm(.Random.seed)
set.seed(as.numeric(Sys.getenv('SLURM_JOB_ID')))

# load shared object if precompiled


files <- list.files(pattern = '.so$')
for (f in files) dyn.load(f)

# List of variablevalues
var_values_1=c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11,12,12,13,13,14,14)
var_values_2=c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11,12,12,13,13,14,14)
var_values_3=c('left','right','left','right','left','right','left','right','left','right','left','right','left','right','left','right','left','right','left','right','left','right','left','right','left','right','left','right')

# Define variable values per run
var_1=var_values_1[(as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID')) + 1)]
var_2=var_values_2[(as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID')) + 1)]
var_3=var_values_3[(as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID')) + 1)]

# Fixed parameters
node_ID = Sys.getenv('SLURM_ARRAY_TASK_ID')
job_ID = Sys.getenv('SLURM_JOB_ID')
jobname = 'profilesModel_27-04-2026'



# Paste function call
cluster_result <- try({
    pargs <- profargs
    pargs$whichPar <- pargs$whichPar[var_1:var_2]
    if (!is.null(var_3)) 
        pargs$side <- var_3
    do.call(dMod::profile, pargs)
})
save(cluster_result, file = paste0(jobname,'_', node_ID, '_result.RData'))