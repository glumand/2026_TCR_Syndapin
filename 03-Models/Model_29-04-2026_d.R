# -------------------------------------------------------------------------#
# T cell receptor levels  ---- Modelling main()
# -------------------------------------------------------------------------#
#
# [PURPOSE]
# 
# Adaptation of Model for TCR processes. Describing the influence of Syndapin on TCR levels.
#
# Info:
#
# [AUTHOR]
# model by Mio Heinrich, based on a model by Simon Beyer and Viviane Timmermann
# data by Jakob Diesner
#
# [Date]
# Tue, 14 April 2026
# 
# 14-04-2026 - initial model
# 17-04-2026 - add steady state data / no ligand data
# 24-04-2026 - remove TCR-Endo feedback
# 27-04-2026 - based on 24-04-2026, begin reduction
# 29-04-2026 - leave out ER, currently keep feedback
# 29-04-2026_b - start reduction, set scale_Total to 1
# 29-04-2026_c - start reduction, TCR_ENDO_CTRL to 0
# 29-04-2026_d - start reduction, k_B_lig, k_P and TCR_Surface to 0

library(dplyr)
library(dMod)
library(ggplot2)
library(deSolve)
library(rootSolve)
library(tidyverse)
library(dplyrExtras)


rm(list = ls(all.names = TRUE))


# Create and set a specific working directory inside your project folder
.workingDir <- file.path(purrr::reduce(1:2, ~dirname(.x), .init = rstudioapi::getSourceEditorContext()$path), "wd")
if (!dir.exists(.workingDir)) dir.create(.workingDir)
setwd(.workingDir)

.modelname <- "Model_29-04-2026_d"
.build = T

# set a specific working directory and define folders
.originaldataFolder       <- normalizePath(file.path(.workingDir,"../00-OriginalData"))
.dataFolder               <- normalizePath(file.path(.workingDir,"../01-Data"))
.scriptFolder             <- normalizePath(file.path(.workingDir,"../02-Scripts"))
.modelFolder              <- normalizePath(file.path(.workingDir,"../03-Models"))
.resultsFolder            <- normalizePath(file.path(.workingDir,"../04-Results", .modelname))
.plotFolder               <- normalizePath(file.path(.workingDir,"../05-Plots", .modelname))

source(file.path(.scriptFolder, "Useful_functions_TCR.R"))


for (x in c(.originaldataFolder, .dataFolder, .scriptFolder, .modelFolder, .resultsFolder, .plotFolder))
  if (!dir.exists(x)) dir.create(x, recursive = TRUE)

unlink(list.files(full.names = TRUE), recursive = TRUE)

# # Data Preprocessing
# source(file.path(.scriptFolder, "S01-DataPreprocessing.R"))




# Plot Data
data <- read.csv(file.path(.dataFolder, "Data_15042026_formatted.csv")) %>% 
  rename(time = time_minutes) %>% 
  rename(name = observable) %>% 
  # choose only cellline shCTRL
  filter(cellline == "shCTRL") %>%
  # filter(ligand == 1) %>% 
  filter(!is.na(value)) %>% 
  # delete time points zero for ligand = 0
  filter(!(time %in% c(0, 120) & ligand == 0))

# add a datapoint at time_minutes = 300 for celline shCTRL, ligand 1 and total_TCR that copies the 240 timepoint
data <- data %>% 
  bind_rows(data %>% filter(time == 240, cellline == "shCTRL", ligand == 1, name == "total_TCR", replicate == "Probe_1") %>% mutate(time = 300)) %>% 
  bind_rows(data %>% filter(time == 360, cellline == "shCTRL", ligand == 1, name == "total_TCR", replicate == "Probe_2") %>% mutate(time = 300)) %>% 
  bind_rows(data %>% filter(time == 240, cellline == "shCTRL", ligand == 1, name == "total_TCR", replicate == "Probe_3") %>% mutate(time = 300))

# calculate standart deviation for same name, cellline, time, date and ligand for the three replicates. add as column sigma to all three
data <- data %>%
  group_by(name, cellline, time, date, ligand) %>%
  mutate(sigma = sd(value)) %>%
  ungroup() %>% 
  filter(!is.na(sigma)) %>% 
  mutate(sigma = sigma/(value*log(2))) %>% 
  mutate(value = log2(value))


# plot data, time_minutes on x, value on y, color by celline and ligand, facet by observable
# add sigma to plot
data_plot <- ggplot(data, aes(x = time, y = value, color = cellline)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = value - sigma, ymax = value + sigma), width = 0.2) +
  facet_wrap(name~ligand, scales = "free_y") +
  theme_bw() +
  labs(x = "Time (minutes)", y = "log2(Value)", color = "Cell Line") +
  theme(legend.position = "bottom")

print(data_plot)

# Define Observables
observables <- c(total_TCR = "log2(TCR_Endo + TCRL_Endo + TCR_Surface + TCRL_Surface + offset_Total) + scale_Total",
                 surface_TCR_sek = "log2(TCRL_Surface + offset_Surface_sec) + scale_Surface",
                 surface_TCR_full = "log2(TCRL_Surface + TCR_Surface + offset_Surface_full) + scale_Surface")


reactions <- NULL %>%
  
  # Transcription, Translation, Assembly and Degradation in ER
  # TCR creation is inhibited by TCR on the surface
  addReaction("", "TCR_Surface", rate = "k_P / (1 + k_IM_1 * (TCR_Surface + TCRL_Surface * ligand))") %>% 
  addReaction("TCR_Surface", "", rate = "k_B * TCR_Surface") %>%  
  addReaction("TCRL_Surface", "", rate = "k_B_lig * TCRL_Surface * ligand") %>%  
  
  
  # Internalization to Endosom and Recycling 
  addReaction("TCR_Surface", "TCR_Endo", rate = "k_I * TCR_Surface") %>%
  addReaction("TCR_Endo", "TCR_Surface", rate = "k_PM_4 * (TCR_Surface + TCRL_Surface*ligand) * TCR_Endo") %>%
  
  addReaction("TCRL_Surface", "TCRL_Endo", rate = "k_I_lig * TCRL_Surface * ligand") %>%
  addReaction("TCRL_Endo", "TCRL_Surface", rate = "k_PM_4_lig * (TCRL_Surface + TCR_Surface) * TCRL_Endo * ligand") %>%
  
  # Degradation in Endosom
  addReaction("TCR_Endo", "", rate = "k_D * TCR_Endo") %>% 
  
  addReaction("TCRL_Endo", "", rate = "k_D_lig * TCRL_Endo * ligand") %>% 
  
  
  
  {.}



# mysteadieslate <- steadyStates(r) produces negative results!
# create SS constraints
# Setting up SS constraint 
deqns <- as.eqnvec(reactions) %>% as.character()
names(deqns) <- paste0("d", reactions$states)


ssCnstrDatalist <- data.frame(time = rep(seq(0, 1000, length.out = 10), length(deqns)) %>% sort())
ssCnstrDatalist$name <- rep(names(deqns), 5)
ssCnstrDatalist$value <- 0
ssCnstrDatalist$sigma <- 1e-5
ssCnstrDatalist$condition <- "0_shCTRL"
ssCnstrDatalist <- as.datalist(ssCnstrDatalist)

datalist <- data %>% as.data.frame() %>% as.datalist(split.by = c("ligand", "cellline"))

cond.grid <- attr(datalist, "condition.grid")




innerpars <- c(getParameters(reactions), getSymbols(observables)) %>% unique()
trafo <- NULL %>% 
  define("x~x", x = innerpars) %>% 
  # insert("inner~steadyEqn", inner = names(mysteadies), steadyEqn = mysteadies) %>%
  branch(table = cond.grid) %>%
  insert("x~y", x="ligand", y=ligand) %>% 
  insert("x~exp(x * log(10))", x = .currentSymbols[!grepl("nhill|scale|ligand", .currentSymbols)]) %>%
  # set experementally given initial conditions
  # define("x~0", x="TCR_Surface", conditionMatch = "1_shCTRL") %>%
  define("x~0", x = "TCRL_Surface", conditionMatch = "0_shCTRL") %>%
  define("x~0", x = "TCRL_Endo") %>%
  
  # set condition specific inits
  # insert("x~y", x="TCR_Endo", y="TCR_Endo_CTRL", conditionMatch = "0_shCTRL") %>% 
  
  # set offset_Surface_full to 0
  define("x~0", x = "offset_Surface_full") %>%
  
  # set scale_Total to 1
  define("x~1", x = "scale_Total") %>%
  
  # set TCr_Endo_CTRL to 0
  define("x~0", x = "TCR_Endo", conditionMatch = "0_shCTRL") %>%
  
  # k_B_lig, k_P and TCR_Surface to 0
  define("x~0", x = "k_B_lig") %>%
  define("x~0", x = "k_P") %>%
  define("x~0", x = "TCR_Surface") %>%
  
  
  {.}

# fixed <- c("TCR_Golgi_1", "TCR_Golgi_2", "TCR_Endo", "TCR_Surface")
fixed <- c()

if (.build) {
  
  if (!dir.exists(file.path(.modelFolder, "wd"))) {
    dir.create(file.path(.modelFolder, "wd"), recursive = TRUE)
  }
  setwd(file.path(.modelFolder, "wd"))
  
  # Construct ODE Predictor
  x <- odemodel(reactions, modelname = "x", compile = F, fixed = fixed) %>% Xs()
  
  # Construct Prediction Function and Error Model function
  # g <- Y(c(observables, deqns), f = x, attach.input = T, compile = F, modelname = "g")
  g <- Y(c(observables, deqns), f = reactions, attach.input = T, compile = F, modelname = "g")
  
  # Construct p2p function
  p <- P(trafo, compile = F, modelname = "p")
  
  # and Compile
  compile(g, x, p, output = .modelname, cores = 6, verbose = F)
  save(x, g, p, file = paste0(.modelname, ".rds"))
  
  # Archive build artifacts before deleting
  build_files <- list.files(pattern = "(\\.o|\\.c|\\.cpp)$")
  if (length(build_files) > 0) {
    tar(paste0(.modelname, "_src.tar.gz"), files = build_files, compression = "gzip")
  }
  unlink(build_files)
  
  # Copy model files to workingDir
  model_files <- list.files(pattern = .modelname)
  if (length(model_files) > 0) {
    file.copy(from = model_files, to = .workingDir, overwrite = TRUE)
    message("All matching files successfully copied!")
  } else {
    stop("Error: No files containing '.modelname' were found in the current directory.")
  }
  setwd(.workingDir)
  .archive <- paste0(.modelname, "_src.tar.gz")
  untar(.archive, exdir = .workingDir)
  loadDLL(x)
  
} else {
  setwd(file.path(.modelFolder, "wd"))
  if (!file.exists(paste0(.modelname, ".rds")) || !file.exists(paste0(.modelname, ".so"))) {
    stop("Run the script again with .build = TRUE")
  }
  
  model_files <- list.files(pattern = .modelname)
  if (length(model_files) > 0) {
    file.copy(from = model_files, to = .workingDir, overwrite = TRUE)
    message("All matching files successfully copied!")
  } else {
    stop("Error: No files containing '.modelname' were found in the current directory.")
  }
  
  # Extract build artifacts archive into workingDir
  .archive <- paste0(.modelname, "_src.tar.gz")
  if (file.exists(.archive)) {
    untar(.archive, exdir = .workingDir)
    message("Build artifacts extracted.")
  } else {
    warning("No build artifacts archive found. Source files (.o, .c, .cpp) are not available.")
  }
  
  setwd(.workingDir)
  load(paste0(.modelname, ".rds"))
  loadDLL(x)
}

mytimes <- c(seq(0, max(data$time) *1.1, 1))
# mytimes <- c(seq(0, 200, 1))

outerpars <- getParameters(p)
ini <- structure(rep(-1, length(outerpars)), names = outerpars)
parlb <- structure(rep(-10, length(outerpars)), names = outerpars)
parub <- structure(rep(10, length(outerpars)), names = outerpars)
prd <- Reduce('*', list(g, x, p))

plotPrediction(prd(mytimes, ini))
plotCombined(prd(mytimes, ini), datalist)




obj_data <- normL2(datalist, g*x*p)
obj_prior <- constraintL2(ini, sigma = 10, attr.name = "prior")
obj_ssconstr <- normL2(ssCnstrDatalist, g*x*p, attr.name = "ssCnstr")

obj <- Reduce("+", list(obj_data, obj_prior, obj_ssconstr))
# obj <- Reduce("+", list(obj_data, obj_prior))
obj(ini)

# out_mstrust <- mstrust(
#   objfun=obj,
#   center=ini,
#   studyname = "multistart",
#   rinit = 0.1, 
#   rmax = 8, 
#   iterlim=1000, #
#   sd = 3,
#   # parupper = upperB, 
#   # parlower = lowerB, 
#   fits = 10, 
#   cores = 4, 
#   fixed = NULL,
#   # printIter = T,
#   # output = T
#   cautiousMode = T
# )
# myframe <- as.parframe(out_mstrust)
# bestfit <- as.parvec(myframe, 1)
# plotPrediction(prd(mytimes, bestfit))
# plotCombined(prd(mytimes, bestfit), datalist)

fit.name <- "2904"

outFitHelix <- distributed_computing(
  {
    mstrust(
      objfun=obj,
      center=ini,
      studyname = "multistart",
      rinit = 0.1,
      rmax = 8,
      iterlim=5000, #
      sd = 3,
      parupper = parub,
      parlower = parlb,
      fits = 1,
      cores = 16,
      fixed = NULL,
      # printIter = T,
      # output = T
      cautiousMode = T
    )
  },
  jobname = fit.name,
  partition = "cpu-single",
  cores = 1,
  nodes = 1,
  mem_per_core = 2,
  walltime = "00:30:00",
  # ssh_passwd = "password",
  machine = "helix",
  var_values = NULL,
  no_rep = 1000,
  recover = F,
  compile = F
)

outFitHelix$check()
outFitHelix$get()
# print(cluster_result)

# outFitHelix$purge()
flat_results <- do.call(c, cluster_result)
my_parlist <- as.parlist(flat_results)
outframe <- as.parframe(my_parlist)
outframe$value <- outframe$value + 848
waterfall <- plotValues(outframe[1:100], tol=0.01)
print(waterfall)
# save waterfall to plot folder
ggsave(paste0(.plotFolder, "/waterfall.pdf"), plot = waterfall, width = 15, height = 11)

bestfit <- as.parvec(outframe, 1)
fitvalue <- outframe$value[1]
plotcombined <- plotCombined(prd(mytimes, bestfit), datalist)
print(plotcombined)
ggsave(paste0(.plotFolder, "/plotcombined.pdf"), plot = plotcombined, width = 15, height = 11)


plottimes <- c(seq(0, max(data$time) *1.1, 1))

full_plotcombined <- plotCombined(prd(plottimes, bestfit), datalist)
print(full_plotcombined)
ggsave(paste0(.plotFolder, "/full_plotcombined.pdf"), plot = full_plotcombined, width = 15, height = 11)

saveRDS(outframe, file = file.path(.resultsFolder, paste0("24042026_outframe.rds")))



outHelix_profs <- SendProfilesHelix(obj, 
                                    bestfit, 
                                    names(bestfit), 
                                    jobname= paste0("profiles", .modelname), 
                                    argsDC=list(recover = F, link = T, walltime="00:50:00", cores = 1),
                                    method = "integrate",
                                    stepControl = list(stepsize = 1e-4, min = 1e-4, max = Inf, atol = 1e-2, rtol = 1e-2, limit = 1e2, stop = "data"),
                                    algoControl = list(reoptimize = T),
                                    optControl = list(rinit = 0.1, rmax = 1, iterlim = 1e3, fterm = 1e-5, mterm = 1e-5))

# outHelix_profs$check()
profiles <- outHelix_profs$get()
plotProfile(profiles, mode=="data")
profile_plot <- plotProfile(profiles, mode==c("data", "ssCnstr"))
print(profile_plot)
ggsave(paste0(.plotFolder, "/profiles.pdf"), plot = profile_plot, width = 15, height = 11)

save(file = file.path(.resultsFolder, "profile_results.Rdata"), profiles)
plotProfilesAndPaths(profiles, whichpars = c("scale_Surface", "scale_Total", "offset_Surface_sec", "offset_Total"))
plotProfilesAndPaths(profiles, whichpars = c("k_B", "k_B_lig", "k_P", "TCR_Surface"))
plotProfilesAndPaths(profiles, whichpars = c("TCR_Surface", "TCRL_Surface", "TCR_Endo", "TCR_Endo_CTRL"))
plotProfilesAndPaths(profiles, whichpars = c("k_I", "k_IM_1"))
