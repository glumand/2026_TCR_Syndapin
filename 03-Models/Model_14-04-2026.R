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

.modelname <- "Model_14-04-2026"
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
  filter(ligand == 1) %>% 
  filter(!is.na(value))

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
observables <- c(total_TCR = "log2(TCR_ER + TCR_Endo + TCRL_Endo + TCR_Golgi_1 + TCR_Golgi_2 + TCR_Surface + TCRL_Surface + offset_Total) + scale_Total",
                 surface_TCR_sek = "log2(TCRL_Surface + offset_Surface_sec) + scale_Surface_sec",
                 surface_TCR_full = "log2(TCRL_Surface + TCR_Surface + offset_Surface_full) + scale_Surface_full")


reactions <- NULL %>%
  
  # Transcription, Translation, Assembly and Degradation in ER
  # TCR creation is inhibited by TCR on the surface
  addReaction("", "TCR_ER", rate = "k_P / (1 + k_IM_1 * (TCR_Surface + TCRL_Surface))") %>% 
  addReaction("TCR_ER", "", rate = "k_B * TCR_ER") %>%  
  
  
  # Export and transport from ER to Cell Surface
  # linear chain for time delay, export is inhibited by TCR on the surface
  addReaction("TCR_ER", "TCR_Golgi_1", rate = "k_E * TCR_ER / (1 + k_IM_3 * (TCR_Surface + TCRL_Surface))") %>% 
  addReaction("TCR_Golgi_1", "TCR_Golgi_2", rate = "k_E * TCR_Golgi_1") %>% 
  addReaction("TCR_Golgi_2", "TCR_Surface", rate = "k_E * TCR_Golgi_2") %>% 
  
  # Internalization to Endosom and Recycling 
  addReaction("TCR_Surface", "TCR_Endo", rate = "k_I * TCR_Surface") %>%
  addReaction("TCR_Endo", "TCR_Surface", rate = "k_PM_4 * (TCR_Surface + TCRL_Surface) * TCR_Endo") %>% 
  
  addReaction("TCRL_Surface", "TCRL_Endo", rate = "k_I_ligand * TCRL_Surface") %>%
  addReaction("TCRL_Endo", "TCRL_Surface", rate = "k_PM_4_ligand * (TCRL_Surface + TCR_Surface) * TCRL_Endo") %>% 
  
  # Degradation in Endosom
  addReaction("TCR_Endo", "", rate = "k_D * TCR_Endo") %>% 
  
  addReaction("TCRL_Endo", "", rate = "k_D_ligand * TCRL_Endo") %>% 
  
  {.}

# mysteadieslate <- steadyStates(r) produces negative results!

datalist <- data %>% as.data.frame() %>% as.datalist(split.by = c("ligand", "cellline"))

cond.grid <- attr(datalist, "condition.grid")




# Setting up SS constraint 
# mysteadies <- c(k_B = "k_P/TCR_ER")
# 
# deqns <- as.eqnvec(r) %>% as.character()
# names(deqns) <- paste0("d", r$states)
# 
# 
# ssCnstrDatalist <- data.frame(time = rep(seq(6000, 20160, length.out = 20), length(deqns)) %>% sort())
# ssCnstrDatalist$name <- rep(names(deqns), 5)
# ssCnstrDatalist$value <- 0
# ssCnstrDatalist$sigma <- 1e-5
# ssCnstrDatalist$condition <- "OPEN"
# ssCnstrDatalist <- as.datalist(ssCnstrDatalist)




innerpars <- c(getParameters(reactions), getSymbols(observables)) %>% unique()
trafo <- NULL %>% 
  define("x~x", x = innerpars) %>% 
  # insert("inner~steadyEqn", inner = names(mysteadies), steadyEqn = mysteadies) %>%
  branch(table = cond.grid) %>%
  insert("x~exp(x * log(10))", x = .currentSymbols[!grepl("nhill|scale", .currentSymbols)]) %>%
  define("x~0", x="TCR_Surface") %>%
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
  g <- Y(observables, f = reactions, attach.input = T, compile = F, modelname = "g")
  
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

outerpars <- getParameters(p)
ini <- structure(rep(-1, length(outerpars)), names = outerpars)
parlb <- structure(rep(-10, length(outerpars)), names = outerpars)
parub <- structure(rep(10, length(outerpars)), names = outerpars)
prd <- Reduce('*', list(g, x, p))

plotPrediction(prd(mytimes, ini))
plotCombined(prd(mytimes, ini), datalist)




obj_data <- normL2(datalist, g*x*p)
obj_prior <- constraintL2(ini, sigma = 5, attr.name = "prior")
# obj_ssconstr <- normL2(ssCnstrDatalist, g*x*p, attr.name = "ssCnstr")

# obj <- Reduce("+", list(obj_data, obj_prior, obj_ssconstr))
obj <- Reduce("+", list(obj_data, obj_prior))
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

fit.name <- "1604"

outFitHelix <- SendFitsHelix(
  obj,
  ini,
  jobname = paste0(fit.name,.modelname),
  argsDC = list(cores = 1, no_rep = 1000, recover = F, link = T),
  iterlim = 5e3,
  sd = 4,
  parlower = parlb,
  parupper = parub,
  fterm = 1e-5,
  mterm = 1e-5
)
# #
# outFitHelix$check()

results <- outFitHelix$get()

# save(file = file.path(.resultsFolder, "Fit_results.Rdata"), results)
# load(file.path(.resultsFolder, "Fit_results.Rdata"))
class(results) <- "parlist"
outframe <- as.parframe(results)
plotValues(outframe[1:500])

bestfit <- as.parvec(outframe, 1)
fitvalue <- outframe$value[1]
plotCombined(prd(mytimes, bestfit), datalist)

bestfit <- as.parvec(outframe)