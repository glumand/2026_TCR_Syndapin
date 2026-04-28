# -------------------------------------------------------------------------#
# T cell receptor levels  ---- Modelling main()
# -------------------------------------------------------------------------#
#
# [PURPOSE]
# 
# Modelling of TCR levels: expression, assembly & transport system
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
rm(list = ls(all.names = TRUE))
# Create and set a specific working directory inside your project folder
.workingDir <- file.path(purrr::reduce(1:2, ~dirname(.x), .init = rstudioapi::getSourceEditorContext()$path), "wd")
if (!dir.exists(.workingDir)) dir.create(.workingDir)
setwd(.workingDir)

.modelname <- "FinalModel"
.build = T

# set a specific working directory and define folders
.originaldataFolder       <- normalizePath(file.path(.workingDir,"../00-OriginalData"))
.dataFolder               <- normalizePath(file.path(.workingDir,"../01-Data"))
.scriptFolder             <- normalizePath(file.path(.workingDir,"../02-Scripts"))
.modelFolder              <- normalizePath(file.path(.workingDir,"../03-Models"))
.resultsFolder            <- normalizePath(file.path(.workingDir,"../04-Results", .modelname))
.plotFolder               <- normalizePath(file.path(.workingDir,"../05-Plots", .modelname))


for (x in c(.originaldataFolder, .dataFolder, .scriptFolder, .modelFolder, .resultsFolder, .plotFolder))
  if (!dir.exists(x)) dir.create(x, recursive = TRUE)

unlink(list.files(full.names = TRUE), recursive = TRUE)

# # Data Preprocessing
# source(file.path(.scriptFolder, "S01-DataPreprocessing.R"))

library(openxlsx)
library(dplyr)
library(dMod)


# Get Plotting function
source(file.path(.scriptFolder, "S02-Plotting.R"))



# Plot Data
data <- read.xlsx(file.path(.dataFolder, "aligned_data.xlsx")) %>% 
  filter(condition != "OPEN_PURO")

# custom_cond_names <- c(CTRL = "w/o", OPEN = "with biotin") # Custom condition names
custom_cond_names <- c(OPEN = "with biotin") # Custom condition names
custom_names <- c(Total = "[TCR] total", Surface = "[TCR] surface") # Custom names
plot_data <- plotdata(subset(data, condition == "OPEN"), 
                      custom.condition.names = custom_cond_names, 
                      custom.names = custom_names, 
                      split.cond = F, 
                      show.legend = FALSE,
                      base_size = 18)
plot_data

r <- NULL %>%
  
  # Transcription, Translation, Assembly and Degradation in ER
  addReaction("", "TCR_ER", rate = "k_P / (1 + k_IM_1 * TCR_Surface)") %>% 
  addReaction("TCR_ER", "", rate = "k_B * TCR_ER") %>%  
  
  
  # Export and transport from ER to Cell Surface
  addReaction("TCR_ER", "TCR_Golgi_1", rate = "k_E * TCR_ER / (1 + k_IM_3 * TCR_Surface)") %>% 
  addReaction("TCR_Golgi_1", "TCR_Golgi_2", rate = "k_E * TCR_Golgi_1") %>% 
  addReaction("TCR_Golgi_2", "TCR_Surface", rate = "k_E * TCR_Golgi_2") %>% 
  
  # Internalization to Endosom and Recycling 
  addReaction("TCR_Surface", "TCR_Endo", rate = "k_I * TCR_Surface") %>%
  addReaction("TCR_Endo", "TCR_Surface", rate = "k_PM_4 * TCR_Surface * TCR_Endo") %>% 
  
  # Degradation in Endosom
  addReaction("TCR_Endo", "", rate = "k_D * TCR_Endo") %>% 
  {.}

# mysteadieslate <- steadyStates(r) produces negative results!

datalist <- as.datalist(data)
cond.grid <- attr(datalist, "condition.grid")

mysteadies <- c(k_B = "k_P/TCR_ER")

# Define Observables
observables <- c(Total = "log2(TCR_ER + TCR_Endo + TCR_Golgi_1 + TCR_Golgi_2 + TCR_Surface + offset_Total) + scale_Total",
                 Surface = "log2(TCR_Surface + offset_Surface) + scale_Surface")

# Setting up SS constraint 
deqns <- as.eqnvec(r) %>% as.character()
names(deqns) <- paste0("d", r$states)


ssCnstrDatalist <- data.frame(time = rep(seq(6000, 20160, length.out = 20), length(deqns)) %>% sort())
ssCnstrDatalist$name <- rep(names(deqns), 5)
ssCnstrDatalist$value <- 0
ssCnstrDatalist$sigma <- 1e-5
ssCnstrDatalist$condition <- "OPEN"
ssCnstrDatalist <- as.datalist(ssCnstrDatalist)




innerpars <- c(getParameters(r), getSymbols(observables)) %>% unique()
trafo <- NULL %>% 
  define("x~x", x = innerpars) %>% 
  insert("inner~steadyEqn", inner = names(mysteadies), steadyEqn = mysteadies) %>%
  branch(table = cond.grid) %>%
  insert("scale_Surface~0") %>%
  insert("x~exp10(x)", x = .currentSymbols[!grepl("nhill|scale", .currentSymbols)]) %>%
  define("x~0", x = r$states, conditionMatch = "CTRL") %>% 
  define("x~0", x = "k_P", conditionMatch = "CTRL") %>% 
  define("x~0", x = setdiff(r$states, "TCR_ER"), conditionMatch = "OPEN") %>% 
  {.}

fixed <- c("TCR_Golgi_1", "TCR_Golgi_2", "TCR_Endo", "TCR_Surface")

if (.build) {
  
  if (!dir.exists(file.path(.modelFolder, "wd"))) {
    dir.create(file.path(.modelFolder, "wd"), recursive = TRUE)
  }
  setwd(file.path(.modelFolder, "wd"))
  
  # Construct ODE Predictor
  x <- odemodel(r, modelname = "x", compile = F, fixed = fixed) %>% Xs()
  
  # Construct Prediction Function and Error Model function
  g <- Y(c(observables, deqns), f = x, attach.input = T, compile = F, modelname = "g")
  
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

times <- 10^seq(0, log10(20160 + 1), length.out = 300) - 1

datalist <- data %>% dplyr::select(time, name, value, condition, sigma) %>% as.datalist()
outerpars <- getParameters(p)
pars <- structure(rep(-1, length(outerpars)), names = outerpars)
parlb <- structure(rep(-5, length(outerpars)), names = outerpars)
parub <- structure(rep(5, length(outerpars)), names = outerpars)
prd <- (g*x*p)(times, pars, condition = "CTRL") %>% as.data.frame()


# (x*p)(times, pars) %>% plot()
plotPred(prd, data, bands = NULL, custom.condition.names = custom_cond_names, custom.names = custom_names, errorbars = T)


obj_data <- normL2(datalist, g*x*p)
obj_prior <- constraintL2(pars, sigma = 20, attr.name = "prior")
obj_ssconstr <- normL2(ssCnstrDatalist, g*x*p, attr.name = "ssCnstr")

obj <- Reduce("+", list(obj_data, obj_prior, obj_ssconstr))
obj(pars)


source(file.path(.scriptFolder, "S03-DistributedComputing.R"))
outFitHelix <- SendFitsHelix(
  obj,
  pars,
  jobname = paste0("msfitting_TREG_",.modelname),
  argsDC = list(cores = 1, no_rep = 200, recover = F, link = T),
  iterlim = 5e3,
  sd = 4,
  parlower = parlb,
  parupper = parub,
  fterm = 1e-5,
  mterm = 1e-5
)
# #
# outFitHelix$check()
# results <- outFitHelix$get()
# save(file = file.path(.resultsFolder, "Fit_results.Rdata"), results)
load(file.path(.resultsFolder, "Fit_results.Rdata"))
class(results) <- "parlist"
outframe <- as.parframe(results)
outframe$value <- outframe$value + 2119
plotValues_alt(outframe[1:300], tol = 0.1, value<1e4) + theme_TCR(base_size = 18)
bestfit <- as.parvec(outframe)

prd <- (g*x*p)(times, bestfit) %>% as.data.frame() %>% 
  filter(condition == "OPEN")
custom_cond_names <- c(OPEN = "with biotin") # Custom condition names
custom_names <- c(Total = "[TCR] total", Surface = "[TCR] surface") # Custom names
plotPred(prd, subset(data,condition == "OPEN") , bands = NULL, custom.condition.names = custom_cond_names, 
                custom.names = custom_names, errorbars = T, show.legend = FALSE, base_size = 18)
trajectories <- (x*p)(times, bestfit) %>% as.data.frame()
plot_x <- plotPred(trajectories, custom.condition.names = custom_cond_names)
plot_x


# ggsave(
#   file.path(.plotFolder,"plotPred.pdf"),
#   plotPred(prd, subset(data,condition == "OPEN") , bands = NULL, 
#                   custom.condition.names = custom_cond_names, 
#                   custom.names = custom_names, errorbars = T) + 
#     theme_TCR(base_size = 11) + 
#     theme(legend.position = "none"),
#   width = 21, height = 10, units = "cm"
# )

obj <- Reduce("+", list(obj_data, obj_ssconstr))
obj(bestfit)

refit <- trust(obj, bestfit, rinit = 1, rmax = 10)
bestfit <- refit$argument %>% as.parvec()

# outHelix_profs <- SendProfilesHelix(obj, bestfit, names(bestfit), jobname= paste0("profiling_TCRReg", .modelname), argsDC=list(recover = F, link = T, walltime="00:30:00", cores = 1),
#                                     method = "integrate",
#                                     stepControl = list(stepsize = 1e-4, min = 1e-4, max = Inf, atol = 1e-2, rtol = 1e-2, limit = 5e2, stop = "data"),
#                                     algoControl = list(reoptimize = T),
#                                     optControl = list(rinit = 0.1, rmax = 1, iterlim = 1e3, fterm = 1e-5, mterm = 1e-5))
# 
# outHelix_profs$check()
# profiles <- outHelix_profs$get()
# save(file = file.path(.resultsFolder, "profile_results.Rdata"), profiles)
load(file.path(.resultsFolder, "profile_results.Rdata"))
plotProfile(profiles, mode == "data") + theme(legend.position = "none")
plotProfile(subset(profiles, whichPar %in% c("k_NFB_1", "k_NFB_3", "k_PFB_4")), mode == "data", ncol = 3) + theme(legend.position = "none")
plotProfilesAndPaths(profiles, c("k_NFB_1", "k_NFB_3", "k_PFB_4"), ncols = 3)
plotProfilesAndPaths(profiles, c("offset_Surface", "offset_Total", "scale_Surface", "scale_Total"), ncols = 4)

plotProfilesAndPaths(profiles, c("k_I", "k_D", "TCR_ER"), ncols = 3)
