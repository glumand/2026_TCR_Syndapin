### Multistart Cluster Function
SendFitsHelix <- function(objfun, pars, jobname, argsDC=list(), ...) {
  
  opts <- modifyList(list(partition="cpu-single", cores=20, nodes=1, mem_per_core=2,
                          walltime="00:30:00", machine="helix", no_rep=10, link = FALSE,
                          recover=FALSE, compile=FALSE, resetSeeds=TRUE, returnAll=TRUE),
                     argsDC)
  message(if(opts$recover) "Recovering..." else "Submitting ...")
  
  mstrustargs <<- c(list(objfun=objfun, center=pars, studyname=jobname,
                         cores=opts$cores, fits=opts$cores), list(...))
  on.exit(rm("mstrustargs", envir=.GlobalEnv), add=TRUE)
  
  iface <- tryCatch(eval(substitute(
    dMod::distributed_computing(
      { do.call(dMod::mstrust, mstrustargs) },
      jobname=JOB, partition=PART, cores=CORES, nodes=NODES, mem_per_core=MEM,
      walltime=WALL, machine=MACH, no_rep=NOREP, recover=REC,
      compile=COMP, link = LINK, resetSeeds=RESET, returnAll=RET),
    list(JOB=jobname, PART=opts$partition, CORES=opts$cores, NODES=opts$nodes,
         MEM=opts$mem_per_core, WALL=opts$walltime, MACH=opts$machine,
         NOREP=opts$no_rep, REC=opts$recover, COMP=opts$compile, LINK=opts$link,
         RESET=opts$resetSeeds, RET=opts$returnAll))),
    error=function(e){ message("Error in distributive_computing: ", e$message); NULL })
  
  ## save original methods
  .check <- iface$check; .get <- iface$get; .purge <- iface$purge
  
  iface$check <- function() !is.null(iface) && .check()
  iface$get   <- function() Filter(function(x)!inherits(x,"Error"),
                                   unlist(.get(), recursive=FALSE))
  iface$purge <- function() invisible(try(.purge(), silent=TRUE))
  
  class(iface) <- "HelixInterface"
  iface
}


### Multistart Cluster Function
SendMiniModelFitsHelix <- function(data, ssCnstrData, prd, pars, MMconds, jobname, argsDC=list(), ...) {
  
  opts <- modifyList(list(partition="cpu-single", cores=32, nodes=1, mem_per_core=2,
                          walltime="00:30:00", machine="helix", link = FALSE,
                          recover=FALSE, compile=FALSE, resetSeeds=TRUE, returnAll=TRUE),
                     argsDC)
  message(if(opts$recover) "Recovering..." else "Submitting ...")
  
  mstrustargs <<- c(list(objfun=NULL, center=pars, studyname=jobname,
                         cores=opts$cores), list(...))
  on.exit(rm("mstrustargs", envir=.GlobalEnv), add=TRUE)
  
  iface <- tryCatch(eval(substitute(
    dMod::distributed_computing(
      { 
        dataL <- as.datalist(mutate(data, condition = ifelse(condition != "CTRL", var_1, condition)))
        ssCnstrData$condition <- var_1
        ssCnstrDataL <- as.datalist(ssCnstrData)
        obj = normL2(dataL, prd) + constraintL2(pars, sigma = 20) + normL2(ssCnstrDataL, prd, attr.name = "ssCnstr")
        mstrustargs$objfun <- obj
        do.call(dMod::mstrust, mstrustargs) 
      },
      jobname=JOB, partition=PART, cores=CORES, nodes=NODES, mem_per_core=MEM,
      walltime=WALL, machine=MACH, var_values = MMconds, no_rep=NULL, recover=REC,
      compile=COMP, link = LINK, resetSeeds=RESET, returnAll=RET),
    list(JOB=jobname, PART=opts$partition, CORES=opts$cores, NODES=opts$nodes,
         MEM=opts$mem_per_core, WALL=opts$walltime, MACH=opts$machine, 
         REC=opts$recover, COMP=opts$compile, LINK=opts$link,
         RESET=opts$resetSeeds, RET=opts$returnAll))),
    error=function(e){ message("Error in distributive_computing: ", e$message); NULL })
  
  ## save original methods
  .check <- iface$check; .get <- iface$get; .purge <- iface$purge
  
  iface$check <- function() !is.null(iface) && .check()
  iface$get <- function() Filter(function(x)!inherits(x,"Error"), .get())
  iface$purge <- function() invisible(try(.purge(), silent=TRUE))
  
  class(iface) <- "HelixInterface"
  iface
}


### Profile function CLuster
SendProfilesHelix <- function(objfun, pars, whichPar, jobname="profiling", argsDC=list(), ...) {
  
  optsDC <- modifyList(list(partition="cpu-single", cores=1, nodes=1, mem_per_core=2,
                            walltime="00:30:00", machine="helix", no_rep=NULL, split=TRUE, parsPerNode=1,
                            recover=FALSE, compile=FALSE, link = FALSE, resetSeeds=TRUE, returnAll=TRUE), argsDC)
  
  v_list <- if(optsDC$split)
    profile_pars_per_node(pars[whichPar], optsDC$parsPerNode, "split")
  else profile_pars_per_node(pars[whichPar], optsDC$parsPerNode)
  
  if(optsDC$recover) message("Recovering...") else message(paste0("Submitting ", as.character(length(whichPar)), " Profiles to ", optsDC$machine))
  
  profargs <<- c(list(obj=objfun, pars=pars, whichPar=whichPar,
                      cores=optsDC$cores), list(...))
  on.exit(rm("profargs", envir=.GlobalEnv), add=TRUE)
  
  iface <- tryCatch(eval(substitute(
    dMod::distributed_computing(
      { pargs<-profargs; pargs$whichPar<-pargs$whichPar[var_1:var_2];
      if(!is.null(var_3)) pargs$side<-var_3; do.call(dMod::profile, pargs) },
      jobname=JOB, partition=PART, cores=CORES, nodes=NODES, mem_per_core=MEM,
      walltime=WALL, machine=MACH, no_rep=NOREP, var_values=VLIST, recover=REC,
      compile=COMP, link=LINK, resetSeeds=RESET, returnAll=RET),
    list(JOB=jobname, PART=optsDC$partition, CORES=optsDC$cores, NODES=optsDC$nodes,
         MEM=optsDC$mem_per_core, WALL=optsDC$walltime, MACH=optsDC$machine,
         NOREP=optsDC$no_rep, VLIST=v_list, REC=optsDC$recover, COMP=optsDC$compile,
         LINK=optsDC$link, RESET=optsDC$resetSeeds, RET=optsDC$returnAll))),
    error=function(e){message("Error in distributive_computing: ",e$message);NULL})
  
  .c<-iface$check; .g<-iface$get; .p<-iface$purge
  iface$check<-function() !is.null(iface)&&.c()
  iface$get <- function() {
    if(!iface$check()) message("Profile computation not finished. Returning partial results...")
    cluster_result <<- .g()
    if(!length(cluster_result)) return(NULL)
    profiles <- NULL
    for(i in cluster_result)
      if(!inherits(i[[1]],"Error"))
        profiles <- rbind(profiles, i)
    profiles
  }
  iface$purge<-function() invisible(try(.p(),silent=TRUE))
  structure(iface, class="HelixInterface")
}


### L1 Cluster function
SendL1FitsHelix <- function(objfun, pars, mu=0*pars, log10.lambda=seq(-2,2,len=10),
                            one.sided=FALSE, jobname, argsDC=list(), ...) {
  
  interfaces <- vector("list", length(log10.lambda)); cleanup_vars <- character(0)
  optsDefault <- list(partition="cpu-single", cores=20, nodes=1, mem_per_core=2,
                      walltime="00:30:00", machine="helix", no_rep=10, link = TRUE,
                      recover=FALSE, compile=FALSE, resetSeeds=TRUE, returnAll=TRUE)
  argsDC <- modifyList(optsDefault, argsDC)
  
  if(!argsDC$recover) message("Submitting ...") else message("Recovering...")
  
  cleanup_vars <- c("mstrustargs","dotArgs")
  
  on.exit({ for(v in cleanup_vars) if(exists(v, envir=.GlobalEnv, inherits=FALSE)) rm(list=v, envir=.GlobalEnv) }, add=TRUE)
  
  lambda <- 10^log10.lambda
  for(i in seq_along(lambda)) {
    message("Helix L1 job ", i, "/", length(lambda), " (lambda = ", lambda[i], ")")
    mstrustargs <<- c(list(objfun=objfun, center=pars, studyname=jobname, optmethod="trustL1", lambda = lambda[i],
                           mu=mu, one.sided=one.sided, cores = argsDC$cores, fits = argsDC$cores),list(...))
    
    this_job <- paste(jobname,i,sep="_")
    
    dc_call <- substitute(
      dMod::distributed_computing({do.call(dMod::mstrust, mstrustargs)},
                                  jobname=JOB, partition=PART, cores=CORES, nodes=NODES, mem_per_core=MEM,
                                  walltime=WALL, machine=MACH, no_rep=NOREP, recover=REC, compile=COMP,
                                  link=LINK, resetSeeds=RESET, returnAll=RET),
      list(JOB=this_job, PART=argsDC$partition, CORES=argsDC$cores, NODES=argsDC$nodes,
           MEM=argsDC$mem_per_core, WALL=argsDC$walltime, MACH=argsDC$machine,
           NOREP=argsDC$no_rep, REC=argsDC$recover, COMP=argsDC$compile,
           LINK=argsDC$link, RESET=argsDC$resetSeeds, RET=argsDC$returnAll)
    )
    
    interfaces[[i]] <- tryCatch(eval(dc_call),
                                error=function(e){ message("Error in distributed computing: ", e$message); NULL })
  }
  
  HelixL1Interface <- list(interfaces=interfaces, lambda=lambda)
  
  HelixL1Interface$check <- function() all(vapply(interfaces, function(x)!is.null(x)&&x$check(), logical(1)))
  
  HelixL1Interface$get <- function() {
    lapply(interfaces, function(iface) {
      if(is.null(iface)) return(NULL)
      res <- iface$get()
      class(res) <- "parlist"
      reps <- unlist(res, recursive = FALSE)            # flatten rep-level
      Filter(function(x) !inherits(x,"Error"), reps)
    })
  }
  
  HelixL1Interface$purge <- function(){ for(i in seq_along(interfaces))
    if(!is.null(interfaces[[i]]))
      tryCatch(interfaces[[i]]$purge(),
               error=function(e) warning("Purge failed for L1 job ",i,": ",e$message))
    invisible(TRUE) }
  
  class(HelixL1Interface) <- "HelixL1Interface"; 
  return(HelixL1Interface)
}


SendFitsKnechte <- function(objfun, pars, machines, jobname,
                            argsDC=list(), ...) {
  
  opts <- modifyList(list(maxcores=24, recover=FALSE,
                          compile=FALSE, wait=FALSE, studyname="modelfit"),
                     argsDC)
  
  message(if(opts$recover) "Recovering..." else "Submitting ...")
  
  mstrustargs <<- list(
    objfun=objfun,
    center=pars,
    studyname=opts$studyname,
    cores=min(dMod::detectFreeCores(), opts$cores)
  )
  dotArgs <<- list(...)
  on.exit(rm(list=c("mstrustargs","dotArgs"), envir=.GlobalEnv), add=TRUE)
  
  iface <- tryCatch(
    dMod::runbg(
      { do.call(dMod::mstrust, c(mstrustargs, dotArgs)) },
      machine=machines, filename=jobname,
      recover=opts$recover, compile=opts$compile, wait=opts$wait
    ),
    error=function(e){ message("runbg error: ", e$message); NULL }
  )
  
  if (is.null(iface)) return(NULL)
  
  iface$filename <- jobname
  iface$machines <- machines
  class(iface) <- "KnechteInterface"
  iface
}