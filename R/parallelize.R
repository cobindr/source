# wrapper function for parallel support
# this function is inspired by flowClust
# internal function

parallelize <- function(func.name, object, func.opt, n.cpu=NA){
	# test whether any parallelization is possible
	# if not, use the serial version of cobindR  

  n.cpu <- determine.cores.option(n.cpu)
  if ( (!is.element("parallel",installed.packages()[,1])) && 
      (!is.element("snowfall",installed.packages()[,1])) || n.cpu==1){
    message("Using the serial version of cobindR - function ", func.name)
    res <- lapply(object, func.name, func.opt)
  }

	# testing if parallel support via multicore is available
  else if(require(parallel)){
    message("Using the parallel (multicore) version of cobindR - function ", func.name," with ", n.cpu," cores")
    res <- mclapply(object, get(func.name), func.opt, mc.preschedule=TRUE, mc.cores=n.cpu)
  }
  
	# testing if parallel support via snowfall is available
	else if(suppressWarnings(require(snowfall))) {
			# Number of clusters
		if(is.na(n.cpu)){
			n.cpu<-sfCpus()
		}
		message("Using the parallel (snowfall) version of cobindR function \"", func.name, "\" with ", n.cpu, " cpus or cores")
		sfInit(parallel=TRUE , cpus = n.cpu)
		res <- sfLapply(object, get(func.name), func.opt)
		sfStop()
		return(res)
	} 
	
	# if no parallelization package is available
	else { 
		message("No parallel support available - using the serial version of cobindR - function ", func.name)
		res <- lapply(object, func.name, func.opt)
	} 

# 	else if(suppressWarnings(require(snow))) {
# 		# Number of clusters
#     if(is.na(n.cpu)){
#       n.cpu<-sfCpus()
#     }
#     message("Using the parallel (snow) version of cobindR function \"", func.name, "\" with ", n.cpu, " cpus or cores")
# 	cl <- makeCluster(n.cpu, type = "SOCK")
# 	res <- parLapply(cl, object, get(func.name), func.opt)
# 	stopCluster(cl)
#     return(res)
#   }
  
}

##############################
# validate number of specified cores
# 2013-01-09: This function does not belong to any classes.
#             This is an internal function.
determine.cores.option <- function(n.cpu=NA) {
  core.avail = 1 
  if(Sys.info()['sysname']=="Windows") return(core.avail)
  if(is.numeric(try(parallel::detectCores(), silent=T))) core.avail = parallel::detectCores()
  if(is.na(n.cpu)) return(core.avail)
  else if(n.cpu > core.avail) {
    warning('requested number of cores (',n.cpu,') exceeds available cores (',core.avail,'), using max. possible...')
    return(core.avail)
  } else if(n.cpu < 1) {
    warning('invalid requested number of cores (',n.cpu,'), setting to single core usage...')
    return(core.avail)
  } else
    return(n.cpu)
}
