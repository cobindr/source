# TODO: IO functions for jaspar and transfac position frequency matrices (pfm) and positions weight matrices (pwm)  which can be used later in the analysis 
# 
# Author: rob <r.lehmann@biologie.hu-berlin.de>
# Author: stefan <kroeger@informatik.hu-berlin.de>
###############################################################################
##for conversion from PFM to PWM (MEME format)

# require("gmp") # for factorize() function


## function to read the pfm files from data folder 
## and convert it to MEME PWM format (http://meme.sdsc.edu/meme/examples/sample-dna-motif.meme-io)
# 2013-01-09: Is this an internal function?
setMethod("read.pfm", signature(x="configuration"),
function(x) {
    file_path = x@pfm_path
	
	if(any(grep('^motifdb$',file_path, ignore.case=TRUE))) 
		return(query.motifDb(x))
	
	# find all files in the given folder
	pfm.files = list.files(file_path, pattern ="*.*", full.names=F, ignore.case = TRUE)
	if(length(pfm.files) < 1) {
		warning('no files in pfm_path')
		return(list())
	}
	# create empty list of pwms
	pwms.list = list()
    count.pwms = 1
	pfm.read.errors = NULL
	
	cat('\n reading pfm files:',file_path,'...\n')
	for(i in 1:length(pfm.files)){
		tmp.path = file.path(file_path,pfm.files[i])
		# test which type of file this is
		ending = tail(unlist(strsplit(pfm.files[i], '\\.')), 1)
		switch(ending,
			# case one: it's a *.pfm jaspar file (which has only one matrix as entry)
			pfm = {
				tmp.pfm = apply(as.matrix(read.table(tmp.path)),c(1,2), as.numeric);
				pfm.id = sub(".pfm","", pfm.files[i])
				rownames(tmp.pfm) = c("A", "C", "G", "T")
				pwms.list = c(pwms.list,list(tmp.pfm))
				names(pwms.list)[length(pwms.list)] = pfm.id
			},
			# case two: it's a *.cm jaspar files (which also has the row names within)	
			cm  = {
				tmp.pfm = apply(as.matrix(read.table(tmp.path), row.names=1), c(1,2), as.numeric);
				pfm.id = sub(".cm","", pfm.files[i])
				rownames(tmp.pfm) = c("A", "C", "G", "T")
				pwms.list = c(pwms.list,list(tmp.pfm))
				names(pwms.list)[length(pwms.list)] = pfm.id
			},
			# case three: it's a *.sites jaspar file (not yet implemented)
			sites = {warning('The option of reading *.sites files is not yet implemented. The file is ignored ...')
			},
			# case: it's a *.tfpfm file that contains multiple matrices
			tfpfm ={
				tmp.pfm = read.transfac.pfm(tmp.path)
# 				pfm.id = sub(".tfpfm","", pfm.files[i])
				pwms.list = c(pwms.list,tmp.pfm)
			},
			# the default case - it's a transfac-file, which might have multiple matrices 
			{
				tmp.pfm = read.transfac.pfm(tmp.path)
				cat(pfm.files[i],' ignored - unknown PFM file ending.' )
				pfm.read.errors = c(pfm.read.errors, c(pfm.files[i],"\n"))
				next; #advance loop to next file
			}
		)
	}
	cat('ignored files:\n', pfm.read.errors,'\n')
	return(pwms.list)
}
)

# function queries motifDb package for pfms which are specified
# in pairs of the configuration file 
# call motif by unique name in MotifDb via names(MotifDb)
query.motifDb <- function(x) {
	if(!require('MotifDb')) {
		warning('package MotifDb is reqired. Aborting...')
		return(list())
	}
	mat.names <- unique(unlist(sapply(x@pairs, function(p) strsplit(p,' '))))
	cat('accessing MotifDb package for the requested ',length(mat.names),' PFMs...\n')
	mat.idx <- c()
	for(ma in mat.names) {
		mat.idx <- c(mat.idx, grep(ma, names(MotifDb)))
	}
	return(as.list(MotifDb[unique(mat.idx)]))
# 	return(as.list(MotifDb[unique(match(mat.names, names(MotifDb)))]))
}

# function reads transfac's pfm format files and returns corresponding matrix
# 2013-01-09: This function does not belong to any class.
#             Is this an internal function?             
read.transfac.pfm = function(path) {
	con = file(path) 
	open(con)
	# list of pfms read
	res = list()
	pfm = NA
	# iterate through file
	in.pfm = FALSE
	while (length(line <- readLines(con, n = 1, warn = FALSE)) > 0) {
          if(any(grep('^//$', line, ignore.case=T))){ # new pfm section starts
            if(any(!is.na(pfm))) # save previous pfm
              res[[id]] = pfm
            pfm = matrix(NA, nrow=4, ncol=0) # initialize new pfm
			rownames(pfm) = c('A','C','G','T')
          }
          
          if(any(grep('^P0', line, ignore.case=T))) {#pfm section starts
            in.pfm = TRUE
            next
          }
          if(any(grep('^ID', line, ignore.case=T))){ #pfm name
            id = unlist(strsplit(line, split=" "))[2]
            print(id)
          }
          if(in.pfm & any(grep('^XX$', line, ignore.case=T))) #pfm section ends
            in.pfm = FALSE
          if(in.pfm)
            pfm = cbind(pfm, as.numeric(unlist(strsplit(line, split=" "))[2:5]))
	}

        # save the last/only pfm (otherwise, they are not in the list)
        res[[id]] <- pfm
        
	close(con)
# 	if(ncol(pfm) > 0) # if last pfm not terminated with //, add it from the buffer to the result
#           res[[id]] = pfm
	res
}

  
##function normalized matrix: all columns must have same column sum;
##nessecary to be appliable to other functions like Biostrings::PWM()
# 2013-01-09: This function does not belong to any class.
#             Is this an internal function?
normalize.matrix.colSum = function(matrix){
  if(length(unique(colSums(matrix))) > 1){
    lcm = as.numeric(lcm.vector(colSums(matrix)))
    matrix = t((lcm/colSums(matrix))*t(matrix))
  }
  return(matrix)
}

# #w = log2 ( ( f + sqrt(N) * p ) / ( N + sqrt(N) ) / p )
# w - is a weight for the current nucleotide we are calculating
# f - is a number of occurences of the current nucleotide in the current column (e.g., "1" for A in column 1, "8" for C etc)
# N - total number of observations, the sum of all nucleotides occurences in a column (13 in this example)
# p - [prior] [background] frequency of the current nucleotide; this one usually defaults to 0.25 (i.e. one nucleotide out of four)
#          
# 2013-01-09: This function does not belong to any class.
#             Is this an internal function?          
pfm2pwm = function(pfm) {
	pfm <- normalize.matrix.colSum(pfm)
	pwm = pfm
	p = .25
	for(x in 1:ncol(pfm)){
		for(y in 1:4){
			f = pfm[y,x]
			n = colSums(pfm)[x]
			pwm[y,x] = log2 ( ( f + sqrt(n) * p ) / ( n + sqrt(n) ) / p )
		}
	}
	if(any(is.na(pwm)))
		stop('PFM-PWM conversion produced NAs !')
	else if(any(which(pwm == Inf)))
		stop('PFM-PWM conversion produced Infs!')
	else if(any(which(pwm == -Inf)))
		stop('PFM-PWM conversion produced -Infs!')
	else if(any(which(pwm == NaN)))
		stop('PFM-PWM conversion produced NaN!')
	return(pwm)
}
###############################################################################
# function converts for each input PWM the predicted TFBS hits into a PWM. function
# is intended to be used together with the sequence logo creation function 
# 'plot.tfbslogo' in plot-methods.R 
# as.pfm - if TRUE, the plain base counts of the predicted binding sites is returned. If FALSE, normalized pwm is returned
# 2013-01-09: Is this an internal function?

setMethod("predicted2pwm", signature(x="cobindr"),
function(x, as.pfm=FALSE) {
	predPwm = list() 
	bs = x@binding_sites
	pwms = setdiff(unique(bs$pwm),'0')

	for(pwm.id in pwms) {
		#not using grep due to reserved chars in pwm ids...
		# make reverse complement of the results from the '-' strand
		hits = c(as.character(bs$seq[which(bs$pwm == pwm.id & bs$strand == 1)]), unlist(lapply(as.character(bs$seq[which(bs$pwm == pwm.id & bs$strand == 2)]), function(seq) as.character(reverseComplement(DNAString(seq))))))

		if(length(hits)<1) {
			cat('cannot make logo: PWM',pwm.id,'has no hits in any sequence...\n')
		}

		pwm = consensusMatrix(DNAStringSet(hits), baseOnly=TRUE, as.prob=FALSE)

		pwm.s = matrix(NA,nrow=0, ncol=ncol(pwm))
		ps = c('A','C','G','T')

		for(p in ps)
			pwm.s = rbind(pwm.s, pwm[grep(paste('^',p,'$',sep=''),rownames(pwm),ignore.case=T),])
		rownames(pwm.s) = ps

		if(as.pfm) { predPwm[[pwm.id]] = pwm.s
		} else {
			for(i in 1:ncol(pwm.s))
				pwm.s[,i] = pwm.s[,i]/colSums(pwm.s)[i]
			predPwm[[pwm.id]] = pfm2pwm(pwm.s)
		}
	}
	return(predPwm)
}
)
          
# function returns the input PWM in the format required by the sequence logo 
# creation function 'plot.tfbslogo' in plot-methods.R
# as.pfm - if TRUE, the plain base counts of the predicted binding sites is returned. If FALSE, normalized pwm is returned
setMethod("input.pwm", signature(x="cobindr"),
function(x, as.pfm=FALSE) {
	#get list of input pwms used
	pwms = x@pfm
	if(length(pwms)<1) stop('object is missing the input id\'s')
	# iterate over input pwms and convert to required format
	inPwm = list() 
	for(pwm.id in names(pwms)) {
		pwm = pwms[[pwm.id]]
		if(as.pfm) {
			inPwm[[pwm.id]] <- pwm
		} else {
			inPwm[[pwm.id]] <- pfm2pwm(pwm / matrix(colSums(pwm),ncol=ncol(pwm),nrow=nrow(pwm),byrow=T))
		}
	}
	return(inPwm)
}
)

#function calculates least common multiple for a numeric vector
# 2013-01-09: This function does not belong to any class.
#             Is this an internal function?

lcm.vector <- function(numvec){
        tmp_lcm = numvec[1]
        for(i in 2:length(numvec)){
                tmp_lcm <- lcm.bigz(tmp_lcm,numvec[i])
        }
        return(tmp_lcm)
}

# checks all combinations of PFMs for high similarity (i.e. redundancy)
# thr - kullback-leibler distance threshold below which two matrices are considered similar
# typical values: 0.2 high, 0.3 moderate, 0.4 low stringency
# 2013-01-09: This function does not belong to any class.
#             Is this an internal function?
input.pfm.similarity <- function(pfms, thr=0.3) {
	combs <- combn(1:length(pfms),2)
	similar.matrices <- list()
	for(i in 1:ncol(combs)) {
		if(any(pfms[[combs[1,i]]]!=round(pfms[[combs[1,i]]]))) { 
			warning('matrix',names(pfms)[combs[1,i]],'is not a PFM (i.e. contains non-integers) and will be ignored...')
			next
		}
		if(any(pfms[[combs[2,i]]]!=round(pfms[[combs[2,i]]]))) { 
			warning('matrix',names(pfms)[combs[2,i]],'is not a PFM (i.e. contains non-integers) and will be ignored...')
			next
		}
		
		kl.mean.dist <- (kl.dist(pfms[[combs[1,i]]], pfms[[combs[2,i]]]) +
				kl.dist(pfms[[combs[2,i]]], pfms[[combs[1,i]]])) / 2
		if(kl.mean.dist < thr) {
			mn <- c(names(pfms)[combs[1,i]], names(pfms)[combs[2,i]])
			similar.matrices[[paste(mn, collapse='_')]] <-  
					list(mat1=mn[1], mat2=mn[2], mean.kl=kl.mean.dist)
		}
	}
	return(similar.matrices)
}

# computes kullback-leibler distance between two matrices
# mat1 - count matrix 1
# mat2 - count matrix 2
# 2013-01-09: This function does not belong to any class.
#             Is this an internal function?
kl.dist <- function(mat1, mat2) {
	# 	if(!require('gplots')) {
	# 		warning('package gplots required for this function, returning...')
	# 		return()
	# 	}
	# functions only works with position frequency matrices, so check if only ints 
	for(i in 1:length(mat1)) {
		if(mat1[i] != round(mat1[i])){
			warning('problem with mat1: function only applicable to PFMs (integer counts)...')
			return()
		}
	}
	for(i in 1:length(mat2)) {
		if(mat2[i] != round(mat2[i])){
			warning('problem with mat2: function only applicable to PFMs (integer counts)...')
			return()
	  }
	}
	# check if pseudocounts are necessary
	if(any(mat1 == 0)) mat1 <- mat1 + 1
	if(any(mat2 == 0)) mat2 <- mat2 + 1	
	# make probabilities
	mat1 <- mat1 / matrix(rep(colSums(mat1),each=nrow(mat1)), nrow=nrow(mat1), ncol=ncol(mat1), byrow=FALSE)
	mat2 <- mat2 / matrix(rep(colSums(mat2),each=nrow(mat2)), nrow=nrow(mat2), ncol=ncol(mat2), byrow=FALSE)
	
	bases <- rownames(mat1)
	l1 <- ncol(mat1) # matrix lengths
	l2 <- ncol(mat2)
	w <- min(l1, l2)
	dists <- NULL
	
	if(l1 == l2) {
		mat1s <- list(mat1)
		mat2s <- list(mat2)
	} else if(l1 > l2) { # second mat is shorter
		A <- seq(0,l1-l2,by=1)
		mat1s <- lapply(A, function(x) mat1[,x:(x+l2)])
		mat2s <- rep(list(mat2),length(A))
	} else {  # first mat is shorter
		A <- seq(0,l2-l1,by=1)
		mat1s <- rep(list(mat2),length(A))
		mat2s <- lapply(A, function(x) mat2[,x:(x+l1)])
	}
	for(i in 1:length(mat1s))
		dists <- c(dists, 1 / w * sum(sapply(1:w, function(j) {
							sum(sapply(bases, function(b){mat1s[[i]][b,j] * log2(mat1s[[i]][b,j]/mat2s[[i]][b,j])} ))
						}))
			)
	max(dists)
}

          
# returns kullback-leibler distance matrix between two sets of matrices
# intended usage: 
# compare input and predicted motifs - pfm.similarity(input.pwm(runObj, as.pfm=TRUE), predicted2pwm(runObj, as.pfm=TRUE))
# compare input motifs: 

#heatmap.2( pfm.similarity(input.pwm(runObj, as.pfm=TRUE), input.pwm(runObj, as.pfm=TRUE)) , dendrogram='row')
# 2013-01-09: This function does not belong to any class.
#             Is this an internal function?    
pfm.similarity <- function(set1, set2) {
	rm.idx = c()
	for(i in 1:length(set1))
		if(any(set1[[i]]!=round(set1[[i]]))) { 
			warning('matrix',names(set1)[i],'is not a PFM (i.e. contains non-integers) and will be ignored...')
			rm.idx = c(rm.idx, i)
		}
	if(length(rm.idx)>0) set1 <- set1[-rm.idx]
	rm.idx = c()
	for(i in 1:length(set2))
		if(any(set2[[i]]!=round(set2[[i]]))) { 
			warning('matrix',names(set2)[i],'is not a PFM (i.e. contains non-integers) and will be ignored...')
			rm.idx = c(rm.idx, i)
		}
	if(length(rm.idx)>0) set2 <- set2[-rm.idx]
	
	sims <- matrix(NA,nrow=length(set1),ncol=length(set2), dimnames=list(names(set1), names(set2)))
	for(i in 1:length(set1))
		for(j in 1:length(set2))
			sims[i,j] <- (kl.dist(set1[[i]], set2[[j]]) + kl.dist(set2[[j]], set1[[i]])) / 2
	return(sims)
}
