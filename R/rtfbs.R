# rtfbs.R
#
# Author: Yue-Hien, stefan <kroeger@informatik.hu-berlin.de>
#
# description: uses RTFBS (http://compgen.bscb.cornell.edu/rtfbs/) package to do a binding site search with FDR thresholding
###############################################################################

# Transcription factor binding site prediction method using RTFBS
# if append=T, predicted hits are appended to the hits in the input object

setMethod("rtfbs", signature(x="cobindr"),
function(x, append=F, background_scan = FALSE, n.cpu=NA) {

	warn = NA
	if(!require(rtfbs)) warn <-  'package rtfbsrequired for this function, returning...'
	if(length(x@pfm) < 1) warn <-  'No PFMs provided, returning...'
	if(!background_scan & length(x@sequences) < 1)  'No input sequences provided, returning...'
	if(background_scan & length(x@bg_sequences) < 1)  'No input background sequences provided, returning...'
	if(!is.na(warn)) {
		warning(warn)
		return(x)
	}
		
        # checking for fdrThresholding
	fdrThreshold = x@configuration@fdrThreshold > 0

	# convert pfm to rtfbs pwm format if necessary
	pwms = lapply(x@pfm, function(x) {
				# determine if integer or float values provided
				is.float = any(x != apply(x,c(1,2), as.integer))
				if(is.float) {
					return(x)
				} else {
					return(t(log2((x+1)/colSums(x+1)))) }
			})
	
	# it is only a normal check
	if(any(sapply(pwms,function(x) (Inf %in% x | -Inf %in% x | NA %in% x) )))
          stop('PFM - PWM conversion yielded Inf/-Inf/NA...')

	
	# check if background sequences should be scanned
	if (background_scan)
          target_sequences = x@bg_sequences
	else
          target_sequences = x@sequences
	
	# set number of cores to use, if not specified
	n.cpu <- determine.cores.option(n.cpu)
	
	# loop over all sequenceObjects
	a = Sys.time()
	cat('Starting binding site search with', length(pwms), 'PWMs.\n')
        func.name <- 'rtfbs.intern'
        func.opt <- list(pwms, fdrThreshold)
        results <- parallelize(func.name, target_sequences, func.opt, n.cpu)
        
        
        binding_sites = do.call(rbind, results)
	
	if(nrow(binding_sites) > 0) {
        # converting to dataframe
	    binding_sites.df = data.frame(uid = as.numeric(1:nrow(binding_sites)))
	    binding_sites.df = cbind(binding_sites.df, seqObj_uid = as.numeric(binding_sites[, 2]))
	    binding_sites.df = cbind(binding_sites.df, pwm = factor(binding_sites[, 3]))
	    binding_sites.df = cbind(binding_sites.df, start = as.numeric(binding_sites[, 4]))
	    binding_sites.df = cbind(binding_sites.df, end = as.numeric(binding_sites[, 5]))
	    binding_sites.df = cbind(binding_sites.df, score = as.numeric(binding_sites[, 6]))
	    binding_sites.df = cbind(binding_sites.df, seq = binding_sites[, 7], stringsAsFactors=FALSE)
	    binding_sites.df = cbind(binding_sites.df, strand = factor(binding_sites[, 8]))
	    binding_sites.df = cbind(binding_sites.df, source = factor(binding_sites[, 9]))
    	rm(binding_sites)

	    cat('Found', nrow(binding_sites.df), 'binding sites.\n')
	    
	    if (background_scan) {
	        if(append)
		        x@bg_binding_sites = rbind(x@bg_binding_sites, binding_sites.df)
	        else
		        x@bg_binding_sites = binding_sites.df
	    }
	    else {
	        if(append)
		        x@binding_sites = rbind(x@binding_sites, binding_sites.df)
	        else
		        x@binding_sites = binding_sites.df	
	    }
	}
	else {
	    cat('No binding sites have been found.\n')
	}
	
	print(Sys.time() - a)
	x
}
)


# ------------ #
# rtfbs intern #
# ------------ #
# 2013-01-09: This is an internal function.
setMethod("rtfbs.intern", signature(x="SeqObj"),
function(x, opts) {
	suppressWarnings(require(rtfbs))
	pwms <- opts[[1]]
	fdrThreshold <- opts[[2]]

	# using matrix to save results
	binding_sites_cache = matrix('0', nrow=0, ncol=9)
    
	# convert sequence to ms object for RTFBS
	seq = as.character(x@sequence)
	msSeq = ms(seq)

	# define a Markov model of order 3
    markovModel = build.mm(msSeq, 3)
	
	# ---------- FDR THRESHOLDING START ---------- #
	if (fdrThreshold) {
	    # create random sequence for FDR thresholding
        simSeq = simulate.ms(markovModel, lengths.ms(msSeq))
    }
    # ---------- FDR THRESHOLDING END ---------- #
    
    
	# searching for binding sites with each PWM using the Markov model
	l = 0
	for (i in 1:length(pwms)) {
		pwm = pwms[[i]]

        l = l + 1

		################################
        ### Search on the '+' strand ###
        ################################
        evaluate = TRUE
        
		score = score.ms(msSeq, pwm, markovModel, TRUE, 0, '+')

        # ---------- FDR THRESHOLDING START ---------- #
		if (fdrThreshold) {
		    simScore = score.ms(simSeq, pwm, markovModel, TRUE, 0, '+')
		    
		    if (nrow(simScore) > 0) {
			    fdrMap = calc.fdr(msSeq, score, simSeq, simScore)

		        if(nrow(fdrMap) < 1) {
			        evaluate = FALSE
		        }
		    }
		    else {
                evaluate = FALSE
		    }
		}
		# ---------- FDR THRESHOLDING END ---------- #

		if(nrow(score) > 0 && evaluate) {
		    bindingSites = output.sites(score, scoreThreshold=0)
		    
		    # ---------- FDR THRESHOLDING START ---------- #
    		if (fdrThreshold) {
	    	    bindingSites = output.sites(score, fdrScoreMap=fdrMap, fdrThreshold=fdrThreshold)
	    	}
		    # ---------- FDR THRESHOLDING END ---------- #
		    
		    bs = matrix(NA, nrow=nrow(bindingSites), ncol=9)
		    
	        # save binding sites to data.frame
	        if(nrow(bindingSites) > 0) {
	            for (i in 1:nrow(bindingSites)) {
	                # save binding site in matrix
        		    binding_sites_cache = rbind(binding_sites_cache, c(1, x@uid, names(pwms[l]), bindingSites[i, 'start'], bindingSites[i, 'end'], bindingSites[i, 'score'], substr(seq, bindingSites[i, 'start'], bindingSites[i, 'end']), '1', 'rtfbs'))
        		}
        	}
        } 
	    
	    ################################
        ### Search on the '-' strand ###
        ################################
		evaluate = TRUE
        
		score = score.ms(msSeq, pwm, markovModel, TRUE, 0, '-')

        # ---------- FDR THRESHOLDING START ---------- #
		if (fdrThreshold) {
		    simScore = score.ms(simSeq, pwm, markovModel, TRUE, 0, '-')
		    
		    if (nrow(simScore) > 0) {
			    fdrMap = calc.fdr(msSeq, score, simSeq, simScore)

		        if(nrow(fdrMap) < 1) {
			        evaluate = FALSE
		        }
		    }
		    else {
                evaluate = FALSE
		    }
		}
		# ---------- FDR THRESHOLDING END ---------- #

		if(nrow(score) > 0 && evaluate) {
		    bindingSites = output.sites(score, scoreThreshold=0)
		    
		    # ---------- FDR THRESHOLDING START ---------- #
    		if (fdrThreshold) {
	    	    bindingSites = output.sites(score, fdrScoreMap=fdrMap, fdrThreshold=fdrThreshold)
	    	}
	    	# ---------- FDR THRESHOLDING END ---------- #
		    
		    bs = matrix(NA, nrow=nrow(bindingSites), ncol=9)
		    
	        # save binding sites to data.frame
	        if(nrow(bindingSites) > 0) {
	            for (i in 1:nrow(bindingSites)) {
	                # save binding site in matrix
                    binding_sites_cache = rbind(binding_sites_cache, c(1, x@uid, names(pwms[l]), bindingSites[i, 'start'], bindingSites[i, 'end'], bindingSites[i, 'score'], substr(seq, bindingSites[i, 'start'], bindingSites[i, 'end']), '2', 'rtfbs'))
        		}
        	}
        } 
	}
	
	binding_sites_cache
}
)
