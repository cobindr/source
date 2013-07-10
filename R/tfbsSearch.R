# Author: rob <r.lehmann@biologie.hu-berlin.de>
###############################################################################

# function perpforms TFBS prediction denovo or based on transfac / jaspar matrices pwms
# if append=T, predicted hits are appended to the hits in the input object
setMethod("search.gadem", signature(x="cobindr"),
function(x, deNovo=FALSE, append=F, background_scan = FALSE) {
	# check input validity
	warn = NA
	if(!require(rGADEM)) warn <- 'package rGADEM required for this function, returning...'
	if(!deNovo & length(x@pfm) < 1) warn <- 'No PWMs provided, returning...'
	if(!background_scan & length(x@sequences) < 1) warn <- 'No input sequences provided, returning...'
	if(background_scan & length(x@bg_sequences) < 1) warn <- 'No input background sequences provided, returning...'
	if(!is.na(warn)) {
		warning(warn)
		return(x)
	}
	
	# collect pwms
	if(!deNovo) pwms = lapply(x@pfm, function(x) apply(x, c(1,2), as.numeric))
	else pwms = list()

	#convert seqs to DNAString
	if (!background_scan) {
		res.slot = 'binding_sites'
		seq.slot = 'sequences'
	} else {
		res.slot = 'bg_binding_sites'
		seq.slot = 'bg_sequences'
	}
	#convert sequence objects to DNAStringSet
	sequences = DNAStringSet(sapply(slot(x,seq.slot), function(x) as.character(x@sequence)))
	
	# TODO: provide background for analysis
	#perform prediction
	pVal = x@configuration@pValue

	if( (!deNovo) & length(pwms) > 0 )
		gadem <- lapply(pwms, function(x) GADEM(Sequences = sequences,verbose=1,Spwm=list(x),
							pValue=pVal))
	else
		gadem <- GADEM(sequences,verbose=1,pValue=pVal)
	
	binding.sites = matrix(NA,nrow=0,ncol=8)
	for(pwm.i in 1:nMotifs(gadem)) {
		pwm.id = names(gadem)[pwm.i]
		pwm.g = gadem[pwm.i]
		if(! 'motifList' %in% slotNames(pwm.g))
			next
		for(i in 1:length(pwm.g@motifList)) {
			if(length(pwm.g@motifList)<1) next
			seqObj_uid = sapply(pwm.g@motifList[[i]]@alignList, slot,'seqID')
			pwm = replicate(length(pwm.g@motifList[[i]]@alignList), paste(pwm.id,'(',pwm.g@motifList[[i]]@consensus,')',sep='')) 
			#TODO: find why returned start/end=0, use pos instead 
			start = sapply(pwm.g@motifList[[i]]@alignList, slot,'pos')
			end = start + nchar(pwm.g@motifList[[i]]@consensus)
			score = sapply(pwm.g@motifList[[i]]@alignList, slot,'pval')
			#convert strand mark to internal standard
			strand = sapply(pwm.g@motifList[[i]]@alignList, slot,'strand')
			strand[which(strand=='+')] = 1
			strand[which(strand=='-')] = 2
			#convert all sequence hits into leading strand for compatibility
			seq = sapply(pwm.g@motifList[[i]]@alignList, slot,'seq')
			for(j in which(strand=='-')) seq[j] = reverseComplement(DNAString(seq[j]))
			source = replicate(length(pwm.g@motifList[[i]]@alignList), 'rGADEM')
			bs = data.frame(seqObj_uid=seqObj_uid,pwm=pwm,start=start,end=end,score=score,seq=seq,strand=strand,source=source)
			binding.sites = rbind(binding.sites, bs)
		}
	}
	# add hits to object as appropriate and return
	if(nrow(binding.sites) > 0) {
		cat('found',nrow(binding.sites), 'hits in total.\n')
		if(append)
			slot(x, res.slot) = as.data.frame(rbind(slot(x, res.slot), cbind(uid=c(1:nrow(binding.sites)), binding.sites)))
		else
			slot(x, res.slot) = as.data.frame(cbind(uid=c(1:nrow(binding.sites)), binding.sites))
	} else
		warning('no hits were found!')
	x
}
)

# applies package Biostrings provided method matchPWM for TFBS prediction
setMethod("search.pwm", signature(x="cobindr"),
function(x, min.score = '80%', append=FALSE, background_scan = FALSE, n.cpu=NA) {
	
	warn = NA

	if(length(x@pfm) < 1) warn <- 'No PFMs provided, returning...'
	if(!background_scan & length(x@sequences) < 1) warn <- 'No input sequences provided, returning...'
	if(background_scan & length(x@bg_sequences) < 1) warn <- 'No input background sequences provided, returning...'
	if(!is.na(warn)) {
		warning(warn)
		return(x)
	}
	
	# set number of cores to use, if not specified
	n.cpu <- determine.cores.option(n.cpu)
	
	# convert PFMs into PWMs
	pwms = lapply(1:length(x@pfm), function(y) {
				pwm = x@pfm[[y]]
# 				# determine if integer or float values provided
				is.float = any(pwm != apply(pwm,c(1,2), as.integer))
				if(is.float) {
					return(pwm)
				} else {
				res = try(pfm2pwm(pwm))
				if(class(res)=='try-error')
					warning('PFM-PWM conversion failed for matrix ',names(x@pfm)[x],'!')
				return(res) }
			})

	names(pwms) = names(x@pfm)
	# sanity check
	if(any(lapply(pwms, class) == 'try-error'))
		error('it appears there were errors in the PFM-PWM conversion, returning...')
	cat('finding binding sites for',length(pwms),'PWMs...\n')
	
	#convert seqs to DNAString
	if (!background_scan) {
		res.slot = 'binding_sites'
		seq.slot = 'sequences'
	} else {
		res.slot = 'bg_binding_sites'
		seq.slot = 'bg_sequences'
	}

	#convert sequence objects to DNAStringSet
# 	sequences = DNAStringSet(sapply(slot(x,seq.slot), function(x) as.character(x@sequence)))
	
	func.name <- 'make.character'
	func.opts <- list()
	res <- parallelize(func.name, slot(x, seq.slot), func.opts, n.cpu)
	if(is.null(res)) { 
		warning('Could not extract sequences.')
		return(x)
	}
	sequences <- DNAStringSet(unlist(res))

	sequences = IRanges::as.list(append(sequences,reverseComplement(sequences))) # add reverseComp for search
	str = c(rep(1,length(slot(x,seq.slot))), rep(2,length(slot(x,seq.slot)))) # strand info
	uid = as.numeric(rep(sapply(slot(x,seq.slot), attr, 'uid'),2)) # ids
	
	#iterate over PWMs
        func.name <- 'find.hits'
        func.opts <- list(pwms, sequences, uid, min.score, str)
        bs.hits <- parallelize(func.name, 1:length(pwms), func.opts, n.cpu)

	#collect results from different PWMs
	binding.sites = matrix(NA,nrow=0,ncol=8)
	for (i in 1:length(bs.hits)) 
		if(nrow(bs.hits[[i]]) > 0) 
			binding.sites = rbind(binding.sites, bs.hits[[i]])
	
	if(nrow(binding.sites) > 0) {
		cat('found',nrow(binding.sites), 'hits in total.\n')
		binding.sites$uid = 1:nrow(binding.sites)
		if(append)
			slot(x, res.slot) = rbind(slot(x, res.slot), binding.sites) 
		else
			slot(x, res.slot) = binding.sites
	} else
		warning('no hits were found!')
	x
        
      }
)

# 2013-03-11 internal function
make.character <- function(seq, opts){
	return(as.character(seq@sequence))
}

# 2013-03-11 internal function
find.hits <- function(count, opts) {
# 	suppressWarnings(require(IRanges))
	pwms <- opts[[1]]
	seq <- opts[[2]]
	uid <- opts[[3]]
	min.score <- opts[[4]]
	str <- opts[[5]]

	cat('finding hits for PWM',names(pwms)[count],'...\n')
	# matching
	hits = lapply(seq, matchPWM, pwm=pwms[[count]], min.score=min.score)  
	## Post-calculate the scores of the hits:
											
  scores = lapply(hits, function(h) PWMscoreStartingAt(pwm=pwms[[count]], subject=subject(h), starting.at=start(h)) )
		
	#iterate over sequence results
	res = matrix(NA,nrow=0,ncol=8)
	for(h.i in 1:length(hits)) {
		h = hits[[h.i]]
		l = length(h)
		if(l < 1) next
		if(str[h.i]==1 | str[h.i]==0) {
		start.coord = start(h)
		end.coord = end(h)
		new.seq = as.character(h)
		} else { # convert hit sequences and reverseComp
		new.seq = as.character(reverseComplement(h))
		# need to subtract the sequence length and switch start and end coord
		start.coord = length(seq[[h.i]]) - end(h) + 1
		end.coord = length(seq[[h.i]]) - start(h) + 1
		}
		bs = data.frame(seqObj_uid = rep(uid[h.i], l), pwm = rep(names(pwms[count]), l), start = start.coord, end = end.coord, score = scores[[h.i]], seq = new.seq, strand = rep(str[h.i], l), source = rep('matchPWM', l), stringsAsFactors=FALSE)

		res = rbind(res, bs)
	}
	cat('found ',nrow(res),'hits for PWM',names(pwms)[count],'...\n')
  return(res)
}


