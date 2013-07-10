# findPairs.R
#
# author: stefan <kroeger@informatik.hu-berlin.de>
# author: Yue-Hien
# description: find pairs of binding sites for every seqObject in a given runObject
# get.pairs returns all pairs and its distances as list

###############################################################################

setMethod("find.pairs", signature(x="cobindr"),
function(x, background_scan, n.cpu = NA) {
	
	if(any(x@configuration@pairs == 'none')) { # abort if no pair finding specified in configuration
		cat('no pairs specified in configuration: \'none\'\n')
		return(x)
	}
	
	warn = NA
	if(!background_scan & nrow(x@binding_sites) < 1) warn <- 'provided input object does not contain any predicted binding sites - please make sure to run prediction first (and see if there are hits)...'
	if(background_scan & nrow(x@bg_binding_sites) < 1) warn <- 'provided input object does not contain any predicted binding sites- please make sure to run prediction first (and see if there are hits)...'
	if(!is.na(warn)) {
		warning(warn)
		return(x)
	}
        
	# set number of cores to use, if not specified
	n.cpu <- determine.cores.option(n.cpu)
        
	a = Sys.time()

    # using matrix to save results
    pairs = matrix('0', nrow=0, ncol=8)

	if (background_scan) {
	    target_sequences = x@bg_sequences
	    binding_sites = x@bg_binding_sites
	} else {
		target_sequences = x@sequences
		binding_sites = x@binding_sites
	}

	if(any(x@configuration@pairs == 'all')) {
		cat('Searching for all combinations of motifs can take long...\n')
		poi = apply(combn(sort(names(x@pfm)),2),2,paste, collapse=' ')
	} else { poi = x@configuration@pairs }
	
	## new <fast> version 
	cat('Searching for pairs...\n')
	func.name <- 'find.pairs_intern'
	func.opt <- list(target_sequences, binding_sites, x@configuration@max_distance)
	res <- parallelize(func.name, poi, func.opt, n.cpu)
	pairs = do.call(rbind, res)
	## old <slow> version
	# ####
	# lapply version of slow
	# 	res = lapply(poi,function(pair){
	# 		pwm1 = unlist(strsplit(pair, ' '))[1]
	# 		pwm2 = unlist(strsplit(pair, ' '))[2]
	# 		pair_of_interest_sorted = paste(sort(c(pwm1, pwm2)), collapse = ' ')
	# 		cat('Searching for pair', pair_of_interest_sorted, '.\n')
	# 		# get only binding sites from pair of interest
	# 		binding_sites_filtered = binding_sites[binding_sites$pwm == pwm1 | binding_sites$pwm == pwm2, ]
	# 		# loop through each sequence
	# 		func.name <- 'find.pairs_intern'
	# 		func.opt <- list(pair_of_interest_sorted, binding_sites_filtered, x)
	# 		results <- parallelize(func.name, target_sequences, func.opt, n.cpu)
	# 		return(do.call(rbind, results))
	# 	})
	# 	pairs = do.call(rbind, res)
	# ###
	# for loop version of slow
	# 	pairs = NULL
	# 	for (pair_of_interest in poi) {
	# 		pwm1 = unlist(strsplit(pair_of_interest, ' '))[1]
	# 		pwm2 = unlist(strsplit(pair_of_interest, ' '))[2]
	# 		
	# 		pair_of_interest_sorted = paste(sort(c(pwm1, pwm2)), collapse = ' ')
	# 		cat('Searching for pair', pair_of_interest_sorted, '.\n')
	# 		# get only binding sites from pair of interest
	# 		binding_sites_filtered = binding_sites[binding_sites$pwm == pwm1 | binding_sites$pwm == pwm2, ]
	# 		
	# 		# loop through each sequence
	# 		func.name <- 'find.pairs_intern'
	# 		func.opt <- list(pair_of_interest_sorted, binding_sites_filtered, x)
	# 		results <- parallelize(func.name, target_sequences, func.opt, n.cpu)
	# 		
	# 		pairs <- rbind(pairs, do.call(rbind, results))
	# 	}
	
	if(nrow(pairs) > 0) {
		pairs.df = data.frame(uid = as.numeric(1:nrow(pairs)))
		pairs.df = cbind(pairs.df, seqObj_uid = as.numeric(pairs[, 2]))
		pairs.df = cbind(pairs.df, pair = factor(pairs[, 3]))
		pairs.df = cbind(pairs.df, strand1 = factor(pairs[, 4]))
		pairs.df = cbind(pairs.df, strand2 = factor(pairs[, 5]))
		pairs.df = cbind(pairs.df, bs_uid1 = as.numeric(pairs[, 6]))
		pairs.df = cbind(pairs.df, bs_uid2 = as.numeric(pairs[, 7]))
		pairs.df = cbind(pairs.df, distance = as.numeric(pairs[, 8]))
		cat('Found', nrow(pairs), 'pairs.\n')
		rm(pairs)
		
		# check if background sequences have been scanned
		if (background_scan) x@bg_pairs = pairs.df
		else x@pairs = pairs.df
		
	} else cat('No pairs have been found.\n')

	print(Sys.time() - a)
    x
}
)

# ----------------- #
# Find Pairs intern #
# ----------------- #
# 2013-01-09: How does multiple binding work?
# setMethod("find.pairs_intern", signature(x="SeqObj"),
# function(x, opts) {
#   # creating variables from the list
# 	pair_of_interest_sorted <- opts[[1]]
# 	binding_sites_filtered <- opts[[2]]
# 	runObj <- opts[[3]]
# 
#     pairs_cache = matrix('0', nrow=0, ncol=8)
# 
#     seqObj_bs = binding_sites_filtered[binding_sites_filtered$seqObj_uid == x@uid, ]
# 
#     # check if enough binding sites are available
#     if (nrow(seqObj_bs) > 1) {
#     
#         # sort binding sites after start position, very important !!!
#         seqObj_bs = seqObj_bs[with(seqObj_bs, order(start, pwm)), ]
# 
#         # get parameters for finding pairs
#         max_distance = runObj@configuration@max_distance
# 
#         # find pairs
#         for (i in 1:(nrow(seqObj_bs) - 1)) {
#             for (k in (i+1):nrow(seqObj_bs)) {
#                 distance_between_bs = seqObj_bs[k, 'start'] - seqObj_bs[i, 'end']
#                 # check if pair is within range
#                 if (distance_between_bs >= max_distance && max_distance > 0) {
#                     break
#                 }
# 
#                 # build pair
#                 pair = paste(sort(c(as.character(seqObj_bs[i, 'pwm']), as.character(seqObj_bs[k, 'pwm']))), collapse = ' ')
# 
#                 # check if pair is of interest
#                 if (pair == pair_of_interest_sorted) {
#                     # save pair to matrix
#                     pairs_cache = rbind(pairs_cache, c(1, seqObj_bs[k, 'seqObj_uid'], paste(c(as.character(seqObj_bs[i, 'pwm']), as.character(seqObj_bs[k, 'pwm'])), collapse=' '), seqObj_bs[i, 'strand'], seqObj_bs[k, 'strand'], seqObj_bs[i, 'uid'], seqObj_bs[k, 'uid'], distance_between_bs))
#                 }
#             }
#         }
#     }
# 
#     pairs_cache
# }
# )

# ----------------- #
# Find Pairs intern 
# optimized for large ncpu and large # of pairs #
# ----------------- #
# 2013-03-17: How does multiple binding work?
# setMethod("find.pairs_intern", signature(x="character"),
find.pairs_intern <- function(x, opts) {
  # creating variables from the list
	target_sequences <- opts[[1]]
	binding_sites <- opts[[2]]
	max_distance <- opts[[3]]
	
	pairs_cache = matrix('0', nrow=0, ncol=8)

	pwm1 = unlist(strsplit(x, ' '))[1]
	pwm2 = unlist(strsplit(x, ' '))[2]
	pair_of_interest_sorted = paste(sort(c(pwm1, pwm2)), collapse = ' ')
	cat('Searching for pair', pair_of_interest_sorted, '.\n')
	# get only binding sites froms pair of interest

	bs_filtered1 = binding_sites[binding_sites$pwm == pwm1,]
	bs_filtered2 = binding_sites[binding_sites$pwm == pwm2,]
	
	if(nrow(bs_filtered1) == 0 || nrow(bs_filtered2) == 0 ) return(pairs_cache)
	
	inter_seqObj = intersect(bs_filtered1$seqObj_uid, bs_filtered2$seqObj_uid)
# 	bs_filtered = unique(rbind(bs_filtered1[match(inter_seqObj,bs_filtered1$seqObj_uid),],bs_filtered2[match(inter_seqObj,bs_filtered2$seqObj_uid),]))
	bs_filtered = unique(rbind(bs_filtered1[bs_filtered1$seqObj_uid %in% inter_seqObj, ], bs_filtered2[bs_filtered2$seqObj_uid %in% inter_seqObj, ]))

	# check if pwm's share bs in same sequences if not return empty set
	if (nrow(bs_filtered) == 0) return(pairs_cache)

	uid_vec = sapply(target_sequences, function(seq) slot(seq,"uid"))
	seqObj_bs = target_sequences[match(unique(bs_filtered$seqObj_uid), uid_vec)]

	found = lapply(seqObj_bs, function(seq) {
		bs = bs_filtered[bs_filtered$seqObj_uid == slot(seq, "uid"),]
		
		# sort binding sites after start position, very important !!!
# 		bs = bs[with(bs, order(start, pwm)), ]
#		bs[with(bs, order(get('start'), get('pwm'))), ]
		bs = bs[with(bs, order(get('start'), get('pwm'))), ]

		 # find pairs
		if(nrow(bs) < 2) return(pairs_cache)
        for (i in 1:(nrow(bs) - 1)) {
            for (k in (i+1):nrow(bs)) {
                distance_between_bs = bs[k, 'start'] - bs[i, 'end']
                # check if pair is within range
                if (distance_between_bs >= max_distance && max_distance > 0) {
                    break
                }

                # build pair
                pair = paste(sort(c(as.character(bs[i, 'pwm']), as.character(bs[k, 'pwm']))), collapse = ' ')

                # check if pair is of interest
                if (pair == pair_of_interest_sorted) {
                    # save pair to matrix
                    pairs_cache = rbind(pairs_cache, c(1, bs[k, 'seqObj_uid'], paste(c(as.character(bs[i, 'pwm']), as.character(bs[k, 'pwm'])), collapse=' '), bs[i, 'strand'], bs[k, 'strand'], bs[i, 'uid'], bs[k, 'uid'], distance_between_bs))
                }
            }
        }
		return(pairs_cache)
	})
	return(do.call(rbind, found))
}
#)


# ----------------- #
# 
# returns all found pairs as data.frame in bed format for later 
# ----------------- #
setMethod ("get.pairs", signature(x="cobindr"), 
function(x, background=FALSE) {
	if(exists(as.character(substitute(x@pairs))) && length(x@pairs)>0){
		seq_slot = "sequences"
		pair_slot = "pairs"
		bs_slot = "binding_sites"
		seq_type = "sequence_type"
		if (background) {
			seq_slot = "bg_sequences"
			pair_slot = "bg_pairs"
			bs_slot = "bg_binding_sites"
			seq_type = "bg_sequence_type"
		}
		
		# get locations from orig. sequences and extract position of TFBS:
		seqUIDs = slot(x, pair_slot)$seqObj_uid

		bed_locs = lapply(slot(x, seq_slot), function(i) c(i@uid,i@location))
		bed_locs = matrix(unlist(bed_locs),ncol=2, byrow=T)
		bed_locs = bed_locs[match(seqUIDs, bed_locs[,1]),]
# 		bed_locs = matrix(unlist(strsplit(bed_locs[match(seqUIDs, bed_locs[,1]),2],":")), ncol=3, byrow=T)
# 		bed_locs = data.frame(as.character(bed_locs[,1]),as.numeric(bed_locs[,2]),as.numeric(bed_locs[,3]))
		bed_locs = t(apply(bed_locs, 1,function(loc) {
			if(length(grep(".+:.+:.+",loc[2]))> 0) return(unlist(strsplit(loc[2],":")))
			else return(c(NA,NA,NA))
		}))
		# names
		bed_names = paste(slot(x, pair_slot)$seqObj_uid,gsub(' ' ,'-',slot(x, pair_slot)$pair),sep="-")
		# score = distance between pairs
		bed_distance = slot(x, pair_slot)$distance
		# strand - always set "+" 
		bed_strand = rep("+",dim(slot(x, pair_slot))[1])
		# thickstart
		bs1 = slot(x, bs_slot)[slot(x, pair_slot)$bs_uid1,] 
		bs2 = slot(x, bs_slot)[slot(x, pair_slot)$bs_uid2,] 
		bed_thickstart = sapply(1:dim(bs1)[1],function(i) min(bs1$start[i], bs1$end[i], bs2$start[i], bs2$end[i])+bed_locs[i,2])
		# thickend
		bed_thickend = sapply(1:dim(bs1)[1],function(i) max(bs1$start[i], bs1$end[i], bs2$start[i], bs2$end[i])+bed_locs[i,2])
		# blockcounts
		bed_bc = 2
		# blockSizes
		bed_bs =apply(data.frame(
			slot(x, bs_slot)[slot(x, pair_slot)$bs_uid1,]$end - slot(x, bs_slot)[slot(x, pair_slot)$bs_uid1,]$start,
# 				slot(x,pair_slot)$distance,
			slot(x, bs_slot)[slot(x, pair_slot)$bs_uid2,]$end - slot(x, bs_slot)[slot(x, pair_slot)$bs_uid2,]$start),1, function(i) paste(i, collapse=","))
		# blockStarts
		bed_bstarts = apply(data.frame(
# 				bed_locs[,2] + 	
			slot(x, bs_slot)[slot(x, pair_slot)$bs_uid1,]$start,
# 				slot(x,pair_slot)$distance,
# 				bed_locs[,2] + 
			slot(x, bs_slot)[slot(x, pair_slot)$bs_uid2,]$start ),1, function(i) paste(i, collapse=","))
		bed_rgb = rep("255,255,255", nrow(bed_locs))
		
		# compose final BEDfile styled table
		bed.output = data.frame(bed_locs, 
					bed_names, 
					bed_distance, 
					bed_strand, 
					bed_thickstart, 
					bed_thickend, 
					bed_rgb)
		colnames(bed.output) = c("chr","start", "end", "name", "distance", "strand","areaOfInterest_start", "areaOfInterest_end", "color")
		
		#if somethings runs out of defined ranges:
		bed.output[which(bed.output$end < bed.output$areaOfInterest_end),]$areaOfInterest_end = bed.output[which(bed.output$end < bed.output$areaOfInterest_end),]$end

		return(bed.output)
	}
})

# get.significant.pairs returns all binding sites of found significant pairs 
# tsv file. If locations for the sequences are given, binding site locations
# will be mapped to these locations.
setMethod ("get.significant.pairs", signature(x="cobindr"),
function(x, pwm1, pwm2, bin_length=20, z_value=3, overlap=0, abs.distance=FALSE) {
    results = detrending(x, pwm1, pwm2, bin_length, overlap, abs.distance)
    if(is.null(results)) return() #if detrending failed stop here
    
    x_coord = results[[1]]
    occurrences = results[[2]]
    bg_occurrences = results[[3]]
    yy = results[[4]]
    y = results[[5]]
    
    # check if candidate pairs were found
    #if (TRUE %in% (y > z_value*sd(y)) || TRUE %in% (y < -z_value*sd(y))) {
    if (TRUE %in% (y > z_value*sd(y))) {
	    pair_bs = matrix(nrow=0, ncol=9)
		colnames(pair_bs) = c('seqUID', paste('seq1(', pwm1, ')', sep=''), 'loc1', 'strand1', paste('seq2(', pwm2, ')', sep=''), 'loc2', 'strand2', 'distance', 'source')

		candidates = matrix(nrow=0, ncol=5)    
		colnames(candidates) = c('Pair', 'start', 'end', 'Z-value', 'source')

		for (i in 1:length(x_coord)) {
			if (y[i] > z_value*sd(y)) {
				# foreground
				more_results = get.detrending_intern(pwm1, pwm2, x@sequences, abs.distance, x@pairs, x_coord[i], bin_length, overlap, x@binding_sites, pair_bs, candidates, y[i]/sd(y), 'FOREGROUND')
				
				pair_bs = more_results[[1]]
				candidates = more_results[[2]]
			}
		}
		# return binding_sites
		if (nrow(candidates) > 0) {
			results = list("pair_bs","candidates")
			results[[1]] = pair_bs
			results[[2]] = candidates
			return(results)
		}
	} else {
        cat('No overrepresented distances were found for pair', pwm1, pwm2, '.\n')
		return(NULL)
    }
}
)

# 2013-01-08: What does this function do???
#             This should be an internal function.
#             Does it belong to a class?
get.detrending_intern = function(pwm1, pwm2, sequences, abs.distance, pairs, x, bin_length, overlap, binding_sites, pair_bs, candidates, p_z_value, source) {
    pair_positive = paste(c(pwm1, pwm2), collapse=' ')
    pair_negative = paste(c(pwm2, pwm1), collapse=' ')
    
    # convert seqObj to list
    seqUID = list()
    for (seq in sequences) seqUID[as.character(seq@uid)] = seq@location
    
    if (abs.distance) {
        pairs_cache = pairs[pairs$pair == pair_positive & pairs$distance >= x - bin_length + 1 - overlap & pairs$distance < x + 1 - overlap, ]
        
        pairs_cache = rbind(pairs_cache, pairs[pairs$pair == pair_negative & pairs$distance >= x - bin_length + 1 - overlap & pairs$distance < x + 1 - overlap, ])
    }
    else {
        if (x > 0) {
            pairs_cache = pairs[pairs$pair == pair_positive & pairs$distance >= x - bin_length + 1 - overlap & pairs$distance < x + 1 - overlap & pairs$strand1 == 1, ]
            
            pairs_cache = rbind(pairs_cache, pairs[pairs$pair == pair_negative & pairs$distance >= x - bin_length + 1 - overlap & pairs$distance < x + 1 - overlap & pairs$strand2 == 2, ])
        }
        else {
            pairs_cache = pairs[pairs$pair == pair_positive & pairs$distance >= -x - bin_length + 1 - overlap & pairs$distance < -x + 1 - overlap & pairs$strand1 == 2, ]
            
            pairs_cache = rbind(pairs_cache, pairs[pairs$pair == pair_negative & pairs$distance >= -x - bin_length + 1 - overlap & pairs$distance < -x + 1 - overlap & pairs$strand2 == 1, ])
        }
    }
    
    for (tt in 1:nrow(pairs_cache)) {
        bs1 = binding_sites[binding_sites$uid == pairs_cache[tt, 'bs_uid1'], ]
        bs2 = binding_sites[binding_sites$uid == pairs_cache[tt, 'bs_uid2'], ]

        if (seqUID[as.character(pairs_cache[tt, 'seqObj_uid'])] != 'unknown') {
            location = unlist(strsplit(unlist(seqUID[as.character(pairs_cache[tt, 'seqObj_uid'])]), ':'))

            loc1 = paste(location[[1]], as.numeric(location[[2]]) + bs1$start - 1, as.numeric(location[[2]]) + bs1$end - 1, sep=':')
            loc2 = paste(location[[1]], as.numeric(location[[2]]) + bs2$start - 1, as.numeric(location[[2]]) + bs2$end - 1, sep=':')
        }
        else {
            loc1 = 'unknown'
            loc2 = 'unknown'
        }

        if (x > 0) {
            pm = 1
        }
        else {
            pm = -1
        }

        if (bs1$pwm == pwm1) {
            pair_bs = rbind(pair_bs, c(pairs_cache[tt, 'seqObj_uid'], bs1$seq, loc1, bs1$strand, bs2$seq, loc2, bs2$strand, pm * pairs_cache[tt, 'distance'], source))
        }
        else {
            pair_bs = rbind(pair_bs, c(pairs_cache[tt, 'seqObj_uid'], bs2$seq, loc2, bs2$strand, bs1$seq, loc1, bs1$strand, pm * pairs_cache[tt, 'distance'], source))
        }
    }
    
    if (x > 0) {
        p_start = x - bin_length + 1
    }
    else {
        p_start = x + bin_length - 1    
    }
    
    p_end = x
    
    cat('Found candidate pair in', source, pwm1, pwm2, 'in distance', p_start, '-', p_end, 'bp. Z-value is', p_z_value, '\n')
    
    candidates = rbind(candidates, c(pair_positive, p_start, p_end, p_z_value, source))
    
    list(pair_bs, candidates)
}