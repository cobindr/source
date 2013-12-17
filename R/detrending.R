# Using a background to detrend the foreground curve to find candidate pairs in certain distances.
# pwm1 is used as the reference for the position.
# 
# Author: yuehien
#######################################################################
# 2013-01-09: Is this an internal function?
setMethod("detrending", signature(x="cobindr"),
function(x, pwm1, pwm2, bin_length=20, overlap=0, abs.distance=FALSE, background.on=TRUE) {
    #######################################
    ### Check parameter for plausibilty ###
    #######################################    
    warn <- NA
    if (bin_length < 1) warn <- paste('bin_length ',bin_length,' is not appropriate.')
    if(nrow(x@pairs) < 1) warn <- 'The provided object does not contain any TFBS pairs, please make sure to call prediction and pair finding prior to this method.'
    if (overlap < 0) warn <- paste('overlap must be positive or zero.')

    if(!is.na(warn)){
      warning(warn)
      return(list(FALSE))
    }

    # length of binding site
    pwm1_length = ncol(x@pfm[[pwm1]])
    pwm2_length = ncol(x@pfm[[pwm2]])
    
    # include possible overlap
    if (overlap > pwm1_length || overlap > pwm2_length) {
      warning(paste('overlap ', overlap, ' is longer than binding motif.'))
      return(NULL)
    }
    
    # calculate number of bins
    if (x@configuration@max_distance > 0)
        bins = as.integer(x@configuration@max_distance / bin_length)
    else bins = as.integer(200 / bin_length)

    if (bins < 2) {
        warning('Not enough bins (only ', bins + 1, 
		') for detrending analysis. Please adjust either the bin_length parameter here or the max_distance parameter in the config file!')
        return(NULL)
    }
    #######################################
    ### --------- count pairs --------- ###
    #######################################
	pair_positive = paste(c(pwm1, pwm2), collapse=' ')
	pair_negative = paste(c(pwm2, pwm1), collapse=' ')
	occurrences_p = bg_occurrences_p = x_p = rep(x@configuration@pseudocount, bins) # added pseudo-count (instead of rep(0, bins))
	occurrences_n = bg_occurrences_n = x_n = rep(x@configuration@pseudocount, bins) # here as well
	# filter for pair of interest
	pairs_positive = x@pairs[J(pair_positive), nomatch=0]
	pairs_negative = x@pairs[J(pair_negative), nomatch=0]
	
	if (background.on) {
        bg_pairs_positive = x@bg_pairs[J(pair_positive), nomatch=0]
	    bg_pairs_negative = x@bg_pairs[J(pair_negative), nomatch=0]
	}
	
	if (nrow(pairs_positive) == 0 && nrow(pairs_negative) == 0) {
	    warning(paste('No pairs found for PWMs ',pwm1,'/',pwm2,'.'))
        return(NULL)
	}

	lapply(1:bins, function(jj) {
	#for (jj in 1:bins) {
	    ###########################################################
        ### --------- count pairs where pwm1 is first --------- ###
        ###########################################################
        # foreground
        occurrences_p[jj] <<- occurrences_p[jj] + nrow(pairs_positive[distance >= (jj-1) * bin_length + 1 - overlap & distance < jj * bin_length + 1 - overlap & strand1 == 1, nomatch=0])
        
        if (pwm1 != pwm2) {
            occurrences_p[jj] <<- occurrences_p[jj] + nrow(pairs_negative[distance >= (jj-1) * bin_length + 1 - overlap & distance < jj * bin_length + 1 - overlap & strand2 == 2, nomatch=0])
        }

        x_p[jj] <<- jj * bin_length

        # background
        if (background.on) {
            bg_occurrences_p[jj] <<- bg_occurrences_p[jj] + nrow(bg_pairs_positive[distance >= (jj-1) * bin_length + 1 - overlap & distance < jj * bin_length + 1 - overlap & strand1 == 1, nomatch=0])
            if (pwm1 != pwm2) {
                bg_occurrences_p[jj] <<- bg_occurrences_p[jj] + nrow(bg_pairs_negative[distance >= (jj-1) * bin_length + 1 - overlap & distance < jj * bin_length + 1 - overlap & strand2 == 2, nomatch=0])
            }
        }
        
        ###########################################################
        ### --------- count pairs where pwm2 is first --------- ###
        ###########################################################
        # foreground
        occurrences_n[jj] <<- occurrences_n[jj] + nrow(pairs_positive[distance >= (jj-1) * bin_length + 1 - overlap & distance < jj * bin_length + 1 - overlap & strand1 == 2, nomatch=0])
        
        if (pwm1 != pwm2) {
            occurrences_n[jj] <<- occurrences_n[jj] + nrow(pairs_negative[distance >= (jj-1) * bin_length + 1 - overlap & distance < jj * bin_length + 1 - overlap & strand2 == 1, nomatch=0])
        }

        x_n[jj] <<- jj * bin_length

        # background
        if (background.on) {
            bg_occurrences_n[jj] <<- bg_occurrences_n[jj] + nrow(bg_pairs_positive[distance >= (jj-1) * bin_length + 1 - overlap & distance < jj * bin_length + 1 - overlap & strand1 == 2, nomatch=0])
            if (pwm1 != pwm2) {
               bg_occurrences_n[jj] <<- bg_occurrences_n[jj] + nrow(bg_pairs_negative[distance >= (jj-1) * bin_length + 1 - overlap & distance < jj * bin_length + 1 - overlap & strand2 == 1, nomatch=0])
            }
        }
    })
    
	# if specified to not distinguish between positive and negative distances
	if(abs.distance) {
		occurrences_n = occurrences_p = occurrences_n + occurrences_p
		bg_occurrences_n = bg_occurrences_p = bg_occurrences_n + bg_occurrences_p
	}
	
    ##############
    # Detrending #
    ##############
	if (pwm1 == pwm2 || abs.distance) {
	    x = x_p
	    occurrences = occurrences_p
	    
	    if (!background.on) {
    	    bg_occurrences_p = predict(lm(occurrences ~ x), data.frame(x = x))
    	}
    	
	    bg_occurrences = bg_occurrences_p
	} else {
	    x = c(rev(-x_n), x_p)
	    occurrences = c(rev(occurrences_n), occurrences_p)
	    
	    if (!background.on) {
	        bg_occurrences_p = predict(lm(occurrences_p ~ x_p), data.frame(x_p = x_p))
            bg_occurrences_n = predict(lm(occurrences_n ~ x_n), data.frame(x_n = x_n))
        }
        
	    bg_occurrences = c(rev(bg_occurrences_n), bg_occurrences_p)
    }
	
	# normalised background
	yy = bg_occurrences / mean(bg_occurrences) * mean(occurrences)
	# detrending
	y = occurrences - yy
	
    list(x, occurrences, bg_occurrences, yy, y)
}
)
