# Using a background to detrend the foreground curve to find candidate pairs in certain distances.
# pwm1 is used as the reference for the position.
# 
# Author: yuehien
#######################################################################
# 2013-01-09: Is this an internal function?
setMethod("detrending", signature(x="cobindr"),
function(x, pwm1, pwm2, bin_length=20, overlap=0, abs.distance=FALSE) {
    #######################################
    ### Check parameter for plausibilty ###
    #######################################    
    warn <- NA
    if (bin_length < 1) warn <- paste('bin_length ',bin_length,' is not appropriate.')
    if(nrow(x@pairs) < 1) warn <- 'The provided object does not contain any TFBS pairs, please make sure to call prediction and pair finding prior to this method.'
    if (overlap < 0) warn <- paste('overlap must be positive or zero.')
    if (length(x@binding_sites) > 0) pwms = setdiff(unique(x@binding_sites$pwm),'0')
    if(!is.na(warn)){
      warning(warn)
      return(list(FALSE))
    }

    for(pwm in c(pwm1, pwm2)) {
        if(!pwm %in% pwms) {
    			warning('specified PWM ',pwm,'is not valid! returning...\n')
          return(NULL)
        }
    }
    if(length(which(paste(pwm1,pwm2 ) == x@pairs | paste(pwm2,pwm1 ) == x@pairs)) < 1) {
	    warning(paste('No pairs found for PWMs ',pwm1,'/',pwm2,', returning...'))
      return(NULL)
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

	for (jj in 1:bins) {
	    ###########################################################
        ### --------- count pairs where pwm1 is first --------- ###
        ###########################################################
        # foreground
        occurrences_p[jj] = occurrences_p[jj] + nrow(x@pairs[x@pairs$pair == pair_positive & x@pairs$distance >= (jj-1) * bin_length + 1 - overlap & x@pairs$distance < jj * bin_length + 1 - overlap & x@pairs$strand1 == 1, ])
        occurrences_p[jj] = occurrences_p[jj] + nrow(x@pairs[x@pairs$pair == pair_negative & x@pairs$distance >= (jj-1) * bin_length + 1 - overlap & x@pairs$distance < jj * bin_length + 1 - overlap & x@pairs$strand2 == 2, ])

        x_p[jj] = jj * bin_length

        # background
        bg_occurrences_p[jj] = bg_occurrences_p[jj] + nrow(x@bg_pairs[x@bg_pairs$pair == pair_positive & x@bg_pairs$distance >= (jj-1) * bin_length + 1 - overlap & x@bg_pairs$distance < jj * bin_length + 1 - overlap & x@bg_pairs$strand1 == 1, ])
        bg_occurrences_p[jj] = bg_occurrences_p[jj] + nrow(x@bg_pairs[x@bg_pairs$pair == pair_negative & x@bg_pairs$distance >= (jj-1) * bin_length + 1 - overlap & x@bg_pairs$distance < jj * bin_length + 1 - overlap & x@bg_pairs$strand2 == 2, ])
        
        ###########################################################
        ### --------- count pairs where pwm2 is first --------- ###
        ###########################################################
        # foreground
        occurrences_n[jj] = occurrences_n[jj] + nrow(x@pairs[x@pairs$pair == pair_positive & x@pairs$distance >= (jj-1) * bin_length + 1 - overlap & x@pairs$distance < jj * bin_length + 1 - overlap & x@pairs$strand1 == 2, ])
        occurrences_n[jj] = occurrences_n[jj] + nrow(x@pairs[x@pairs$pair == pair_negative & x@pairs$distance >= (jj-1) * bin_length + 1 - overlap & x@pairs$distance < jj * bin_length + 1 - overlap & x@pairs$strand2 == 1, ])

        x_n[jj] = jj * bin_length

        # background
        bg_occurrences_n[jj] = bg_occurrences_n[jj] + nrow(x@bg_pairs[x@bg_pairs$pair == pair_positive & x@bg_pairs$distance >= (jj-1) * bin_length + 1 - overlap & x@bg_pairs$distance < jj * bin_length + 1 - overlap & x@bg_pairs$strand1 == 2, ])
        bg_occurrences_n[jj] = bg_occurrences_n[jj] + nrow(x@bg_pairs[x@bg_pairs$pair == pair_negative & x@bg_pairs$distance >= (jj-1) * bin_length + 1 - overlap & x@bg_pairs$distance < jj * bin_length + 1 - overlap & x@bg_pairs$strand2 == 1, ])
    }
    
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
	    bg_occurrences = bg_occurrences_p
	} else {
	    x = c(rev(-x_n), x_p)
	    occurrences = c(rev(occurrences_n), occurrences_p)
	    bg_occurrences = c(rev(bg_occurrences_n), bg_occurrences_p)
    }
	
	# normalised background
	yy = bg_occurrences / mean(bg_occurrences) * mean(occurrences)
	# detrending
	y = occurrences - yy
	
    list(x, occurrences, bg_occurrences, yy, y)
}
)
