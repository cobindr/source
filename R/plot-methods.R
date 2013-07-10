# methods to visualize TFBS prediction result and pair analysis
# 
# Author: rob <r.lehmann@biologie.hu-berlin.de>
# 		  stefan <kroeger@informatik.hu-berlin.de>
###############################################################################

###############################################################################
# plot.tfbslogo
#
# produces sequence Logo of specified pfms
# if a file path is specified in pdf.name, sequences logos will be written into the specified file
# currently uses seqLogo::makePWM() methd to generate pwm's from x@pfm
setMethod("plot.tfbslogo", signature(x="cobindr"),
function(x, pwms) { 

	if(!require('seqLogo')) {
		warning('package seqLogo required for this function, returning...')
    return()
  }
	
	if(is.null(x@pfm) || is.null(pwms) || length(x@pfm) < 1 || length(pwms) < 1) {
		warning('bad parameters  - check input')
		return()
	}
	pwm.no = match(pwms, names(x@pfm))
	sapply(pwm.no, function(pwm.id) {
		cat('plotting pwm',pwm.id,'\n')

		#normalize pwm to colSum 1 --> requirement of seqLogo
		tmp.pwm <- makePWM(normalize.matrix.colSum(x@pfm[[pwm.id]]) / colSums(normalize.matrix.colSum(x@pfm[[pwm.id]])))
		
		#make Logo plot and add pwm name to it
		seqLogo(tmp.pwm, ic.scale=TRUE, xaxis=TRUE, yaxis=TRUE, xfontsize=15, yfontsize=15)
		grid.text(pwm.id,.5,.9,gp=gpar(fontsize=20))
	})
}
)

###############################################################################
# plot.tfbs.venndiagram
#
# the distribution of PWM hits over the sequences is visualized as Venn diagram
# if alist of PWM names is provided, only these PWMs are included in the Venn diagram
# if include.empty.seqs == TRUE, sequences without hits of the specified PWMs are included in the diagram  
# if a file path is specified in pdf.name, the diagram will be written into the specified file

setMethod("plot.tfbs.venndiagram", signature(x="cobindr"),
function(x, pwms, include.empty.seqs=FALSE) {
	
	if(!require(VennDiagram)) {
    warning('package VennDiagram required for this function, returning...')
    return()
}
   

	bs = x@binding_sites
	if(nrow(bs) < 1) {
	  warning('provided object contains no predicted TFBS hits, please ensure calling prediction prior to this method...')
          return()
	}
    
  sequence.ids = sapply(x@sequences,attr,'uid')
  
	if(missing(pwms)) {
		pwms = unique(bs$pwm)
		# remove empty entries from list of pwms (0)
		pwms = setdiff(pwms,'0')
		cat('PWM id list not specified, including all PWMs with hits (',paste(pwms,sep=','),')...\n')
	} else {
		# check if list of specified pwms contains only valid pwm ids
		pwm.names = names(x@pfm)
		if( length(setdiff(pwms,pwm.names)) > 0 ) {
			warning('specified list of PWMs contains unknown PWM identifiers: ',paste(setdiff(pwms,pwm.names),sep=','),'! returning...')
		return()
		}
		cat('specified PWM id list! Including ',paste(pwms,sep=','),'...\n')
	}

	if (length(pwms) > 4) {
    warning('venn diagrams cannot be shown for more than 4 pwms. use plot.tfbs.heatmap() or use less pwms, returning ...')
    return()
	}
	
  #holds sequence ids used for creating the venn diagram
	hits.lst = vector(mode='list')
	
  #get list of hits for each pwm
	for(pwm.id in pwms)
		hits.lst[[pwm.id]] = unique(bs$seqObj_uid[which(bs$pwm == pwm.id)])
  
  #include sequences without hits as 'None'
	if(include.empty.seqs)
		hits.lst[['None']] =  setdiff(sequence.ids, unique(bs$seqObj_uid))

	if(length(hits.lst) == 0) {
		warning('no hits for binding sites')
    return()
	}
	venn.plot = venn.diagram(hits.lst, main=paste(pwms, collapse=" : ") ,filename=NULL, scaled=TRUE)
	grid.draw(venn.plot)
	return(venn.plot)
}
)

##############################################################################
# plot.tfbs.heatmap
#
# plots a heatmap of overlaps between all specified PWMs
# for each overlap, the significance is determined via Fisher's exact test
# if a file path is specified in pdf.name, the diagram will be written into the specified file

setMethod("plot.tfbs.heatmap", signature(x="cobindr"),
function(x, pwms, include.empty.seqs=FALSE) {
  if(!require(gplots)) {
    warning('package gplots required for this function, returning...')
    return()
  }
  
	bs = x@binding_sites
	if(nrow(bs) < 1) {
	  warning('provided object contains no predicted TFBS hits, please ensure calling prediction prior to this method...')
	  return()
	}
	sequence.ids = sapply(x@sequences,attr,'uid')
  
	if(missing(pwms)) {
		pwms = unique(bs$pwm)
    # remove empty entries from list of pwms (0)
		pwms = setdiff(pwms,'0')
		cat('PWM id list not specified, including all PWMs with hits (',paste(pwms,sep=','),')...\n')
	} else {
		# check if list of specified pwms contains only valid pwm ids
		pwm.names = names(x@pfm)
		if( length(setdiff(pwms,pwm.names)) > 0 ) {
			warning('specified list of PWMs contains unknown PWM/PFM identifiers: ',paste(setdiff(pwms,pwm.names),sep=','),'! returning...')
      return()
		}
		cat('specified PWM id list! Including ',paste(pwms,sep=','),'...\n')
	}

	if(length(pwms) < 2) {
	  warning('At least 2 PWMs are required to draw a heatmap, provided PWMS:', paste(pwms, collapse=', '))
    return()
	}
  #holds sequence ids used for creating the venn diagram
	hits.lst = vector(mode='list')
	
  #get list of hits for each pwm
	for(pwm.id in pwms)
		hits.lst[[pwm.id]] = unique(bs$seqObj_uid[which(bs$pwm == pwm.id)])
  
	#include sequences without hits as 'None'
	if(include.empty.seqs)
		hits.lst[['None']] =  setdiff(sequence.ids, unique(bs$seqObj_uid))

	if(length(hits.lst) == 0) {
		warning('no hits for binding sites')
    return()
	}

	overlap.mat = matrix(0, nrow=length(pwms), ncol=length(pwms), dimnames=list(pwms, pwms))
	overlap.sig = matrix(0, nrow=length(pwms), ncol=length(pwms), dimnames=list(pwms, pwms))
  
	for(pwm1 in pwms)
		for(pwm2 in pwms){
		overlap.mat[pwm1, pwm2] = length(intersect(hits.lst[[pwm1]], hits.lst[[pwm2]]))
		overlap.sig[pwm1, pwm2] = 1 - phyper(overlap.mat[pwm1, pwm2] - 1, length(hits.lst[[pwm1]]), length(sequence.ids), length(hits.lst[[pwm2]]), lower.tail=F)
		}

	overlap.char = overlap.sig
	overlap.char[overlap.sig < .05] = "*"
	overlap.char[overlap.sig < .01] = "**"
	overlap.char[overlap.sig >= .05] = ""

	heatmap.2(overlap.mat, col = heat.colors(256), Rowv=FALSE, Colv=FALSE, dendrogram="none", trace="none", cellnote=overlap.char, margins=c(10,10))

	return(overlap.sig)
}
)

###############################################################################
# plot.detrending
#
# Using a background to detrend the foreground curve to find candidate pairs in
# certain distances
# if a file path is specified in pdf.name, the diagram will be written into the specified file

setMethod("plot.detrending", signature(x="cobindr"),
function(x, pwm1, pwm2, bin_length=20, z_value=3, overlap=0, abs.distance=FALSE) {
  
  results = detrending(x, pwm1, pwm2, bin_length, overlap, abs.distance)
  
 # if detreding could not be preformed, return with NULL
  if(is.null(results)) {
    warning('provided object contains no predicted hits, please ensure calling prediction prior to this method...')
    return()
  }
	
  x = results[[1]]
	occurrences = results[[2]]
	bg_occurrences = results[[3]]
	yy = results[[4]]
	y = results[[5]]
	
	# check if candidate pairs were found
	#if (TRUE %in% (y > z_value*sd(y)) || TRUE %in% (y < -z_value*sd(y))) {
	if (TRUE %in% (y > z_value*sd(y))) {
		cat('Overrepresented distances were found for pair', pwm1, pwm2, '. Use method write.significant.pairs() to get details.\n')
	}
	else {
		cat('No overrepresented distances were found for pair', pwm1, pwm2, '.\n')
	}
	########################
	# ------- Plot ------- #
	########################
	par(oma=c(0,0,2,0))
	par(mfrow=c(2,2))

	# target
	plot(x, occurrences, type='o', col='blue', xlab='distance in bp', ylab='# of pairs', main='Foreground')
	ylim = c(par("usr")[3], par("usr")[4])
	lines(c(0, 0), ylim, type='l', col='black', lty=2)

	# background
	plot(x, bg_occurrences, type='o', col='darkgreen', xlab='distance in bp', ylab='# of pairs', main='Background')
	ylim = c(par("usr")[3], par("usr")[4])
	lines(c(0, 0), ylim, type='l', col='black', lty=2)

	# target + normalised background
	plot(x, occurrences, type='o', col='blue', xlab='distance in bp', ylab='# of pairs', ylim=range(c(occurrences, yy)), main='Foreground with normalized background')
	lines(x, yy, type='o', col='darkgreen', lty=2)
	ylim = c(par("usr")[3], par("usr")[4])
	lines(c(0, 0), ylim, type='l', col='black', lty=2)

	# detrending
	plot(x, y, type='o', col='blue', xlab='distance in bp', ylab='normalized # of pairs', ylim=range(c(y, 3*sd(y), -3*sd(y))), main='Detrended distance distribution')
	xx = c(par("usr")[1], par("usr")[2])
	lines(xx, c(0, 0), type='l', col='black', lty=2)
	lines(xx, c(z_value*sd(y), z_value*sd(y)), type='l', col='red', lty=1)
	ylim = c(par("usr")[3], par("usr")[4])
	lines(c(0, 0), ylim, type='l', col='black', lty=2)

	title(main=paste('Pair (', pwm1, ', ', pwm2, ')', sep = ''), outer=TRUE)
}
)

setMethod("plot.pairdistribution", signature(x="cobindr"),
function(x, pwm1, pwm2){
	
	pwms = setdiff(unique(x@binding_sites$pwm),'0')
	for(pwm in c(pwm1, pwm2)) {
    if(!pwm %in% pwms) {
	  warning('specified PWM ',pwm,'is not valid! returning...\n')
      return()
    }
  }

	# gets list of PWM pairs
	pairs = x@pairs
	if(length(x@pairs) == 0) {
		warning('no hits for pairs of ', pwm1, ' and ',pwm2, ' found. returning... \n')
		return()
	}

	char.pairs = as.character(pairs$pair)
	# find pwm pair occurrences
	pair.idx = grep(gsub('\\$','\\\\\\$',paste(pwm1,pwm2)), char.pairs)
	# and in opposite
	if(pwm1 != pwm2)
		pair.idx = c(pair.idx, grep(gsub('\\$','\\\\\\$',paste(pwm2,pwm1)),char.pairs))

	# get list of Seq_uids for plotting
	seq.idx = pairs$seqObj_uid[pair.idx]
	seq.fct = as.factor(summary(as.factor(seq.idx), maxsum=length(seq.idx) + 1)) # the counts of pairs should be factors

	plot(seq.fct, main=paste('Pair distribution for \n PWM',pwm1,'and',pwm2), xlab="#pairs per sequence", ylab="#sequences")
}
)

###############################################################################
# plot.pairdistance
# creates histogram plot of distances between pairs of TFs as specified by pwm1 and pwm2
# pwm1/pwm2 - names of TFs
# pdf.name - outpuf pdf file if specified
# breaks - number of breaks to separate the distance distribution into
# main - figure main
# xlab - figure xlab 
# ylab - figure ylab 
# background - flag allowing to plot forground or background distance distribution
setMethod("plot.pairdistance", signature(x="cobindr"),
function(x, pwm1, pwm2, breaks=50, main=NA, xlab=NA, ylab=NA, background=FALSE){
	pwms = setdiff(unique(x@binding_sites$pwm),'0')
	for(pwm in c(pwm1, pwm2)) {
    if(!pwm %in% pwms) {
			warning('specified PWM identifier ',pwm,' is not valid! returning...\n')
      return()
    }
	}
	
	# gets list of PWM pairs
	if(!background)
		pairs = x@pairs
	else
		pairs = x@bg_pairs
	
	if(length(pairs) < 1) {
		warning('no hits for pairs of ', pwm1, ' and ',pwm2, ' found. returning... \n')
    return()
	}
	
	char.pairs = as.character(pairs$pair)
	# find pwm pair occurrences
	pair.idx = grep(gsub('\\$','\\\\\\$',paste(pwm1,pwm2)), char.pairs)
	# and in opposite
	if(pwm1 != pwm2)
		pair.idx = c(pair.idx, grep(gsub('\\$','\\\\\\$',paste(pwm2,pwm1)),char.pairs))	
	if(is.na(main)) main = paste('Pair distance distribution for \n PWM',pwm1,'and',pwm2)
	if(is.na(xlab)) xlab = "distance btw. binding sites [bp]"
	if(is.na(ylab)) ylab = "# pairs with specific distance" 
	hist(pairs$distance[pair.idx], breaks = breaks, 
			main=main, xlab=xlab, 
			ylab=ylab)
}
)

###############################################################################
# plot.positionprofile
#
# provides position-wise profile plot over total number of predicted TFBS for each PWM over all input sequences
# A windowing is used to provide a smoother appreance, the window size can be adjusted with the window paramter
setMethod("plot.positionprofile", signature(x="cobindr"),
function(x, wind.len=50) {
	#get pwm names
	pwms = as.character(unique(x@binding_sites$pwm))
	#make position count matrix
	l = max( sapply(x@sequences, function(x) length(x@sequence)) )
	hit.mat = matrix(0,nrow=length(pwms),ncol=l)
	rownames(hit.mat) = pwms
	#iterate over hits and collect for each PWM      
	bs = x@binding_sites
    if (nrow(bs) < 1) {
    	warning('no hits for binding sites. returning ... \n ')
      return()
    }
	for(i in 1:nrow(bs)) {
		p = as.character(bs$pwm[i])	
		hit.mat[p,bs$start[i]] = hit.mat[p,bs$start[i]] + 1
	}
	# window-smoothing
	hit.mat.sm = matrix(0,nrow=length(pwms),ncol=l-wind.len)
	rownames(hit.mat.sm) = pwms
	if(nrow(hit.mat.sm) > 1 ) {
		for(i in 1:(l-wind.len)) hit.mat.sm[,i] = rowMeans(hit.mat[,i:(i+wind.len)])
	} else {
		for(i in 1:(l-wind.len)) hit.mat.sm[,i] = mean(hit.mat[,i:(i+wind.len)]) }
	# remove invalid pwm entry
	pwms.valid = as.character(setdiff(unique(x@binding_sites$pwm),'0'))
	if(length(which(rownames(hit.mat.sm) == '0'))>0) {
		hit.mat.sm.val = hit.mat.sm[-(which(rownames(hit.mat.sm) == '0')),]
	} else {
		hit.mat.sm.val = hit.mat.sm
	}
	#plotting resulting profiles
	cols = rainbow(length(pwms.valid))
	matplot(t(hit.mat.sm.val), type='l', lty=1, col=cols, lwd=2, 
			main='', xlab = 'position (bp)', ylab = 'average number of predicted TFBS')
	legend(x='topleft', legend=pwms.valid, col=cols, lty=1)
}
)

###############################################################################
# plot.positions
# 
# plots hits for each PWM on the individual sequence
# which sequences to plot can be specified by providing list of sequence identifiers seq.ids
# which PWMs to plot can be specified as list pwms
# total height of the plot can be adjusted via argument height 
# TODO: missing legend for PWM color
setMethod("plot.positions", signature(x="cobindr"),
function(x, seq.ids, pwms, main, order.seq=FALSE, wind.size = 400, frac = 10 ) {
  warn <- NA
  if(!require('genoPlotR')) warn <- 'package genoPlotR required for this function, returning...'
  if(!require('vcd')) warn <- 'package vcd required for this function, returning...'
#  if(!require('multicore')) warn <- 'package multicore required for this function, returning...'
  if(!is.na(warn)) {
    warning(warn)
    return(x)
  }
  
	if(missing(seq.ids)) sequence.ids = sapply(x@sequences,attr,'uid')
	else if(length(seq.ids)>length(x@sequences)) {
		warning('specified more sequences for plot (',length(seq.ids),') than used in prediction (',length(x@sequences),')! using all sequences instead...')
		sequence.ids = sapply(x@sequences,attr,'uid')
	} else sequence.ids =  seq.ids
	#check if list of pwms to plot is provided
	if(missing(pwms)) pwms = as.character(setdiff(unique(x@binding_sites$pwm),'0'))
	#check, if hits for pwms are found
    if(nrow(x@binding_sites) < 1) {
      warning('no hits for binding sites. returning ...\n ')
      return()
    }
	#assign colors to each PWM
	cols = rainbow(length(pwms))
	names(cols) = pwms
	#sequence lengths 
	seq.len = sapply(x@sequences, function(x) length(x@sequence))
	names(seq.len) = sapply(x@sequences, attr, 'uid')

  func.name <- 'get.positions'
  opts <- list(x@binding_sites, seq.len, cols)
  hits <- parallelize(func.name, sequence.ids, opts, n.cpu)
  names(hits) = as.character(sequence.ids)
	
	#if selected plot sequences ordered such that similar patterns in TFBS are shown together
	# can take long for many sequences due to clustering of large matrix...
	if(order.seq) {
		cat('ordering sequences by TF hit pattern similarity. This can take long...\n')
		#convert hits into int matrix to cluster sequences by hit pattern similarity
		hit.m = matrix(0, nrow=length(hits),ncol=seq.len)
		for(i in 1:length(hits)) 
			for(j in 1:nrow(hits[[i]])) {
				if(hits[[i]][j,'name']=='seq') next
				hit.m[ i, hits[[i]][j,'start'] ] = which(hits[[i]][j,'name']== pwms) 
			}
		#window hit matrix to speed up clustering
		wind.start = c(1, seq( from = wind.size/frac, to = seq.len[1]-wind.size, by = wind.size/frac ))
		hit.m.w = matrix(0, nrow=length(hits),ncol=0)
		#for every pwms separately do windowing / concat result to window matrix
		for(idx in 1:length(pwms)) {
			hit.tmp = matrix(0, nrow=length(hits),ncol=length(wind.start))
			for(i in 1:nrow(hit.m))
				hit.tmp[i,] = unlist(sapply(wind.start, function(x) length(which(hit.m[i,x:(x+wind.size)]==1)) ))
			hit.m.w = cbind(hit.m.w, hit.tmp)
		}
		#cluster sequences
		hcl.res = hclust(dist(hit.m.w),method='complete')
		hits = hits[order.dendrogram(as.dendrogram(hcl.res))]
		cat('Done reordering...\n')
	}
	
	if(missing(main)) main = 'predicted TFBS positions per sequence'
	grid.newpage()
	pushViewport(viewport(x=0, width=.8, just="left")) 
	plot_gene_map(hits, gene_type = "bars", main=main, plot_new=F,
			dna_seg_label_cex=.8)
	popViewport()
	pushViewport(viewport(x=.9, width=.2, just="left")) 
	gp = gpar()
	gp$cex = .8
	grid_legend(0 ,.8 ,labels=names(cols) , pch=replicate(length(cols),15), col=cols, gp=gp )
	popViewport()
}
)

# There is a NOTE for R CMD check: no visible binding for dna_seg - parallelize problem?

get.positions <- function(s.id, opts) { 
  bs <- opts[[1]]
  seq.length <- opts[[2]]
  colors <- opts[[3]]
  
  idx = which(bs$seqObj_uid==s.id)
  names = list('seq')
  start = strand = list(1)
  end = list(seq.length[s.id])
  col = list('white')
		
  for (i in idx) {
    if(bs$pwm[i] == '0') next
    names = c(names, as.character(bs$pwm[i]))
    start = c(start, bs$start[i])
    end = c(end, bs$end[i])
    strand = c(strand, bs$strand[i])
    col = c(col, colors[as.character(bs$pwm[i])])
  }
  #TODO: until known why strand can be 0,1,2 - fix strand to 1
  strand = replicate(length(strand),1)
  #make and store dna seg.
  df = data.frame(name = unlist(names), start = unlist(start), end = unlist(end), strand = unlist(strand), col = unlist(col))
  return(genoPlotR::dna_seg(df))
}
          

###############################################################################
# plot.positions
# 
# plots hits for each PWM on the individual sequence
# which sequences to plot can be specified by providing list of sequence identifiers seq.ids
# which PWMs to plot can be specified as list pwms
# total height of the plot can be adjusted via argument height
# 2013-01-09: Is this an internal function?
setMethod("plot.positions.simple", signature(x="cobindr"),          
function(x, seq.ids, pwms, main) {
	#check if list of segments to plot is provided
	if(missing(seq.ids)) sequence.ids = sapply(x@sequences,attr,'uid')
	else if(length(seq.ids)>length(x@sequences)) {
		warning('specified more sequences for plot (',length(seq.ids),') than used in prediction (',length(x@sequences),')! using all sequences instead...')
		sequence.ids = sapply(x@sequences,attr,'uid')
	} else sequence.ids =  seq.ids
	#check if list of pwms to plot is provided
	if(missing(pwms)) pwms = as.character(setdiff(unique(x@binding_sites$pwm),'0'))
	#check, if hits for pwms are found
	if(nrow(x@binding_sites) < 1) {
		warning('no hits for binding sites. returning ...\n ')
    return()
	}
	#assign colors to each PWM
	cols = rainbow(length(pwms))
	names(cols) = pwms
	#sequence lengths 
	seq.len = sapply(x@sequences, function(x) length(x@sequence))
	names(seq.len) = sapply(x@sequences, attr, 'uid') 
	
	if(missing(main)) main = 'predicted TFBS positions per sequence'

	par(mar=c(5, 5, 1, 8)) # make the margin large enough
	plot(seq.len, xlim = c(0,max(seq.len)), ylim = c(0,length(x@sequences)), type='n',
			xlab='Position [bp]', ylab='Input Sequence')
	for(i in 1:length(pwms)) {
		hit.idx = grep(gsub('\\$','\\\\$',pwms[i]), as.character(x@binding_sites$pwm))
		points(x@binding_sites$start[hit.idx], as.numeric(x@binding_sites$seqObj_uid)[hit.idx], 
				col=cols[i], bg=cols[i], pch=21, cex=.5)
	}
	par(xpd=NA)
	tmp.u<-par('usr')  # get the dimensions
	legend(tmp.u[2], tmp.u[3],   # set the legend right to the plot
			xjust=-0.1, yjust=0,  # offset, so it's not  right at the border
			legend=pwms , pch=replicate(length(cols),21), col=cols, pt.bg = cols, cex=.7)
}
)

###############################################################################
# plot.gc 
# 
# visualize GC content or CpG content of input sequences
# cpg - if TRUE, CpG content is shown, otherwise GC content
# sig.test - if TRUE, wilcoxon.test is performed per individual window against all windows in other sequence at the same position
# wind.size - window size for GC/CpG computation
# main - plot heading
# pdf.name - optional name of the pdf file
# height- optional height of pdf to be created
# hm.margin - optional argument providing the margin widths for the heatmap (if sig.test=FALSE)
# frac - determines the overlap between consecutive windows as fraction wind.size/frac
setMethod("plot.gc", signature(x="cobindr"),
function(x, seq.ids, cpg=F, wind.size=50, sig.test=F, hm.margin=c(4,10), frac = 10, n.cpu=NA) {

  warn = NA
#   if(!require('gplots'))  warn <- 'package gplots required for this function, returning...'
#  if(!require('multicore'))  warn <- 'package multicore required for this function, returning...'
  if(!require('RColorBrewer'))  warn <- 'package RColorBrewer required for this function, returning...'
  if(!is.na(warn)) {
    warning(warn)
    return(x)
  }
	# set number of cores to use, if not specified
	n.cpu <- determine.cores.option(n.cpu)
		
	#check if list of segments to plot is provided
	if(missing(seq.ids)) sequence.ids = names(x@sequences)
	else if(length(seq.ids)>length(x@sequences)) {
		warning('specified more sequences for plot (',length(seq.ids),') than used in prediction (',length(x@sequences),')! using all sequences instead...')
		sequence.ids = names(x@sequences)
	} else sequence.ids =  seq.ids
	
	if(length(unique(sapply(x@sequences, function(x) length(x@sequence)))) > 1) {
		warning('function supports only sequence sets with identical length! aborting...')
		return()
	}
	
	seq.len = max(sapply(x@sequences, function(x) length(x@sequence)))
	wind.start = c(1, seq( from = wind.size/frac, to = seq.len-wind.size, by = wind.size/frac ))
	gene.gcs = matrix(NA,nrow=length(sequence.ids),ncol=length(wind.start)) #init gc matrix
	rownames(gene.gcs) = sequence.ids
	colnames(gene.gcs) = wind.start
	pb = txtProgressBar(min = 0, max = length(sequence.ids), style = 3)
	cat('calculating GC / CpG content...\n')
	
	if(!cpg) ylab = 'GC content (%)'
	else ylab = 'CpG content (%)'
	#iterate over sequences
	for(i in 1:length(sequence.ids)) {
          seq.id = sequence.ids[i]
          seq = x@sequences[[seq.id]]
		#iterate over window indices
          if(!cpg){
            func.name <- 'wind.gc.content'
            opts <- list(seq, wind.size)
            gene.gcs[seq.id,] = unlist(parallelize(func.name, wind.start, opts, n.cpu))
          }
          else{
            func.name <- 'wind.cpg.content'
            opts <- list(seq, wind.size)
            gene.gcs[seq.id,] = unlist(parallelize(func.name, wind.start, opts, n.cpu))
            print(gene.gcs[seq.id,][[1]])
          }              
          setTxtProgressBar(pb, i)
	}
	close(pb)
	
	# wilcoxon test if requested (slow for large number of sequences)
	if(sig.test) {
          gc.pvals = gene.gcs
          func.name <- 'compare.samples'
          
          for(i in 1:nrow(gene.gcs)){
            opts <- list(gene.gcs, i)
            gc.pvals[i,] =  unlist(parallelize(func.name, 1:ncol(gene.gcs), opts, n.cpu))
          }
		#convert pvals to useful sizes
		gc.pvals.sc = -log2(gc.pvals)
		gc.pvals.sc = (gc.pvals.sc / (max(gc.pvals.sc)-min(gc.pvals.sc))) + .2
	}
	
	#switch to simpler sequence names
	if(missing(seq.ids))
		rownames(gene.gcs) = sapply(x@sequences,attr,'name')
	
	if(!sig.test) {
		cols = colorpanel(256,low='blue', mid='yellow', high='red')
		heatmap.2(gene.gcs, dendrogram='row', Colv=F, col=cols, trace='none', margins=hm.margin,
				main=ylab, xlab='Position (bp)', ylab='Sequence')
	} else #if requested plot significance 
	{
		x = as.integer(colnames(gene.gcs))
		cols = rainbow(length(sequence.ids))
		#if sig test, plot only dots, otherwise lines
		par(mar=c(5, 5, 1, 8)) # make the margin large enough
		matplot(x,t(gene.gcs),type='n',lty=1,pch=1, col=cols,
				main='', xlab='Location (bp)',ylab=ylab)
		for(i in 1:nrow(gene.gcs))
			points(x,gene.gcs[i,],pch=1, cex=gc.pvals.sc[i,], col=cols[i])
		#legend
		par(xpd=NA)
		tmp.u<-par('usr')  # get the dimensions
		legend(tmp.u[2], tmp.u[3],   # set the legend right to the plot
				xjust=-0.1, yjust=0,  # offset, so it's not  right at the border
				legend=sapply(sequence.ids,function(x) paste(substr(x,1,min(25,nchar(x))),'.',sep='') )
				,col=cols,lwd=1,cex=.5)
	}
}
)

compare.samples <-  function(x, opts){
  gcs <- opts[[1]]
  idx <- opts[[2]]
  wilcox.test(gcs[idx, x], gcs[-idx, x])$p.value
}

wind.gc.content <- function(w.idx, opts) {
  seq <- opts[[1]]
  wind.size <- opts[[2]]
  
  letterFrequency(subseq(seq@sequence, start=w.idx, end=(w.idx+wind.size)), letters="CG", as.prob=T)
}

wind.cpg.content <- function(w.idx, opts) {
  seq <- opts[[1]]
  wind.size <- opts[[2]]
  
  cont = letterFrequency(subseq(seq@sequence, start=w.idx, end=(w.idx+wind.size)), letters=c("C",'G'))
  gccont = dinucleotideFrequency(subseq(seq@sequence, start=w.idx, end=(w.idx+wind.size)))['CG']

  return( (gccont / prod(cont)) * wind.size )
}
