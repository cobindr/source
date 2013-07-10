# file provides IO methods for sequence data
# 
# Author: rob <r.lehmann@biologie.hu-berlin.de>
# Author: stefan <kroeger@informatik.hu-berlin.de>
###############################################################################

#################
# function definition

# reads sequences specified in the provided configuration object
# 2013-01-09: Is this an internal function?
# 2013-01-18: Yes!
setMethod("read.sequences", signature(x="configuration"),
function(x, background_scan = FALSE) {

	res.seqs = list()
	# obtain list of sequence files from which to load
	seq.sources = x@sequence_source
	seq.type = x@sequence_type
	seq.origin = x@sequence_origin
	# check if background sequences should be scanned
	if (background_scan) {
	    seq.sources = x@bg_sequence_source
	    seq.type = x@bg_sequence_type
	    seq.origin = x@bg_sequence_origin
	}
	
	# decide which input method to use:
	# cases: 1. list of fasta-files
	#        2. file with ensembl gene ids
	#        3. file with chipseq reads (?)
	# check of vocabulary in constructor of configuration
	if (seq.type == "fasta"){
          # iterate over input file and rad seqs
          seq.iterator = 0
          for(seqsource in seq.sources) {
            cat('reading file',seqsource,'...\n')
            fasta.seqs = readDNAStringSet(seqsource)
            
	    # iterate over all sequences in current file
	    # TODO: find location based on sequence or fasta information
            pb = txtProgressBar(min = 0, max = length(fasta.seqs), style = 3)

            for(i in 1:length(fasta.seqs)) {
              id = paste(x@species,seq.type,seq.origin,Sys.time(),seq.iterator,sep='_')
              new.seq = new('SeqObj', 
                seq = fasta.seqs[[i]],
                                        #uid     = id,
				id = as.character(seq.iterator),
				name = names(fasta.seqs[i]),
				species = x@species,
				comment = '')
				res.seqs[id] = new.seq
				seq.iterator = seq.iterator+1
				setTxtProgressBar(pb, i)
			}
		}
	}
	else if (seq.type == "geneid") {
		if(!require('BSgenome'))
			stop('package BSGenome must be available when gene-ids are specified in configuration file.')
		#set species information for BSgenome
		av.species = available.genomes()
		species.idx = grep(x@species, av.species, ignore.case=T)
		#load species genome
		if(!require(av.species[tail(species.idx, 1)], character.only=TRUE))
           stop('package', av.species[tail(species.idx, 1)], 'not installed - use biocLite for downloading')# using the most recent version

		tmp.species = unlist(strsplit(av.species[tail(species.idx, 1)],"\\."))[2]
		cat('retrieving sequences for gene ids in ',seq.sources,'...\n')
		gene.ids = readLines(seq.sources)
		# definition of the biomart object 
		species.data = paste(x@species, "gene", "ensembl", sep="_")
		ensembl = useMart(biomart=x@mart, dataset = species.data)
		# get all gene locations
		gene.ids.pos = getBM(attributes=c("ensembl_gene_id","chromosome_name", "start_position", "strand"),filters="ensembl_gene_id", values=gene.ids, mart=ensembl)
		gene.ids.pos$chromosome_name = paste("chr",gene.ids.pos$chromosome_name, sep="")
		gene.ids.pos$strand[(gene.ids.pos$strand == -1)] <-"-"
		gene.ids.pos$strand[(gene.ids.pos$strand == 1)] <-"+"
		
		#get sequences from BSgenome
		gene.ids.seqs = getSeq(get(tmp.species), gene.ids.pos$chromosome_name, start=gene.ids.pos$start_position-x@upstream, end=gene.ids.pos$start_position+x@downstream, strand=gene.ids.pos$strand,  as.character=FALSE)
		names(gene.ids.seqs) = gene.ids.pos$ensembl_gene_id
		genesToremove = NULL
		res.seqs = lapply(1:nrow(gene.ids.pos), function(i)
			if(gene.ids.pos$ensembl_gene_id[i] %in% names(gene.ids.seqs)){
				new.seq = try(new('SeqObj',
					seq = DNAString(as.character(gene.ids.seqs[i])),
					id     = as.character(i),
					species = x@species,
					name = as.character(names(gene.ids.seqs)[i]),
					location = paste(gene.ids.pos$chromosome_name[i], start=gene.ids.pos$start_position[i]-x@upstream, end=gene.ids.pos$start_position[i]+x@downstream, strand=gene.ids.pos$strand[i],sep=":"),
					comment = paste(seq.type,seq.origin,Sys.time(),sep='_')))
			} else { 
				warning('Wasn\'t  able to retrieve sequence of ',gene.ids.pos$ensembl_gene_id[i],'; discarding sequence for analysis...\n')
				
			}
		)
		#delete genes with biomart error from list
		res.seqs[grep("error", res.seqs)] = NULL
	}
	else if (seq.type == "chipseq"){
		## obtain bed file from which to  extract peak location, retrieve fasta seqs and load
		if(!require('BSgenome'))

			stop('package BSGenome must be available when ChipSeq is specified in configuration file.')
		av.species <- available.genomes()
		species.idx <- grep(x@species, av.species, ignore.case=T)
		
		if(length(species.idx) == 0)
			stop('Species \'', x@species, '\' not in package BSgenome available')
		
		if(!require(av.species[tail(species.idx, 1)], character.only=TRUE))
           stop('package', av.species[tail(species.idx, 1)], 'not installed - use biocLite for downloading')# using the most recent version
	
		# read input file and retrieve sequences
		seq.iterator = 0
        for(seqsource in seq.sources) {
			cat('reading file',seqsource,'...\n')
			bed.data = read.table(seqsource, header=FALSE, sep = "\t")
			bed.chr = as.character(bed.data[,1])
			bed.start = bed.data[,2]-x@upstream
			bed.end =  bed.data[,3]+x@downstream

			cat('retrieving sequence data ... may take a while ...\n')
			tmp.species = unlist(strsplit(av.species[tail(species.idx, 1)],"\\."))[2]
			bed.seq = getSeq(get(tmp.species), bed.chr,start=bed.start, end=bed.end, as.character=FALSE)

			#create sequence objects
			cat('creating sequence objects...\n')
			pb = txtProgressBar(min = 0, max = length(bed.seq), style = 3)
			tmp.lst <-list()
			for(i in 1:length(bed.seq)){
				id = paste(x@species,seq.type,seq.origin,Sys.time(),seq.iterator,sep='_')
				new.seq = new('SeqObj', 
				seq = DNAString(as.character(bed.seq[i])),
				id=as.character(seq.iterator),
				name = id,
				location = paste(bed.chr[i], bed.start[i], bed.end[i], sep=':'),
				species = x@species,
				comment = id)
				tmp.lst[id] = new.seq
				setTxtProgressBar(pb, i)
				seq.iterator = seq.iterator+1
			}
			res.seqs = c(res.seqs, tmp.lst)
		}
	}
	cat('ready retrieving sequences!')
	return(res.seqs)
}
)
          
###############################################################################

# testCpG
#
# diagnostical function - GC content and CpG content are clustered using 2D gaussian
# models (Mclust). FALSE is returned if more than one subgroups are found
# using the bayesian information criterion (BIC)
# if do.plot=TRUE, the results are visualized
setMethod("testCpG", signature(x="cobindr"),
function(x, max.clust=4, do.plot=F, n.cpu = NA) {
	best.clust = 1  # as a homogenous set of genes is best to analyse (definition of background)
        
	if(!require('parallel')) {
		warning('package parallel required for this function returning...')
		return(list(result=FALSE))
	}
	# set number of cores to use, if not specified
	n.cpu <- determine.cores.option(n.cpu)
	
	if(length(x@sequences) > 5000) {
		warnings('large dataset (',length(x@sequences),') found: disabling input evaluation due to GC/CpG content! Please perform manually...')
		return(list(result=TRUE))
	}

	# old variant:
	# 	gc = matrix(NA,nrow=length(x@sequences), ncol=2)
	# 	colnames(gc) = c('CpG','GC')
	# 	rownames(gc) = names(x@sequences)

	# 	func.name <- 'cpg.content'
	# 	opts <- list()
	# 	gc[,'CpG'] <- unlist(parallelize(func.name, x@sequences, opts, n.cpu))
	# 	func.name <- 'gc.content'
	# 	gc[,'GC'] <- unlist(parallelize(func.name, x@sequences, opts, n.cpu))
	
	#new variant
	func.name <- 'cpg.gc.content'
	opts <- list()
	gc <- parallelize(func.name, x@sequences, opts, n.cpu)
	gc =matrix(gc,nrow=length(x@sequences), ncol=2, byrow=T)
	
	gc= t(sapply(seq(1:length(x@sequences)), function(i) {
		seq = x@sequences[[i]]
		c(dinucleotideFrequency(slot(seq, "sequence"),as.prob=T)['CG'],
			letterFrequency(slot(seq,"sequence"),letters=c('CG'),as.prob=T))
	}))
	colnames(gc) = c('CpG','GC')
	rownames(gc) = names(x@sequences)

	#perform clustering
	gc.cls = Mclust(gc,G=c(1:max.clust))	
	if(do.plot){
		xhist = hist(gc[,'GC'], plot=FALSE)   # get histogram data
		yhist = hist(gc[,'CpG'], plot=FALSE)
		top = max(c(xhist$density, yhist$density))
		# prepare layout
		nf = layout(matrix(c(2,0,1,3),2,2,byrow=TRUE), c(3,1), c(1,3), TRUE)
		par(mar=c(6,6,1,1))
		mclust2Dplot(data = gc, parameters = gc.cls$parameters, z = gc.cls$z, what = "classification",identify=T) 
		par(mar=c(0,6,1,1))
		barplot(yhist$density, axes=FALSE, ylim=c(0, top), space=0)
		par(mar=c(6,0,1,1))
		barplot(xhist$density, axes=FALSE, xlim=c(0, top), space=0, horiz=TRUE)
	}
	return(list(result=(gc.cls$G <= best.clust),gc=gc))
}
)

cpg.content <- function(x, opts){
  dinucleotideFrequency(x@sequence,as.prob=T)['CG']
}
# 2013-03-11 internal function
gc.content <- function(x, opts){
  letterFrequency(x@sequence,letters=c('CG'),as.prob=T)
}
# 2013-03-15 internal function
#for parallel() it's more efficient to combine methods.
cpg.gc.content <- function(x, opts) {
	return(c(cpg.content(x, opts),gc.content(x,opts) ))
}
