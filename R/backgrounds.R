# backgrounds.R
# provides various methods to create background for subsequent sequence analysis
# 
# Author: rob <r.lehmann@biologie.hu-berlin.de>
# 		  yuehien
###############################################################################



# reads background sequences specified in the provided configuration object
# or generates background according to the specified method

# require('rtfbs')

##############################################################################
# 2013-01-09: Is this an internal function?
setMethod("generate.background", signature(x="cobindr"),
function(x, n.cpu = NA) {
    if(!require('rtfbs')) {
		warning('package rtfbs must be available when specified in configuration...')
		return(x)
	}

	print('creating background sequence...')
	
	res.seqs = list()
	# decide which input method to use:
	# cases: 1. list of fasta-files
	#        2. simulation using markov models
	#        3. simulation using ushuffle permutations
	#        4. list of geneids
	# check of vocabulary in constructor of configuration
	if (x@configuration@bg_sequence_type == "fasta"){
		cat('reading background sequences from specified fasta file(s)\n')
		x@bg_sequences = read.background.fasta(x@configuration)
	} else if (any(grep('markov', x@configuration@bg_sequence_type, ignore.case=TRUE))){
	    counter = 1
		# get model order
		txt = x@configuration@bg_sequence_type
		k = unlist(strsplit(txt,'\\.'))[2]
		n = unlist(strsplit(txt,'\\.'))[3]
		n.input = length(x@sequences)
		if(is.na(k) | nchar(k) == 0) {
			warning('could not determine degree of Markov model to use. Using default value 3!')
			k = 3
		}
		if(is.na(n) | nchar(n) == 0) {
			warning('could not determine number of background sequences. Using default of 5 * input sequences!')
			n = 5 * n.input
		} else 
			n = as.numeric(n)
		if(n < n.input) {
			warning('Number of background sequences too small (', n,')! Setting number of background sequences to equal the input sequences (', n.input, ')...')
			n = n.input
		}
		
		cat('simulating',n ,'background sequences using markov models with degree ',k ,'\n')
		#determine how many background sequences to simulate from each input seq.
		bg.n = rep(floor(n/n.input), n.input)
		bg.n[n.input] = bg.n[n.input] + n - sum(bg.n)  
		for (i in 1:n.input) {
		    # input sequence
			msSeq = ms(as.character(x@sequences[[i]]@sequence))
		    # define a Markov model of order k
            markovModel = build.mm(msSeq, k)
			for(j in 1:bg.n[i]) {
				# create random background sequence with Markov model
                          simSeq = simulate.ms(markovModel, lengths.ms(msSeq))
                          res.seqs[counter] = new('SeqObj', seq = DNAString(sequences.ms(simSeq)), id = as.character(counter), name = as.character(counter), species = 'random', comment = 'background sequence simulated by Markov model')
	            counter = counter + 1
			}
	    }
	    x@bg_sequences = res.seqs
	} else if (any(grep('ushuffle', x@configuration@bg_sequence_type, ignore.case=TRUE))){
		#determine correct ushuffle command
		if(.Platform$OS.type == "unix")
			ushuffle = 'ushuffle'
		else
			ushuffle = 'ushuffle.exe'
		# check if ushuffle is available via command line		
		us.check = try(system(ushuffle,intern=TRUE))
		if(inherits(us.check, 'try-error')) {
			stop('ushuffle not found! Please make sure that ushuffle can be called from the command line using the command "ushuffle"...')
		}
		# get length of elements to shuffle
		txt = x@configuration@bg_sequence_type
		k = unlist(strsplit(txt,'\\.'))[2]
		n = unlist(strsplit(txt,'\\.'))[3]
		n.input = length(x@sequences)
		if(is.na(k) | nchar(k) == 0) {
			warning('could not determine degree of Markov model to use. Using default value 3!')
			k = 3
		} else 
			k = as.numeric(k)
		if(is.na(n) | nchar(n) == 0) {
			warning('could not determine number of background sequences. Using default of 5 * input sequences!')
			n = 5 * n.input
		} else 
			n = as.numeric(n)
		if(n < n.input) {
			warning('Number of background sequences too small (', n,')! Setting number of background sequences to equal the input sequences (', n.input, ')...')
			n = n.input
		}		
		cat('simulating',n ,'background sequences using ushuffle and sequence segments of length',k ,'\n')
		#determine how many background sequences to simulate from each input seq.
		bg.n = rep(floor(n/n.input), n.input)
		bg.n[n.input] = bg.n[n.input] + n - sum(bg.n)

		# generate permutation for each input sequence
                func.name <- 'global.permutation'
                opts <- list(x@sequences, bg.n, k, ushuffle)
                bg.seqs <- parallelize(func.name, 1:n.input, opts, n.cpu)
                bg.seqs <- unlist(bg.seqs)

		for(i in 1:length(bg.seqs)) { 
			bg.seqs[[i]]@uid = as.character(i)
			bg.seqs[[i]]@name = as.character(i)
		}

		x@bg_sequences = bg.seqs
	} else if (x@configuration@bg_sequence_type == "geneid"){
		x@bg_sequences = read.sequences(x@configuration, background_scan=TRUE)
	} else if (x@configuration@bg_sequence_type == "chipseq"){
		x@bg_sequences = read.sequences(x@configuration, background_scan=TRUE)
	} else if (any(grep('local', x@configuration@bg_sequence_type, ignore.case=TRUE))){
        # get parameters of local shuffle
        txt = x@configuration@bg_sequence_type
        k = unlist(strsplit(txt,'\\.'))[2]
        n = unlist(strsplit(txt,'\\.'))[3]
        n.input = length(x@sequences)

        if(is.na(k) | as.numeric(k) <= 1) {
            warning('No or false window length (', k,') given. Window length should be at least 2. Using default value of 10!')
            k = 10
        } else
            k = as.numeric(k)
        if(is.na(n) | as.numeric(n) <= 0) {
            warning('No or false number for background sequences given. Using the number of input sequences!')
            n = n.input
        } else
            n = as.numeric(n)
        if(n < n.input) {
            warning('Number of background sequences too small (', n,')! Setting the number of background sequences to the number of input sequences (', n.input, ').')
            n = n.input
        }		
        cat('Simulating',n ,'background sequences using local shuffling with window length',k ,'\n')

        # start simulation
        func.name <- 'local.permutation'
        opts <- list(k)
        input_sequences = rep(x@sequences, length.out=n)
        
        sequences <- parallelize(func.name, input_sequences, opts, n.cpu)

        for(i in 1:length(sequences)) { 
            sequences[[i]]@uid = as.character(i)
            sequences[[i]]@name = as.character(i)
        }

        x@bg_sequences = sequences
	} else
	    warning('I did not understand the provided background sequence type (field bg_sequence_type): ',x@configuration@bg_sequence_type,
				'! Should be chipseq / geneid / ushuffle / markov - please consult the manual...\n')
	#return list of obtained sequences
	return(x)
}
)

####################################################################################
# 2013-01-09: Is this an internal function?

setMethod ("read.background.fasta", signature(x="configuration"),
function(x) {
	
	res.seqs = list()
	
	# obtain list of sequence files from which to load
	seq.sources = x@bg_sequence_source
	# iterate over input file and rad seqs
	seq.iterator = 0

	for(seqsource in seq.sources) {
		cat('reading file',seqsource,'...\n')
		fasta.seqs =  readDNAStringSet(seqsource)
		# iterate over all sequences in current file
		# TODO: define name from fasta sequence
		#      find location based on sequence or fasta information
		for(i in 1:length(fasta.seqs)) {
			id = paste(x@species,x@bg_sequence_type,x@bg_sequence_origin,Sys.time(),seq.iterator,sep='_')
			new.seq = new('SeqObj', 
					seq = fasta.seqs[[i]],
					id     = id,
					name     = names(fasta.seqs[i]),
					species = x@species,
					comment = '')
			res.seqs[id] = new.seq
			seq.iterator = seq.iterator+1
		}
	}
	return(res.seqs)
}
)

global.permutation <- function(count, opts) { # iterate over input seqs and generate bg seqs
  seq.objs <- list()
  seq.list <- opts[[1]]
  bg.n <- opts[[2]]
  k <- opts[[3]]
  ushuffle.cmd <- opts[[4]]
  
  u.cmd = paste(ushuffle.cmd,'-s',as.character(seq.list[[count]]@sequence),"-k", k, "-n", bg.n[count])
  seq = try(system(u.cmd,intern=TRUE))
  
  if(inherits(seq, 'try-error'))
    stop('cobindr encountered an error while calling ushuffle. Possible reason can be the length of the input sequence.\n')

  for(i in 1:bg.n)
    seq.objs[i] <- new('SeqObj', seq = DNAString(seq[i]), id = '1', name = '1', species = 'random', comment = 'background sequence simulated by ushuffle')

  return(seq.objs)
}

local.permutation <- function(seq, opts) {
  nn = nchar(seq@sequence)
  k = opts[[1]]

  if (nn%%k == 0) {
    rn = c(c(replicate(as.integer(nn/k), sample(1:k))))
  }
  else {
    rn = c(c(replicate(as.integer(nn/k), sample(1:k))), sample(1:(nn%%k)))
  }
  
  de = c(rep(0:(as.integer(nn/k)-1) * k, each = k), rep(as.integer(nn/k) * k, nn%%k))

  seq.obj <- new('SeqObj', seq = DNAString(seq@sequence[rn+de]), id = '1', name = '1', species = 'random', comment = 'background sequence simulated by local shuffling')

  return(seq.obj)
}
