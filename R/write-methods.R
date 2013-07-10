# writeResults.R
# description: provides methods to write obejcts to files
#
# Author: stefan <kroeger@informatik.hu-berlin.de>
#
#
###############################################################################
# library(seqinr)


# generic function for writing a complete configuration object into an folder(=file)
setMethod ("write" , signature("configuration", file="character"),
function (x, file, ncolumns = NULL, append = FALSE, sep = " ") {
 	if(is.null(file)){
		file = removeSpecialCharacters(paste(x@experiment_description,gsub(":","-",gsub(" ","",x@id)),"config","yml",sep="."))
	}
	if(file.exists(file)){
		msg = cat('existing file ',file, ' was overwritten ...\n')
		warning(msg)
	}
	vals.lst = vector(mode='list')
	for(x.val in slotNames(x)) {
		if( length(slot(x, x.val)) > 1 )
			vals.lst = append(vals.lst,list(slot(x, x.val)))
		else if( length(slot(x, x.val)) == 1 )
			vals.lst = append(vals.lst,slot(x, x.val))
		else
			vals.lst = append(vals.lst,NA)
	}
	names(vals.lst) = slotNames(x)
	# convert configuration list to yaml 
	x.yml = as.yaml(vals.lst)
	# write to file
	res = try(write(x.yml, file = file))
	# catch errors
	if(inherits(res,"try-error")){
		cat('could not save configuration to file...\n')
		cat(res,'\n')
	}
})


# generic function for writing a complete run object into an folder(=file)
setMethod ("write" , signature("cobindr", file="character"),
function (x, file , ncolumns = NULL, append = FALSE, sep = " ") {
	# file is used as path that will be the taregt dir bfor all written files
	if(is.null(file) || file==""){
		file = getwd()
	}
# 	print(x)
	if(!file.exists(file)){
		dir.create(file, showWarnings = TRUE, recursive = FALSE, mode = "0755")
	}
	uid = slot(x, "uid")
	name = slot(x, "name")
	f = removeSpecialCharacters(paste(name,gsub(":","-",gsub(" ","",uid)),sep=".")) #general file part for all files
	file = paste(file, f,sep="/") #added path to file

	#write single parts of object:
	print(paste('writing sequences to:', paste(file, "sequences","fasta",sep="."), '...',sep=" "))
	write.sequences(x, "sequences", paste(file, "sequences","fasta",sep="."))
	print(paste('writting bg_sequences to:', paste(file, "bg_sequences","fasta",sep="."), '...',sep=" "))
	write.sequences(x, "bg_sequences", paste(file, "bg_sequences","fasta",sep="."))
	print(paste('writing tfbs pairs to:', paste(file, "pairs","bed",sep="."), '...',sep=" "))
	write.table(get.pairs(x), paste(file, "pairs","bed",sep="."),sep="\t")
	print(paste('writing binding sites to:', paste(file, "bindingsites","bed",sep="."), '...',sep=" "))
	write.bindingsites(x, paste(file, "bindingsites","bed",sep="."))
	print(paste('writing configuration to:',  paste(file, "config","yml",sep="."), '...',sep=" "))
	write(x@configuration, paste(file, "config","yml",sep="."))
	print(paste('writing description to:',  paste(file, "description","txt",sep="."), '...',sep=" "))
	write(x@desc, paste(file, "description","txt",sep="."))
	print(paste('wrote', class(x)[1],':', x@name))
})


# generic function for writing a complete SeqObj into an file
# based on write.fasta of library seqinr - no definition of generic method necessary
setMethod ("write.fasta",  signature(sequence="SeqObj"), 
function (sequences,names, file.out, open, nbchar) { 
	#fasta sequence
	seqs = as.character(slot(sequences, "sequence"))
	
	#additional information for fasta header
	uid = slot(sequences, "uid")
	name = slot(sequences, "name")
	if (exists(as.character(substitute(slot(sequences, "location"))))) {
		tmp_loc = slot(sequences, "location")
	} else {tmp_loc = NA}
	if(is.na(names)) {
		header = paste(slot(sequences, "uid"),
						slot(sequences, "name"),
						slot(sequences, "species"),
						tmp_loc,
						slot(sequences, "comment"),
						sep="|")
	} else {header=names}
	

	#create filename
	if (is.na(file.out)) {
		file.out = removeSpecialCharacters(paste(gsub(":","-",gsub(" ","",uid)), name,"fasta",sep="."))
		file.out = uniqueFilename(file.out)
	}

	# Write the sequences TO file:
	write.fasta(sequences = seqs, names = header, nbchar = nbchar, file.out = file.out, open=open)
})


#write output of findPairs to file
# setMethod ("write.pairs", signature(x="cobindr"),
# function(x, file=NULL, background=FALSE) {
# 	if(is.null(file)){
# 		file = removeSpecialCharacters(paste(x@name,gsub(":","-",gsub(" ","",x@uid)),"pairs","tsv",sep="."))
# 		#check if filename exists and change if so
# 		file = uniqueFilename(file)
# 	}
# 	if(exists(as.character(substitute(x@pairs))) && length(x@pairs)>0){
# 		seq_slot = "sequences"
# 		pair_slot = "pairs"
# 		bs_slot = "binding_sites"
# 		seq_type = "sequence_type"
# 		if (background) {
# 			seq_slot = "bg_sequences"
# 			pair_slot = "bg_pairs"
# 			bs_slot = "bg_binding_sites"
# 			seq_type = "bg_sequence_type"
# 		}
# 		if (any(grep('\\.bed', file))) {
# 			# get locations from orig. sequences and extract position of TFBS:
# 			seqUIDs = slot(x, pair_slot)$seqObj_uid
# 
# 			bed_locs = lapply(slot(x, seq_slot), function(i) c(i@uid,i@location))
# 			bed_locs = matrix(unlist(bed_locs),ncol=2, byrow=T)
# 			
# 			bed_locs = matrix(unlist(strsplit(bed_locs[match(seqUIDs, bed_locs[,1]),2],":")), ncol=3, byrow=T)
# 			
# 			bed_locs = data.frame(as.character(bed_locs[,1]),as.numeric(bed_locs[,2]),as.numeric(bed_locs[,3]))
# 
# 			# names
# 			bed_names = paste(slot(x, pair_slot)$seqObj_uid,gsub(' ' ,'-',slot(x, pair_slot)$pair),sep="-")
# 			# score = distance between pairs
# 			bed_distance = slot(x, pair_slot)$distance
# 			# strand - always set "+" 
# 			bed_strand = rep("+",dim(slot(x, pair_slot))[1])
# 			# thickstart
# 			bs1 = slot(x, bs_slot)[slot(x, pair_slot)$bs_uid1,] 
# 			bs2 = slot(x, bs_slot)[slot(x, pair_slot)$bs_uid2,] 
# 			bed_thickstart = sapply(1:dim(bs1)[1],function(i) min(bs1$start[i], bs1$end[i], bs2$start[i], bs2$end[i])+bed_locs[i,2])
# 			# thickend
# 			bed_thickend = sapply(1:dim(bs1)[1],function(i) max(bs1$start[i], bs1$end[i], bs2$start[i], bs2$end[i])+bed_locs[i,2])
# 			# blockcounts
# 			bed_bc = 2
# 			# blockSizes
# 			bed_bs =apply(data.frame(
# 				slot(x, bs_slot)[slot(x, pair_slot)$bs_uid1,]$end - slot(x, bs_slot)[slot(x, pair_slot)$bs_uid1,]$start,
# # 				slot(x,pair_slot)$distance,
# 				slot(x, bs_slot)[slot(x, pair_slot)$bs_uid2,]$end - slot(x, bs_slot)[slot(x, pair_slot)$bs_uid2,]$start),1, function(i) paste(i, collapse=","))
# 			# blockStarts
# 			bed_bstarts = apply(data.frame(
# # 				bed_locs[,2] + 	
# 				slot(x, bs_slot)[slot(x, pair_slot)$bs_uid1,]$start,
# # 				slot(x,pair_slot)$distance,
# # 				bed_locs[,2] + 
# 				slot(x, bs_slot)[slot(x, pair_slot)$bs_uid2,]$start ),1, function(i) paste(i, collapse=","))
# 			bed_rgb = rep("255,255,255", nrow(bed_locs))
# 			# compose final BED table
# 			bed.output = data.frame(bed_locs, 
# 						bed_names, 
# 						bed_distance, 
# 						bed_strand, 
# 						bed_thickstart, 
# 						bed_thickend, 
# 						bed_rgb
# # 						, bed_bc, bed_bs, bed_bstarts
# 			)
# 			colsBed = paste("#chr","start", "end"
# 								, "name", "distance", "strand"
# 								,"areaOfInterest_start", "areaOfInterest_end", "color",sep="\t")
# 			
# 			#if somethings runs out of defined ranges:
# 			bed.output[which(bed.output$end < bed.output$areaOfInterest_end),]$areaOfInterest_end = bed.output[which(bed.output$end < bed.output$areaOfInterest_end),]$end
# 			
# 			#first write BED-Header
# 			write(paste("track name=\"pairs_",x@name,"\"  description=\"detected pairs of motifs. ",x@desc,"\" useScore=1\n",colsBed, sep=""),file=file)
# 			#write data to BED file
# 			write.table(bed.output, file=file, sep="\t", quote=F,row.names=F,col.names=F,append=TRUE)
# 			print(paste('wrote pairs to: ',file,sep=""))
# 			
# 
# 		} else {	
# 			pair_tab = cbind(x@pairs, x@uid, x@name)
# 			colnames(pair_tab)[c(dim(pair_tab)[2]-1,dim(pair_tab)[2])] = c("uid","name")
# 			write.table(x@pairs, file=file, col.names=T, row.names=F, sep="\t", dec=".")
# 			print(paste('wrote pairs to: ',file,sep=""))
# 		}
# 	}else{
# 		cat('no pairs to write')
# 	}
# }
# )


#write bindings sites as table
setMethod ("write.bindingsites.table",  signature(x="cobindr"), 
function(x, file=NULL){
	if(length(x@binding_sites)==0){
		stop('object does not contain any binding sites')
	}
	if(is.null(file)){
		file = removeSpecialCharacters(paste(x@name,gsub(":","-",gsub(" ","",x@uid)),"bindingsites","csv",sep="."))
		file = uniqueFilename(file)
	}
	if(exists(as.character(substitute(x$binding_sites)))){
		bs.ranges = get.bindingsite.ranges(x)
		bs_tab = data.frame(as.character(x@binding_sites$seqObj_uid), 
			as.character(seqnames(bs.ranges)),
			as.character(slot(bs.ranges, "elementMetadata")$elementMetadata.pwm),
			start(slot(bs.ranges, "ranges")),
			end(slot(bs.ranges, "ranges")),
			width(slot(bs.ranges, "ranges")),
			as.character(slot(bs.ranges, "strand")),
			slot(bs.ranges, "elementMetadata")$elementMetadata.score,
			as.character(slot(bs.ranges, "elementMetadata")$elementMetadata.seq))
		
		colnames(bs_tab) = c("seqObj_uid", "seqObj_name", "pwm", "hit_start", "hit_end", "hit_width", "hit_strand", "hit_score", "hit_sequence")
		
		write.table(bs_tab, file=file, col.names=T, row.names=F, sep="\t", dec=".", quote=F)
		print(paste('wrote binding sites as table to: ',file,sep=""))
	}else{
		cat('Object contains no binding sites.\n')
	}
}
)


#write binding sites as BED file
setMethod ("write.bindingsites",  signature(x="cobindr"), 
function(x, file=NULL, background=FALSE){
	fname = "bindingsites"
	bs_slot = "binding_sites"
	seq_slot = "sequences"
	seq_type = "sequence_type"
	if (background) {
		fname = "bg_bindingsites"
		bs_slot = "bg_binding_sites"
		seq_slot = "bg_sequences"
		seq_type = "bg_sequence_type"
	}
	if (length(slot(x, bs_slot))==0) {
		print('object does not contain any binding sites')
	}
	else {
		if (is.null(file)) {
			file = removeSpecialCharacters(paste(x@name,gsub(":","-", gsub(" ","",x@uid)),x@configuration@species,fname,"bed",sep="."))
			file = uniqueFilename(file)
		}
		if (exists(as.character(substitute(slot(x, bs_slot)))) && length(slot(x, bs_slot)) > 0) {
			if(slot(slot(x, "configuration"), seq_type) =="geneid"){
				# get locations from orig. sequences and extract position of TFBS:
				locs = lapply(slot(x, seq_slot)[slot(x, bs_slot)$seqObj_uid], function(i) slot(i, "location"))
				locs = matrix(unlist(strsplit(as.character(locs),":")), ncol=4, byrow=T)
				locs = data.frame(locs[,1],as.numeric(locs[,2]),as.numeric(locs[,3]),locs[,4])
				
				# create BED formatted table
				bed.output = data.frame(
				locs[,1], 
				locs[,2] + slot(x, bs_slot)$start,
				locs[,2] + slot(x, bs_slot)$end,
				# name:
				paste(slot(x, bs_slot)$seqObj_uid,
				as.character(slot(x, bs_slot)$pwm),as.character(slot(x, bs_slot)$seq),as.character(slot(x, bs_slot)$source),sep="-"),
				# score ad strand
				as.character(slot(x, bs_slot)$score),
				ifelse (slot(x, bs_slot)$strand==1, "+","-"))
				
				colnames(bed.output) = c("# chr","start", "end", "source_sequence_uid-name","hit_score", "source_strand") #"seq", "start_pos", "end_pos", "identified_by")
				
				#first write BED-Header
				write(paste("track name=\"",x@name,"\"  description=\"",x@desc,"\" useScore=1",sep=""),file=file)
				#write data to BED file
				write.table(bed.output, file=file, sep="\t", quote=F,row.names=F,append=TRUE)
				print(paste('wrote binding sites to: ',file,sep=""))

			}  else if (slot(slot(x, "configuration"), seq_type)  =="chipseq") {
				# get locations from orig. sequences and extract position of TFBS:
				#extract seqObj information
				locs = lapply(slot(x, seq_slot), function(i) c(i@uid,i@location))
				locs = matrix(unlist(locs),ncol=2, byrow=T)
				
				seqs = slot(x, bs_slot)$seqObj_uid
 				locs = matrix(unlist(strsplit(locs[match(seqs, locs[,1]),2],":")), ncol=3, byrow=T)

				locs = data.frame(as.character(locs[,1]),as.numeric(locs[,2]),as.numeric(locs[,3]))
				
				# create BED formatted table
				bed.output = data.frame(
				locs[,1], 
				locs[,2] + slot(x, bs_slot)$start,
				locs[,2] + slot(x, bs_slot)$end,
				# name:
				paste(slot(x, bs_slot)$seqObj_uid,
				as.character(slot(x, bs_slot)$pwm),as.character(slot(x, bs_slot)$seq),as.character(slot(x, bs_slot)$source),sep="-"),
				# score ad strand
				as.character(slot(x, bs_slot)$score),
				ifelse (slot(x, bs_slot)$strand==1, "+","-"))
				
				colsBed = paste("# chr","start", "end", "source_sequence_uid-name","hit_score", "source_strand",sep="\t")
				
				# write output to file:
				write(paste("track name=\"",x@name,"\"\n#description=\"",bs_slot,"; ",x@desc,"\"\n",colsBed, sep=""),file=file)
				# write data to BED file
				write.table(bed.output, file=file, sep="\t", quote=F,row.names=F, col.names=F,append=TRUE)	
				print(paste('wrote binding sites to: ',file,sep=""))
				
			} else {
				bed.output = data.frame(
				slot(x, bs_slot)$seqObj_uid,
				slot(x, bs_slot)$start,
				slot(x, bs_slot)$end,
				as.character(slot(x, bs_slot)$pwm), 
				slot(x, bs_slot)$score,
				ifelse (slot(x, bs_slot)$strand==1, "+","-"),
				as.character(slot(x, bs_slot)$seq),
				as.character(slot(x, bs_slot)$source))

				colnames(bed.output) = c("# source_sequence_uid", "hit_start","hit_end", "pwm", "hit_score", "source_strand", "pwm_sequence", "identified_by")

				# write output to file:
				write(paste("track name=\"",x@name,"\"\n#description=\"",bs_slot,"; ",x@desc,"\"",sep=""),file=file)
				# write data to BED file
				write.table(bed.output, file=file, sep="\t", quote=F,row.names=F, col.names=T,append=TRUE)	
				print(paste('wrote binding sites to: ',file,sep=""))
			}
		}
	}
}
)
           
#write the sequences of a cobindr object (object@sequences) into fasta file:
#writeSequences = function(object, slotname= "sequences", file=NULL){
setMethod ("write.sequences",  signature(x="cobindr"), 
function(x, slotname= "sequences", file=NULL){
	if(length(slot(x, slotname))> 0){
		uid = slot(x, "uid")
		name = slot(x, "name")
		if(is.null(file)){
			file = removeSpecialCharacters(paste(name,gsub(":","-",gsub(" ","",uid)),slotname,"fasta",sep="."))
			file = uniqueFilename(file)
		}
		
		# write all seqObj of list into a single file:
		seqs = lapply(slot(x, slotname), function(seqObj) as.character(slot(seqObj, "sequence")))
		fasta_header = lapply(slot(x, slotname), function(seqObj) {
 			uid = slot(seqObj, "uid")
			name = slot(seqObj, "name")
			if (exists(as.character(substitute(slot(seqObj, "location"))))) {
				tmp_loc = slot(seqObj, "location")
			} else {tmp_loc = NA}
			header = paste(slot(seqObj, "uid"),
							slot(seqObj, "name"),
							slot(seqObj, "species"),
							tmp_loc,
							slot(seqObj, "comment"),
							sep="|")
			} )
		names(seqs) = fasta_header

		# Write the sequences TO file:
		write.fasta(sequences = seqs, names = names(seqs), nbchar = 80, file.out = file, open="w")
	
	} else cat("Error, object is empty!\n") 
}           
)
           

# create a unique file for a given string
# 2013-01-08: This method does not belong to a class
#             This should be an internal method.
uniqueFilename = function(file){
	i=0
	parts = unlist(strsplit(file,"\\."))
	while(file.exists(file)){
		prefix = paste(parts[1:length(parts)-1], collapse = '.')
		suffix = parts[length(parts)]
		i=i+1
		file = paste(prefix,"_",i,".",suffix,sep="")
	}
	return(file)
}


# remove all special character from string to avoid problems
# 2013-01-08: This method does not belong to a class
#             This should be an internal method
removeSpecialCharacters = function(string){
	string = gsub(" ","",string)
	string = gsub("[^[:alnum:]\\._-]","_",string)
	return(string)
}


# function returns binding site as GRanges object (package: GenomicFeatures) which is easy to 
# write to bed or gff files ( export(get.bindingsite.ranges(object), "tfbs_hits.gff3") )
setMethod ("get.bindingsite.ranges",  signature(x="cobindr"), 
function(x) {
  if(sum(sapply(c('GenomicFeatures'), require, character.only=T))!=1) {
    warning('package GenomicFeatures required for this function, returning...')
    return()
  }
   #get binding sites from object
  bs.df = x@binding_sites
  #convert strand into expecetd symbols
  bs.df$strand = rep('+', nrow(bs.df))
  bs.df$strand[which(x@binding_sites$strand==2)] = '-'
  #assign external sequence names
  bs.df$sequence_names = sapply(bs.df$seqObj_uid, function(uid) slot(x@sequences[[which(sapply(x@sequences, slot, 'uid')==uid)]], 'name') 
 )
  # assign sequence lengths
  sequence_lengths =  sapply(x@sequences, function(s) nchar(slot(s, 'sequence')))
  names(sequence_lengths) = sapply(x@sequences, function(s) slot(s, 'name'))
	# create GRanges object
  res = with(bs.df, GRanges(seqnames=sequence_names,
                           ranges=IRanges(start, end, names=uid),
                            strand=strand,seqlengths=sequence_lengths,
                           elementMetadata=data.frame(pwm=pwm, score=score, source = source, seq=seq)) )
  return(res)
}
)


# # write.significant.pairs writes all binding sites of found significant pairs into a
# # tsv file. If locations for the sequences are given, binding site locations
# # will be mapped to these locations.
# setMethod ("write.significant.pairs", signature(x="cobindr"),
# function(x, pwm1, pwm2, out.file, bin_length=20, z_value=3, overlap=0, abs.distance=FALSE) {
#     results = detrending(x, pwm1, pwm2, bin_length, overlap, abs.distance)
#     if(is.null(results)) return() #if detrending failed stop here
#     
#     x_coord = results[[1]]
#     occurrences = results[[2]]
#     bg_occurrences = results[[3]]
#     yy = results[[4]]
#     y = results[[5]]
#     
#     # check if candidate pairs were found
#     #if (TRUE %in% (y > z_value*sd(y)) || TRUE %in% (y < -z_value*sd(y))) {
#     if (TRUE %in% (y > z_value*sd(y))) {
# 	    pair_bs = matrix(nrow=0, ncol=9)
# 		colnames(pair_bs) = c('seqUID', paste('seq1(', pwm1, ')', sep=''), 'loc1', 'strand1', paste('seq2(', pwm2, ')', sep=''), 'loc2', 'strand2', 'distance', 'source')
# 	
#         candidates = matrix(nrow=0, ncol=5)    
#         colnames(candidates) = c('Pair', 'start', 'end', 'Z-value', 'source')
# 	
#         for (i in 1:length(x_coord)) {
#             if (y[i] > z_value*sd(y)) {
#                 # foreground
#                 more_results = writeDetrending_intern(pwm1, pwm2, x@sequences, abs.distance, x@pairs, x_coord[i], bin_length, overlap, x@binding_sites, pair_bs, candidates, y[i]/sd(y), 'FOREGROUND')
#                 
#                 pair_bs = more_results[[1]]
#                 candidates = more_results[[2]]
#             }
#             #else if (y[i] < -z_value*sd(y)) {
#                 # background
#             #    more_results = writeDetrending_intern(pwm1, pwm2, x@bg_sequences, abs.distance, x@bg_pairs, x[i], bin_length, overlap, x@bg_binding_sites, pair_bs, candidates, y[i]/sd(y), 'BACKGROUND')
#                 
#             #    pair_bs = more_results[[1]]
#             #    candidates = more_results[[2]]
#             #}
#         }
# 
#         # write binding_sites
#         if (nrow(candidates) > 0) {
#             filename1 = paste(out.file, 'binding_sites.tsv', sep='_')
#             filename2 = paste(out.file, 'candidate_pairs.tsv', sep='_')
#             write.table(pair_bs, file=filename1, sep="\t", row.names=F)
#             write.table(candidates, file=filename2, sep="\t", row.names=F)
#         }
#     }
#     else {
#         cat('No overrepresented distances were found for pair', pwm1, pwm2, '.\n')
#     }
# }
# )
