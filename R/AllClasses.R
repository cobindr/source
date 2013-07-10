# file provides all class definitions
# includes: the class definition of the configuration object
#           the class definition for the sequence object
#           the class definition for the run object
# 
# Author: rob  <r.lehmann@biologie.hu-berlin.de>
#         manu <manuela.benary@cms.hu-berlin.de>
#         stefan <kroeger@informatik.hu-berlin.de>
###############################################################################

# object definition
# prototyping creates a useful configuration object (using the foxp3.txt file with gene ids

setClass("configuration", representation(
				id="character",               # unique id set when reading a configuration from file 
				experiment_description="character", # for free use
				sequence_source="character", 
				sequence_origin="character",    # source of sequence data, e.g. ensembl
				sequence_type="character",      # chipseq, or whatever
				bg_sequence_source="character", # file path or list of paths 
				bg_sequence_origin="character", # how the background sequence was obtained:simulated, curated or pulled out of ear
				bg_sequence_type="character",   # determines the reading / generation of the background model. values: NA, fasta, geneid 
				species = "character",          # if available
				downstream = "numeric", 
				upstream = "numeric",
				max_distance ="numeric",
				pairs ="character",				# specify list pairs of PWMs (by name as in the provided file(s) / database), separated by space  
                                pfm_path = "character",          # combining input of all pfms in one directory/list of files
				threshold="numeric",            # threshold for tfbs prediction
				fdrThreshold="numeric",         # threshold for rtfbs prediction
				date="character",               # set automatically when reading config file
				path = "character",				# holds path to the configuration file which was used to create the corresponding configuration instance
				mart = "character",              # optional mirror for biomart
				pseudocount = "numeric",               # pseudocount for detrending analysis
                                pValue = "numeric"             # for searching with RGadem
				)
         )

setClass("SeqObj",
         representation(uid      = "character",   # unique id for internal representation
                        name     = "character",   # biological reference name (p.e. gene name)
                        species   = "character",  # reference species
                        location = "character",   # location on the reference genome. e.g. chr1 start=490, end=510 -> location = chr1:490:510
                        comment  = "character",   # comments and notes
                        sequence = "DNAString"),  # the actual sequence
         prototype(uid      = paste(Sys.time(), sep=''),
                   name     = NA_character_,
                   species  = NA_character_,
                   location = 'unknown',
                   comment  = NA_character_,
                   sequence = DNAString('')),
         )

setClass("cobindr",
         representation(uid       = "character", # some unique id 
                        name      = "character", # name of the experiment
                        sequences = "list",      # list of sequence objects to analyze
                        bg_sequences = "list",   # list of background sequence objects
                        desc      = "character", # free chosen experiment description 
                        configuration = 'configuration', # the configuration used
                        pfm = 'list',            # keep only one list of transcription factors
                        binding_sites = "data.frame",    # predicted binding sites are saved here. Data frame structure: uid:integer, seqObj_uid:integer, pwm:factor, start:integer, end:integer, score:double, seq:character, strand:factor, source:factor
                        bg_binding_sites = "data.frame",    # background binding sites
                        pairs     = "data.frame",        # found pairs are saved here. Data frame structure: uid:integer, seqObj_uid:integer, pair:factor, orientation:factor, bs_uid1:integer, bs_uid2:integer, distance:integer
                        bg_pairs     = "data.frame",       # background pairs
                        pairs_of_interest = "factor"    # contains pairs for search
                        ))
