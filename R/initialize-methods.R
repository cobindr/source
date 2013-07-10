# file provides all initializing methods (as suggested by Bioconductor)
#
# Author: manu <manuela.benary@cms.hu-berlin.de>
#         rob  <r.lehmann@biologie.hu-berlin.de>
#######################################################################

# initializing method for the config-file
setMethod("initialize", "configuration",
          function(.Object, fname){
            # generating default-object
            default.fname = system.file('extdata', 'config_default.yml', package='cobindR')
            
            yml = yaml.load_file(default.fname)

            # set values based on default config-file
            for(cfg.val in names(yml))
              if(cfg.val %in% slotNames(.Object))
                if(is.numeric( slot(.Object, cfg.val)))
                  slot(.Object,cfg.val) = as.numeric(yml[[cfg.val]])
                else
                  slot(.Object, cfg.val) = yml[[cfg.val]]
                        
            # set internal parameters - date, uid and path of config-file
            .Object@date = as.character(Sys.Date())
            .Object@id = paste("Configuration", Sys.time(), sep = "_")
            .Object@path = default.fname
                                  
            # if no config-filename, throw warning (return default-object)
            if(missing(fname)){
              msg = sprintf("no config-file defined, generating configuration-object with default values")
              warning(msg)
            }

            # overwrite default-values via config-file
            else{
              cat('Reading the configuration file: ', fname, '\n')
              # checking if YAML-file is existent
              if(!file.exists(fname)) {
                msg = sprintf("specified input %s does not exist or is not readable", fname)
                warning(msg)
              } else {
              
                .Object@path = fname         # set path variable
                yml = yaml.load_file(fname)  # load configuration file

                # set values based on user-given file - throw warning if additional slots are in the config-file
                for(cfg.val in names(yml))
                  if(cfg.val %in% slotNames(.Object)){
                    if(is.numeric( slot(.Object, cfg.val)))
                      slot(.Object,cfg.val) = as.numeric(yml[[cfg.val]])
                    else
                      slot(.Object,cfg.val) = yml[[cfg.val]]
                  } else{
                    msg = sprintf('the configuration file contains a field %s which has no counterpart in the configuration object', cfg.val)
                    warning(msg)
                  }
              } # end if(!file.exists())

              # make bg_sequence_type consistent
              .Object@bg_sequence_type = tolower(.Object@bg_sequence_type)
              if (grepl('fasta', .Object@bg_sequence_type))
                .Object@bg_sequence_type = "fasta"
            
               # make sequence_type consistent
              .Object@sequence_type = tolower(.Object@sequence_type)                  # making lower cases for sequence_type
              .Object@sequence_type = gsub('[[:punct:]]', '', .Object@sequence_type)  # strip special characters
              .Object@sequence_type = gsub('[[:space:]]', '', .Object@sequence_type)  # strip white spaces
              # find "chip" to identify the chipseq call
              if (grepl('chip', .Object@sequence_type))
              .Object@sequence_type = "chipseq"
            }

            miss.slots = c()
            # check, if object has empty slots or the text "NA" -> throw warnings            
            for(obj.val in slotNames(.Object))
              if(length(slot(.Object, obj.val)) == 0)
                miss.slots = cbind(miss.slots, obj.val)
              else
                # if someone writes the text "NA" in his config-file, than this is parsed and the corresponding slot is set to NA
                if(length(slot(.Object, obj.val)) == 1 && slot(.Object, obj.val) == "NA"){
                  miss.slots = cbind(miss.slots, obj.val)
                  na.string = paste('NA', typeof(slot(.Object, obj.val)), '', sep='_')
                  slot(.Object, obj.val) = eval(parse(text=na.string))
                }
           

            if(!(length(miss.slots) == 0)){
              msg = sprintf('the slots %s are empty, use set() or check config-file', paste(miss.slots, collapse=", "))
              warning(msg)
            }
            
            .Object  # return new configuration object
          })

# initializing method for the cobindr
setMethod("initialize", "cobindr",
          function(.Object, conf, name, desc = ''){
            # test if conf is fully defined
            miss.slots = c()
            for(obj.val in slotNames(conf))
              if(length(slot(conf, obj.val))==0) 
                miss.slots = cbind(miss.slots, obj.val)
              else
                if(any(is.na(slot(conf, obj.val))))
                  miss.slots = cbind(miss.slots, obj.val)

            if(!(length(miss.slots) == 0)){
              msg = sprintf('the slots %s are empty, use set() or check config-file, returning ...', paste(miss.slots, collapse=", "))
              stop(msg)
            }  
              
            print("Creating a new experiment!")
            .Object@uid = paste("Experiment", Sys.time(), sep = "_")
            .Object@name = name
            .Object@sequences = read.sequences(conf)

            .Object@desc = desc
            .Object@configuration = conf
            .Object@pfm = read.pfm(conf)
            
            gc.t = testCpG(.Object)
            if(gc.t$result != FALSE)
              if(!is.null(gc.t))
                warning('heterogeneous GC content detected in input sequences - use plot.gc and testCpG(x, do.plot=T)')
            
            .Object
          })


# initializing method for the sequence object
setMethod("initialize", "SeqObj",
          function(.Object, seq, id, species, name, comment = NA, location = 'unknown'){
            .Object@uid  = id
            .Object@name = name
            .Object@species = species
            .Object@comment = comment
            .Object@location = location

            if (length(seq) == 0){
              msg = sprintf('the sequence %s is empty', id)
              stop(msg)
            }
            else
              .Object@sequence = seq

            .Object
          })

# configuration constructor
cobindRConfiguration <- function(fname = NULL) {
	if(is.null(fname))
		new('configuration')
	else
		new('configuration', fname = fname)
}
# cobindR constructor
cobindr <- function(conf = NULL, name = '', desc = '') {
	if(is.null(conf)) # if no configuration provided, make default one
		conf = new('configuration')
	new('cobindr', conf, name=name, desc=desc)
}
# SeqObj constructor
seqObj <- function(seq, id, species, name, comment = NA, location = 'unknown') {
	new('SeqObj', seq, id, species, name, comment, location = 'unknown')
}
