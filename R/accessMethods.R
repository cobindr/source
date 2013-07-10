# file provides all accessor and replacement methods for cobindR objects
#
# Author: rob  <r.lehmann@biologie.hu-berlin.de>
###############################################################################

# accessor methods

# cobindr object
setMethod("experiment_description", signature =  "cobindr", definition = function(x) x@desc)
setMethod("uid", signature =  "cobindr", definition = function(x) x@uid)
setMethod("name", signature =  "cobindr", definition = function(x) x@name)
setMethod("sequences", signature =  "cobindr", definition = function(x) x@sequences)
setMethod("bg_sequences", signature =  "cobindr", definition = function(x) x@bg_sequences)
setMethod("configuration", signature =  "cobindr", definition = function(x) x@configuration)
setMethod("pfm", signature =  "cobindr", definition = function(x) x@pfm)
setMethod("binding_sites", signature =  "cobindr", definition = function(x) x@binding_sites)
setMethod("bg_binding_sites", signature =  "cobindr", definition = function(x) x@bg_binding_sites)
setMethod("pairs", signature =  "cobindr", definition = function(x) x@pairs)
setMethod("bg_pairs", signature =  "cobindr", definition = function(x) x@bg_pairs)
setMethod("pairs_of_interest", signature =  "cobindr", definition = function(x) x@pairs_of_interest)

# SeqObj
setMethod("uid", signature =  "SeqObj", definition = function(x) x@uid)
setMethod("name", signature =  "SeqObj", definition = function(x) x@name)
setMethod("species", signature =  "SeqObj", definition = function(x) x@species)
setMethod("location", signature =  "SeqObj", definition = function(x) x@location)
setMethod("comment", signature =  "SeqObj", definition = function(x) x@comment)
setMethod("sequence", signature =  "SeqObj", definition = function(x) x@sequence)

# configuration
setMethod("id", signature =  "configuration", definition = function(x) x@id)
setMethod("experiment_description", signature =  "configuration", definition = function(x) x@experiment_description)
setMethod("sequence_source", signature =  "configuration", definition = function(x) x@sequence_source)
setMethod("sequence_origin", signature =  "configuration", definition = function(x) x@sequence_origin)
setMethod("sequence_type", signature =  "configuration", definition = function(x) x@sequence_type)
setMethod("bg_sequence_source", signature =  "configuration", definition = function(x) x@bg_sequence_source)
setMethod("bg_sequence_origin", signature =  "configuration", definition = function(x) x@bg_sequence_origin)
setMethod("bg_sequence_type", signature =  "configuration", definition = function(x) x@bg_sequence_type)
setMethod("species", signature =  "configuration", definition = function(x) x@species)
setMethod("downstream", signature =  "configuration", definition = function(x) x@downstream)
setMethod("upstream", signature =  "configuration", definition = function(x) x@upstream)
setMethod("max_distance", signature =  "configuration", definition = function(x) x@max_distance)
setMethod("pairs", signature =  "configuration", definition = function(x) x@pairs)
setMethod("pfm_path", signature =  "configuration", definition = function(x) x@pfm_path)
setMethod("threshold", signature =  "configuration", definition = function(x) x@threshold)
setMethod("fdrThreshold", signature =  "configuration", definition = function(x) x@fdrThreshold)
# setMethod("date", signature =  "configuration", definition = function(x) x@date)
setMethod("path", signature =  "configuration", definition = function(x) x@path)
setMethod("mart", signature =  "configuration", definition = function(x) x@mart)
setMethod("pseudocount", signature =  "configuration", definition = function(x) x@pseudocount)
setMethod("pValue", signature =  "configuration", definition = function(x) x@pValue)

# replacement methods

setReplaceMethod("id", signature = c("configuration", "character"),
      definition = function(x, value){
          x@id <- value
          x
      })
setReplaceMethod("experiment_description", signature = c("configuration", "character"),
    definition = function(x, value){
        x@experiment_description <- value
        x
    })
setReplaceMethod("sequence_source", signature = c("configuration", "character"),
     definition = function(x, value){
         x@sequence_source <- value
         x
     })
setReplaceMethod("sequence_origin", signature = c("configuration", "character"),
      definition = function(x, value){
          x@sequence_origin <- value
          x
      })
setReplaceMethod("sequence_type", signature = c("configuration", "character"),
      definition = function(x, value){
          x@sequence_type <- value
          x
      })
setReplaceMethod("species", signature = c("configuration", "character"),
    definition = function(x, value){
        x@species <- value
        x
    })
setReplaceMethod("downstream", signature = c("configuration", "numeric"),
     definition = function(x, value){
         x@downstream <- value
         x
     })
setReplaceMethod("upstream", signature = c("configuration", "numeric"),
      definition = function(x, value){
          x@upstream <- value
          x
      })

setReplaceMethod("max_distance", signature = c("configuration", "numeric"),
  definition = function(x, value){
      x@max_distance <- value
      x
  })
setReplaceMethod("pairs", signature = c("configuration", "character"),
   definition = function(x, value){
       x@pairs <- value
       x
   })
setReplaceMethod("pfm_path", signature = c("configuration", "character"),
    definition = function(x, value){
        x@pfm_path <- value
        x
    })
setReplaceMethod("threshold", signature = c("configuration", "numeric"),
     definition = function(x, value){
         x@threshold <- value
         x
     })
setReplaceMethod("fdrThreshold", signature = c("configuration", "numeric"),
      definition = function(x, value){
          x@fdrThreshold <- value
          x
      })
setReplaceMethod("bg_sequence_origin", signature = c("configuration", "character"),
     definition = function(x, value){
         x@bg_sequence_origin <- value
         x
     })
setReplaceMethod("bg_sequence_source", signature = c("configuration", "character"),
     definition = function(x, value){
         x@bg_sequence_source <- value
         x
     })
setReplaceMethod("bg_sequence_type", signature = c("configuration", "character"),
     definition = function(x, value){
         x@bg_sequence_type <- value
         x
     })
	

# setReplaceMethod("date", signature = c("configuration", "character"),
#  definition = function(x, value){
#      x@date <- value
#      x
#  })
setReplaceMethod("path", signature = c("configuration", "character"),
  definition = function(x, value){
      x@path <- value
      x
  })
setReplaceMethod("mart", signature = c("configuration", "character"),
   definition = function(x, value){
       x@mart <- value
       x
   })
setReplaceMethod("pseudocount", signature = c("configuration", "character"),
   definition = function(x, value){
       x@mart <- value
       x
   })
setReplaceMethod("pValue", signature = c("configuration", "numeric"),
    definition = function(x, value){
        x@pValue <- value
        x
    })

# # SeqObj
setReplaceMethod("uid", signature = c("SeqObj", "character"),
     definition = function(x, value){
         x@uid <- value
         x
     })
setReplaceMethod("name", signature = c("SeqObj", "character"),
definition = function(x, value){
    x@name <- value
    x
})
setReplaceMethod("species", signature = c("SeqObj", "character"),
 definition = function(x, value){
     x@species <- value
     x
 })
setReplaceMethod("location", signature = c("SeqObj", "character"),
  definition = function(x, value){
      x@location <- value
      x
  })
setReplaceMethod("comment", signature = c("SeqObj", "character"),
   definition = function(x, value){
       x@comment <- value
       x
})
setReplaceMethod("sequence", signature = c("SeqObj", "DNAString"),
     definition = function(x, value){
         x@sequence <- value
         x
     })

# cobindr
setReplaceMethod("uid", signature = c("cobindr", "character"),
definition = function(x, value){
    x@uid <- value
    x
})
setReplaceMethod("name", signature = c("cobindr", "character"),
 definition = function(x, value){
     x@name <- value
     x
 })
setReplaceMethod("sequences", signature = c("cobindr", "list"),
  definition = function(x, value){
      x@sequences <- value
      x
  })
setReplaceMethod("bg_sequences", signature = c("cobindr", "list"),
   definition = function(x, value){
       x@bg_sequences <- value
       x
})
setReplaceMethod("experiment_description", signature = c("cobindr", "character"),
definition = function(x, value){
    x@desc <- value
    x
})

setReplaceMethod("configuration", signature = c("cobindr", "configuration"),
 definition = function(x, value){
     x@configuration <- value
     x
 })
setReplaceMethod("pfm", signature = c("cobindr", "list"),
  definition = function(x, value){
      x@pfm <- value
      x
  })
setReplaceMethod("binding_sites", signature = c("cobindr", "data.frame"),
   definition = function(x, value){
       x@binding_sites <- value
       x
})
setReplaceMethod("bg_binding_sites", signature = c("cobindr", "data.frame"),
definition = function(x, value){
    x@bg_binding_sites <- value
    x
})
setReplaceMethod("pairs", signature = c("cobindr", "data.frame"),
 definition = function(x, value){
     x@pairs <- value
     x
 })
setReplaceMethod("bg_pairs", signature = c("cobindr", "data.frame"),
  definition = function(x, value){
      x@bg_pairs <- value
      x
  })
setReplaceMethod("pairs_of_interest", signature = c("cobindr", "factor"),
   definition = function(x, value){
       x@pairs_of_interest <- value
       x
})
