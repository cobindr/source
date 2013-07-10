# file provides all method definitions
#
# Author: manu <manuela.benary@cms.hu-berlin.de>
#         rob  <r.lehmann@biologie.hu-berlin.de>
###############################################################################

setGeneric (
     name = "get.pairs",
     def  = function(x,...){standardGeneric("get.pairs")}
)

setGeneric (
     name = "write.bindingsites.table",
     def  = function(x,...){standardGeneric("write.bindingsites.table")}
            )

setGeneric (
     name = "write.bindingsites",
     def  = function(x,...){standardGeneric("write.bindingsites")}
            )

setGeneric (
     name = "write.sequences",
     def  = function(x,...){standardGeneric("write.sequences")}
            )

setGeneric (
     name = "get.bindingsite.ranges",
     def  = function(x,...){standardGeneric("get.bindingsite.ranges")}
            )
setGeneric (
     name = "get.significant.pairs",
     def  = function(x,...){standardGeneric("get.significant.pairs")}
            )
# Internal
setGeneric (
     name = "generate.background",
     def  = function(x,...){standardGeneric("generate.background")}
            )
# Internal
setGeneric (
     name = "read.background.fasta",
     def  = function(x,...){standardGeneric("read.background.fasta")}
            )
# Internal
setGeneric (
     name = "detrending",
     def  = function(x,...){standardGeneric("detrending")}
            )

setGeneric (
     name = "find.pairs",
     def  = function(x, background_scan = FALSE, n.cpu = NA){standardGeneric("find.pairs")}
            ) 
# Internal
# setGeneric (
#      name = "find.pairs_intern",
#      def  = function(x,...){standardGeneric("find.pairs_intern")}
#             ) 

setGeneric (
     name = "plot.tfbslogo",
     def  = function(x,...){standardGeneric("plot.tfbslogo")}
            ) 

setGeneric (
     name = "plot.tfbs.venndiagram",
     def  = function(x,...){standardGeneric("plot.tfbs.venndiagram")}
            )

setGeneric (
     name = "plot.tfbs.heatmap",
     def  = function(x,...){standardGeneric("plot.tfbs.heatmap")}
            ) 

setGeneric (
     name = "plot.detrending",
     def  = function(x,...){standardGeneric("plot.detrending")}
            )

setGeneric (
     name = "plot.pairdistribution",
     def  = function(x,...){standardGeneric("plot.pairdistribution")}
            )

setGeneric (
     name = "plot.pairdistance",
     def  = function(x,...){standardGeneric("plot.pairdistance")}
            ) 

setGeneric (
     name = "plot.positionprofile",
     def  = function(x,...){standardGeneric("plot.positionprofile")}
            )

setGeneric (
     name = "plot.positions",
     def  = function(x,...){standardGeneric("plot.positions")}
            ) 

setGeneric (
     name = "plot.positions.simple",
     def  = function(x,...){standardGeneric("plot.positions.simple")}
            ) 

setGeneric (
     name = "plot.gc",
     def  = function(x,...){standardGeneric("plot.gc")}
            )

# Internal
setGeneric (
     name = "read.pfm",
     def  = function(x,...){standardGeneric("read.pfm")}
            ) 

setGeneric (
     name = "predicted2pwm",
     def  = function(x,...){standardGeneric("predicted2pwm")}
            )

# Internal
setGeneric (
     name = "input.pwm",
     def  = function(x,...){standardGeneric("input.pwm")}
            )

setGeneric (
     name = "rtfbs",
     def  = function(x,...){standardGeneric("rtfbs")}
            ) 
# Internal
setGeneric (
     name = "rtfbs.intern",
     def  = function(x,...){standardGeneric("rtfbs.intern")}
            ) 

setGeneric (
     name = "search.gadem",
     def  = function(x,...){standardGeneric("search.gadem")}
            ) 
# TODO: check mistake in example
setGeneric (
     name = "search.pwm",
     def  = function(x,...){standardGeneric("search.pwm")}
            ) 
# Internal
setGeneric (
     name = "read.sequences",
     def  = function(x,...){standardGeneric("read.sequences")}
            ) 

setGeneric (
     name = "testCpG",
     def  = function(x,...){standardGeneric("testCpG")}
            ) 

## generic functions for accessor-methods
# cobindr object
setGeneric("uid", function(x) standardGeneric("uid"))
setGeneric("name", function(x) standardGeneric("name"))
setGeneric("sequences", function(x) standardGeneric("sequences"))
setGeneric("bg_sequences", function(x) standardGeneric("bg_sequences"))
setGeneric("experiment_description", function(x) standardGeneric("experiment_description"))
setGeneric("configuration", function(x) standardGeneric("configuration"))
setGeneric("pfm", function(x) standardGeneric("pfm"))
setGeneric("binding_sites", function(x) standardGeneric("binding_sites"))
setGeneric("bg_binding_sites", function(x) standardGeneric("bg_binding_sites"))
setGeneric("pairs", function(x) standardGeneric("pairs"))
setGeneric("bg_pairs", function(x) standardGeneric("bg_pairs"))
setGeneric("pairs_of_interest", function(x) standardGeneric("pairs_of_interest"))
# # SeqObj
setGeneric("uid", function(x) standardGeneric("uid"))
setGeneric("name", function(x) standardGeneric("name"))
setGeneric("species", function(x) standardGeneric("species"))
setGeneric("location", function(x) standardGeneric("location"))
setGeneric("comment", function(x) standardGeneric("comment"))
setGeneric("sequence", function(x) standardGeneric("sequence"))
# # configuration
setGeneric("id", function(x) standardGeneric("id"))
setGeneric("experiment_description", function(x) standardGeneric("experiment_description"))
setGeneric("sequence_source", function(x) standardGeneric("sequence_source"))
setGeneric("sequence_origin", function(x) standardGeneric("sequence_origin"))
setGeneric("sequence_type", function(x) standardGeneric("sequence_type"))
setGeneric("bg_sequence_source", function(x) standardGeneric("bg_sequence_source"))
setGeneric("bg_sequence_origin", function(x) standardGeneric("bg_sequence_origin"))
setGeneric("bg_sequence_type", function(x) standardGeneric("bg_sequence_type"))
setGeneric("species", function(x) standardGeneric("species"))
setGeneric("downstream", function(x) standardGeneric("downstream"))
setGeneric("upstream", function(x) standardGeneric("upstream"))
setGeneric("max_distance", function(x) standardGeneric("max_distance"))
setGeneric("pairs", function(x) standardGeneric("pairs"))
setGeneric("pfm_path", function(x) standardGeneric("pfm_path"))
setGeneric("threshold", function(x) standardGeneric("threshold"))
setGeneric("fdrThreshold", function(x) standardGeneric("fdrThreshold"))
# setGeneric("date", function(x) standardGeneric("date"))
setGeneric("path", function(x) standardGeneric("path"))
setGeneric("mart", function(x) standardGeneric("mart"))
setGeneric("pseudocount", function(x) standardGeneric("pseudocount"))
setGeneric("pValue", function(x) standardGeneric("pValue"))

# generic functions for replacement-methods
# cobindr object
setGeneric("uid<-", function(x, value) standardGeneric("uid<-"))
setGeneric("name<-", function(x, value) standardGeneric("name<-"))
setGeneric("sequences<-", function(x, value) standardGeneric("sequences<-"))
setGeneric("bg_sequences<-", function(x, value) standardGeneric("bg_sequences<-"))
setGeneric("experiment_description<-", function(x, value) standardGeneric("experiment_description<-"))
setGeneric("configuration<-", function(x, value) standardGeneric("configuration<-"))
setGeneric("pfm<-", function(x, value) standardGeneric("pfm<-"))
setGeneric("binding_sites<-", function(x, value) standardGeneric("binding_sites<-"))
setGeneric("bg_binding_sites<-", function(x, value) standardGeneric("bg_binding_sites<-"))
setGeneric("pairs<-", function(x, value) standardGeneric("pairs<-"))
setGeneric("bg_pairs<-", function(x, value) standardGeneric("bg_pairs<-"))
setGeneric("pairs_of_interest<-", function(x, value) standardGeneric("pairs_of_interest<-"))

# SeqObj
setGeneric("uid<-", function(x, value) standardGeneric("uid<-"))
setGeneric("name<-", function(x, value) standardGeneric("name<-"))
setGeneric("species<-", function(x, value) standardGeneric("species<-"))
setGeneric("location<-", function(x, value) standardGeneric("location<-"))
setGeneric("comment<-", function(x, value) standardGeneric("comment<-"))
setGeneric("sequence<-", function(x, value) standardGeneric("sequence<-"))

# configuration
setGeneric("id<-", function(x, value) standardGeneric("id<-"))
setGeneric("experiment_description<-", function(x, value) standardGeneric("experiment_description<-"))
setGeneric("sequence_source<-", function(x, value) standardGeneric("sequence_source<-"))
setGeneric("sequence_origin<-", function(x, value) standardGeneric("sequence_origin<-"))
setGeneric("sequence_type<-", function(x, value) standardGeneric("sequence_type<-"))
setGeneric("bg_sequence_source<-", function(x, value) standardGeneric("bg_sequence_source<-"))
setGeneric("bg_sequence_origin<-", function(x, value) standardGeneric("bg_sequence_origin<-"))
setGeneric("bg_sequence_type<-", function(x, value) standardGeneric("bg_sequence_type<-"))
setGeneric("species<-", function(x, value) standardGeneric("species<-"))
setGeneric("downstream<-", function(x, value) standardGeneric("downstream<-"))
setGeneric("upstream<-", function(x, value) standardGeneric("upstream<-"))
setGeneric("max_distance<-", function(x, value) standardGeneric("max_distance<-"))
setGeneric("pairs<-", function(x, value) standardGeneric("pairs<-"))
setGeneric("pfm_path<-", function(x, value) standardGeneric("pfm_path<-"))
setGeneric("threshold<-", function(x, value) standardGeneric("threshold<-"))
setGeneric("fdrThreshold<-", function(x, value) standardGeneric("fdrThreshold<-"))
# setGeneric("date<-", function(x, value) standardGeneric("date<-"))
setGeneric("path<-", function(x, value) standardGeneric("path<-"))
setGeneric("mart<-", function(x, value) standardGeneric("mart<-"))
setGeneric("pseudocount<-", function(x, value) standardGeneric("pseudocount<-"))
setGeneric("pValue<-", function(x, value) standardGeneric("pValue<-"))
