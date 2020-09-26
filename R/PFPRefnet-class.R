# PFPRefnet-class
.check.PFPRefnet <- function(object){
  if(!is(object, "PFPRefnet")) stop("object has to be of class \"PFPRefnet\" ")
  errors <- character()
  if(is.null(object@network))
    errors <- c(errors, 'you must create your network first')
  if(!is.list(object@network))
    errors <- c(errors, "network must be a list object")
  if(!length(object@net_info))
    errors <- c(errors, 'net_info slot of PFPRefnet must be specified')
  if(!is.data.frame(object@net_info))
    errors <- c(errors, "net_info must be a dataframe")
  if(length(errors) == 0){
    TRUE
  }else{
    errors
  }
}

#'\code{PFPRefnet-class}
#'
#'An S4 object for storing PFP reference network information.
#'
#'@slot network, object of graphNEL list represents the basic networks.
#'
#'@slot net_info, a dataframe which contains the index, id, name, group and species 
#'information of the networks. Its row number is the same with \emph{network}.
#'

#' #'@section method:
#'    \itemize{
#'      \item{net, \code{signature(object = "PFPRefnet")}:
#'        extract the basic networks}
#'      \item{group, \code{signature(object = "PFPRefnet")}:
#'        extract group information}
#'      \item{subnet, \code{signature(object = "PFPRefnet")}:
#'        subset basic networks, e.g. a group of a networks or same networks of
#'        a given group}
#'      \item{refnet_name, \code{signature(object = "PFPRefnet")}:
#'        the  names of basic networks}
#'      \item{show, \code{signature(object = "PFPRefnet")}:
#'        display methods for S4 classes PFPRefnet, see also
#'        \code{\link[methods]{show}}}
#'    }
#'
#' @name PFPRefnet-class
#' @rdname PFPRefnet-class
#' @exportClass PFPRefnet
#' @seealso \code{\link{show-methods}},
#' \code{\link{net-methods}}, \code{\link{refnet_name-methods}},
#' \code{\link{group-methods}}, \code{\link{subnet-methods}}
#'
setOldClass('igraph')
setClass("PFPRefnet", slot=list(network = "list", net_info = "data.frame"),
         prototype = list(network = NULL, net_info = NULL),
         validity = .check.PFPRefnet)

#' Basic networks of \emph{PFPRefnet} class
#'
#' This function extract the basic networks of PFPRefnet class.
#'
#'@exportMethod net
#'@rdname net-methods
#'@name net-methods
#'@param object, \code{PFPRefnet} class
#'@aliases net net-methods
#'@docType methods
#'@seealso \code{\link{PFPRefnet-class}}
#'@return a igraph list of all basic networks

setGeneric("net",
           function(object){standardGeneric("net")})
#' @rdname net-methods
#' @aliases net net-methods
setMethod("net",signature="PFPRefnet",
          function(object){
            net <- object@network
            return(net)
          }
)

#' group information of \emph{PFPRefnet}
#'
#' This function extract the group information PFP basic networks.
#'
#'@exportMethod group
#'@rdname group-methods
#'@name group-methods
#'@param object, \code{PFPRefnet} class
#'@aliases group group-methods
#'@docType methods
#'@seealso \code{\link{PFPRefnet-class}}
#'@return a list which contains the group number and names of basic networks, as
#'well as the size of each group

setGeneric("group",
           function(object){standardGeneric("group")})
#' @rdname group-methods
#' @aliases group group-methods
setMethod("group",signature="PFPRefnet",
          function(object){
            group_name <- unique(object@net_info$group)
            group_num <- length(group_name)
            group_size <- vapply(group_name,function(x)sum(x==object@net_info$group),0)
            return(list(name = group_name, num = group_num, size = group_size))
          }
)

#' Names of basic networks
#'
#' This function extract names of PFP basic networks.
#'
#'@exportMethod refnet_name
#'@rdname refnet_name-methods
#'@name refnet_name-methods
#'@param object, \code{PFPRefnet} class
#'@aliases refnet_name refnet_name-methods
#'@docType methods
#'@seealso \code{\link{PFPRefnet-class}}
#'@return a data.frame

setGeneric("refnet_name",
           function(object){standardGeneric("refnet_name")})
#' @rdname refnet_name-methods
#' @aliases refnet_name refnet_name-methods
setMethod("refnet_name",signature="PFPRefnet",
          function(object){
            object@net_info[c("id","name","group")]
          }
)

#' Subset the basic networks
#'
#' Extract or Replace parts of the PFP basic networks.
#'
#'@exportMethod subnet
#'@rdname subnet-methods
#'@name subnet-methods
#'@param object, \code{PFPRefnet} class.
#'@param group_name, character, indicating the groups to subset.
#'@param index, numeric, character or NA, indices specifying elements to extract. This
#'parameter only works while \code{group_name} is a length-one character.
#'Default is \emph{NULL}, indicating extract all the networks of a group. See
#'\emph{details} for more information.
#'
#'@details This function help users to extract the specific networks for
#'customized analysis, which could be of entire group networks or some part of
#'a specific group networks.subsequent analysis.
#'
#'Note, the \code{index} argument is only worked while one argument is
#'consideration, which means group_name is a length-one character. And default
#'is \emph{NULL}, indicating extract the entire group basic networks.
#'
#'@aliases subnet subnet-methods
#'@docType methods
#'@seealso \code{\link{PFPRefnet-class}}

setGeneric("subnet",
           function(object, group_name, index_type = c("pathway_id","pathway_name"),index = NULL){standardGeneric("subnet")})
#' @rdname subnet-methods
#' @aliases subnet subnet-methods
setMethod("subnet",signature="PFPRefnet",
          function(object, group_name, index = NULL, index_type = c("slice","pathway_id","pathway_name")){
            if (is.null(group_name)){
              group_name <- group(object)$name
            }
            
            net_info <- object@net_info
            group_vec <- as.vector(object@net_info$group)
            group_select_info <- lapply(X = group_name,function(x)net_info[x==group_vec,])
            
            all_group_names <- unique(group_vec)
            tf <- match(group_name,all_group_names,nomatch = 0) == 0
            if (sum(tf) < length(group_name)){
              stop("Please input right group name(s)! You should choose one or more in the following names.","\n",
                   paste0("\"",all_group_names,collapse = "\","),"\"")
            }
            
            if (is.null(index)){
              group_select_info <- do.call("rbind",group_select_info)
              net_select <- group_select_info
            }else{
              if (index_type == "slice"){
                if(length(group_name) != length(index))
                  stop('When the index_type is slice, the length of index must be equal to the selected group numbers')
                if(!is.list(index))
                  stop('When the index_type is slice, index must be a list with the same length with group_name')
                max_slice <- vapply(index,max,0)
                max_group_select <- group(object)$size[,group_name]
                group_if <- max_slice > max_group_select
                if (sum(group_if) > 0){
                  stop("You input oversize slice!\n",
                       "The max pathway number is in the following!\n",
                       paste0(group_name,","),"\n",
                       paste0(max_group_select,","))
                }
                net_select <- lapply(seq_len(length(group_name)),function(i)group_select_info[[i]][index[[i]],])
                net_select <- do.call('rbind',net_select)
              }else if(index_type == "pathway_id"|index_type == "pathway_name"){
                group_select_info <- do.call("rbind",group_select_info)
                ids_vec <- unlist(index)
                if (index_type == "pathway_id"){
                  match_tf <- match(unlist(index),group_select_info$id,nomatch = 0)
                }else if (index_type == "pathway_name"){
                  match_tf <- match(unlist(index),group_select_info$name,nomatch = 0)
                }
                match_tf <-match_tf[match_tf!=0]
                net_select <- group_select_info[match_tf,]
              }else{
                stop("index_type must be one of the (\"slice\",\"pathway_id\",\"pathway_name\")")
              }
            }
            return(new(Class = "PFPRefnet",network=object@network[as.vector(net_select$id)],net_info=net_select))
          }
)

#' Show an Object
#'
#' show method short for PFPRefnet object, see \code{\link[methods]{show}}
#'
#'@exportMethod show
#'@param object, \code{PFPRefnet} class
#'@docType methods
#'@rdname show_PFPRefnet-methods
#'@aliases show_PFPRefnet show_PFPRefnet-methods
#'
setGeneric("show",
           function(object){standardGeneric("show")})
setMethod("show", "PFPRefnet",
          function(object){
            group_name <- unique(object@net_info$group)
            group_size <- vapply(group_name,function(x)sum(x==object@net_info$group),0)
            print(group_size)
          }
)

save(list = c("net_info","graph_list"),file = "~/文档/PFP/RData/PFPRefnet_data.RData")
