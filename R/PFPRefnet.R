# PFPRefnet-class
.check.PFPRefnet <- function(object){
  if(!is(object, "PFPRefnet")) stop("object has to be of class \"PFPRefnet\" ")
  errors <- character()
  if(!is.list(object@network))
    errors <- c(errors, "network must be a list object")
  if(is.null(object@network))
    errors <- c(errors, 'you must create your network first')
  if(!length(object@group))
    errors <- c(errors, 'group slot of PFPRefnet must be specified')
  if(!is.character(object@group))
    errors <- c(errors, "group must be a character")
  if(is.null(object@name))
    errors <- c(errors,'name slot of PFPRefnet must be specified')
  if(!is.list(object@name))
    errors <- c(errors, "name must be a list")
  if(is.null(object@organism))
    errors <- c(errors, 'organism slot of PFPRefnet must be specified')
  if(!is.character(object@organism))
    errors  <- c(errors, "organism must be a character")
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
#'@slot Refnet, object of graphNEL list represents the basic networks, and each
#'elements contains a group of basic networks.
#'
#'@slot group, a character vector whose length is the same with \emph{Refnet},
#'the group names of basic networks.
#'
#'@slot name, names of the basic networks, with the same data structure
#'with \emph{Refnet}.
#'
#'@slot organism, character, indicating the activation organism of basic networks.

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
setClass("PFPRefnet", slot=list(network = "list", name = "list",
                                group = "character", organism = "character"),
         prototype = list(network = NULL, name = NULL, group = character(0),
                          organism = character(0) ),
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

#' Group information of \emph{PFPRefnet}
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
            group_name <- object@group
            group_num <- length(group_name)
            net_name <- object@name
            group_size <- lapply(object@name,length) %>% unlist
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
#'@return a list

setGeneric("refnet_name",
           function(object){standardGeneric("refnet_name")})
#' @rdname refnet_name-methods
#' @aliases refnet_name refnet_name-methods
setMethod("refnet_name",signature="PFPRefnet",
          function(object){
            object@name
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
           function(object, group_name, index = NULL){standardGeneric("subnet")})
#' @rdname subnet-methods
#' @aliases subnet subnet-methods
setMethod("subnet",signature="PFPRefnet",
          function(object, group_name, index = NULL){
            net <- net(object)
            names(net) <- group(object)$name
            ref_name  <- refnet_name(object)
            names(ref_name) <- group(object)$name
            len <- length(group_name)
            organism <- object@organism
            fun_char2num <- function(i,index,ref_name){
              if(is.character(index[[i]])){
                return(match(index[[i]],ref_name[[i]]))
              }else{
                return(index)
              }
            }
            index <- lapply(X = seq_len(length(index)),FUN = fun_char2num,index,ref_name)
            # for(i in seq_len(length(index))){
            #   if(is.character(index[[i]]))
            #     index[[i]] <- match(index[[i]],ref_name[[i]])
            # }
            #sub_net <- net[[group_name]]
            if (len){
              if(len > 1){
                if(!is.null(index)){
                  if(length(group_name) != length(index))
                    stop('the length of index must be equal to the selected group numbers')
                  if(!is.list(index))
                    stop('index must a list with the same length with group_name')
                  sub_net <- mapply(function(x,y)net[[x]][y], group_name, index, SIMPLIFY = FALSE)
                  net_name <- mapply(function(x,y)ref_name[[x]][y], group_name, index, SIMPLIFY = FALSE)
                }
                else{
                  sub_net <- llply(group_name,"[[", x = net)
                  net_name <- llply(group_name, "[[", x= ref_name)
                }
              }
              else{
                if(is.null(index)){
                  sub_net <- net[[group_name]] %>% list
                  net_name <- ref_name[[group_name]] %>% list
                }
                else{
                  sub_net <- net[[group_name]][index] %>% list
                  net_name <-  ref_name[[group_name]][index] %>% list
                }
              }
              return(new("PFPRefnet",network = sub_net, name = net_name, group = group_name,
                         organism = object@organism))
            }
            else
              return(NULL)
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
setMethod("show", "PFPRefnet",
          function(object){
            organism <- object@organism
            ref_net <- object@network
            group <- object@group
            group_len <- llply(ref_net,length) %>% unlist
            len <- sum(group_len)
            ref_net_name <- refnet_name(object)
            if (length(ref_net_name) == 1)
              show_net_name <- paste0(ref_net_name[[1]][1], ", ...")
            else
              show_net_name <- sapply(ref_net_name,"[", 1) %>% paste0(', ...')
            show_group <- data.frame(group_name = group, net_num = group_len,
                                     net_name = show_net_name)
            row.names(show_group) <- paste0("group",1:nrow(show_group))
            ## cat("Object of class ", class(object), "\n", sep = "")
            ## cat("\n")
            cat( "Basic networks of","organism", organism,"\n")
            cat("\n")
            cat(len, "basic networks;", "classification into",
                length(ref_net),"groups:","\n")
            print(show_group)
          }
)

