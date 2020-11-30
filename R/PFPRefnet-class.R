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
#'@slot net_info, a dataframe which contains the index, id, name, group and
#'species.It contains the information of the pathway networks, whose row number
#'is the same with \emph{network}.
#'
#' #'@section method:
#'    \itemize{
#'      \item{network, \code{signature(object = "PFPRefnet")}:
#'        extract networks of PFPRefnet}
#'      \item{net_info, \code{signature(object = "PFPRefnet")}:
#'        extract net information of PFPRefnet}
#'      \item{group, \code{signature(object = "PFPRefnet")}:
#'        extract group information}
#'      \item{net_names, \code{signature(object = "PFPRefnet")}:
#'        the names of basic networks}
#'      \item{subnet, \code{signature(object = "PFPRefnet")}:
#'        subset basic networks, e.g. a group of a networks or some networks of
#'        some given groups}
#'      \item{show_net, \code{signature(object = "PFPRefnet")}:
#'        display methods for S4 classes PFPRefnet, see also
#'        \code{\link{show_net}}}
#'    }
#'
#'@name PFPRefnet-class
#'@rdname PFPRefnet-class
#'@exportClass PFPRefnet
#'@seealso \code{\link{network-methods}}, \code{\link{net_info-methods}},
#'\code{\link{group-methods}}, \code{\link{net_names-methods}},
#'\code{\link{subnet-methods}}, \code{\link{show_net-methods}},
#'@return a object of PFPRefnet class
#'@examples
#' data(PFPRefnet_hsa)
#' PFPRefnet_hsa
setClass("PFPRefnet", slot=list(network = "list", net_info = "data.frame"),
         prototype = list(network = NULL, net_info = NULL),
         validity = .check.PFPRefnet)


#' Basic pathway networks of \emph{PFPRefnet} class
#'
#' This function extract the basic networks of PFPRefnet class.
#'
#'@exportMethod network
#'@rdname network-methods
#'@name network-methods
#'@param object, \code{PFPRefnet} class
#'@aliases network network-methods
#'@docType methods
#'@seealso \code{\link{PFPRefnet-class}}
#'@return a graphNEL list of all basic networks
#'@examples
#' data(PFPRefnet_hsa)
#' network <- network(PFPRefnet_hsa)
setGeneric("network",
           function(object){standardGeneric("network")})
#' @rdname network-methods
#' @aliases network network-methods
setMethod("network",signature="PFPRefnet",
          function(object){
            object@network
          }
)


#' Basic pathway networks information of \emph{PFPRefnet} class
#'
#' This function extract the basic networks information of PFPRefnet class.
#'
#'@exportMethod net_info
#'@rdname net_info-methods
#'@name net_info-methods
#'@param object, \code{PFPRefnet} class
#'@aliases net_info net_info-methods
#'@docType methods
#'@seealso \code{\link{PFPRefnet-class}}
#'@return a dataframe contains basic networks' information
#'@examples
#' data(PFPRefnet_hsa)
#' net_info <- net_info(PFPRefnet_hsa)
setGeneric("net_info",
           function(object){standardGeneric("net_info")})
#'@rdname net_info-methods
#'@aliases net_info net_info-methods
setMethod("net_info",signature="PFPRefnet",
          function(object){
            net_info <- object@net_info
          }
)


#' group information of \emph{PFPRefnet}
#'
#' This function contains names of basic groups of the networks and group
#' number, as well as the size of each group
#'@exportMethod group
#'@rdname group-methods
#'@name group-methods
#'@param object, \code{PFPRefnet} class
#'@aliases group group-methods
#'@docType methods
#'@seealso \code{\link{PFPRefnet-class}}
#'@return a list contains names of basic groups of the networks and group
#'number, as well as the size of each group
#'@examples
#' data(PFPRefnet_hsa)
#' group <- group(PFPRefnet_hsa)
setGeneric("group",function(object){standardGeneric("group")})
#'@rdname group-methods
#'@aliases group group-methods
setMethod("group",signature="PFPRefnet",
          function(object){
            group_name <- unique(object@net_info$group)
            group_num <- length(group_name)
            group_size <- vapply(group_name,
                                 function(x)sum(x==object@net_info$group),0)
            return(list(name = group_name, num = group_num, size = group_size))
          }
)


#' Names of basic networks
#'
#' This function extract the network names of PFPRefnet.
#'
#'@exportMethod net_names
#'@rdname net_names-methods
#'@name net_names-methods
#'@param object, \code{PFPRefnet} class
#'@aliases net_names net_names-methods
#'@docType methods
#'@seealso \code{\link{PFPRefnet-class}}
#'@return a vector contains pathway names
#'@examples
#' data(PFPRefnet_hsa)
#' net_names <- net_names(PFPRefnet_hsa)
setGeneric("net_names",
           function(object){standardGeneric("net_names")})
#'@rdname net_names-methods
#'@aliases net_names net_names-methods
setMethod("net_names",signature="PFPRefnet",
          function(object){
            object@net_info$name
          }
)


#' Subset the basic networks
#'
#' Extract or Replace parts of the PFPRefnet.
#'
#'@exportMethod subnet
#'@rdname subnet-methods
#'@name subnet-methods
#'@param object, \code{PFPRefnet} class.
#'@param group_name, character, indicating the groups to subset.
#'@param index_type, character, the type pf index, which could be
#'"slice","id","name".
#'@param index, NULL or a list contains slice/numeric, character,
#'specifying elements to extract.
#'This parameter' length must be the same as \code{group_name}.
#'Default is \emph{NULL}, indicating extract all the networks of a group. See
#'\emph{details} for more information.
#'
#'@details This function help users to extract the specific networks for
#'customized analysis, which could be of entire group networks or some part of
#'a specific group networks.
#'
#'Note, the \code{index} argument is only worked while the group_name argument
#'is consideration, which means group_name is not \emph{NULL}. And the length
#'must be the same as \code{group_name}. Default is \emph{NULL}, indicating
#'extract the entire group basic networks.
#'@aliases subnet subnet-methods
#'@docType methods
#'@seealso \code{\link{PFPRefnet-class}}
#'@return sub the network
#'@examples
#' data(PFPRefnet_hsa)
#' subnet <- subnet(PFPRefnet_hsa)
setGeneric("subnet",function(object, group_name = NULL,
                             index = NULL,
                             index_type = c("slice",
                                            "pathway_id",
                                            "pathway_name"))
  {standardGeneric("subnet")})
#' @rdname subnet-methods
#' @aliases subnet subnet-methods
setMethod("subnet",signature="PFPRefnet",
          function(object, group_name = NULL, index = NULL, index_type =
                     c("slice","pathway_id","pathway_name")){
            index_type <- match.arg(index_type, c("slice",
                                                  "pathway_id",
                                                  "pathway_name"))
            if (is.null(group_name)){
              group_name <- group(object)$name
            }

            net_info <- object@net_info
            group_vec <- as.vector(object@net_info$group)
            group_select_info <- lapply(X = group_name,
                                        function(x)net_info[x==group_vec,])

            all_group_names <- unique(group_vec)
            tf <- match(group_name,all_group_names,nomatch = 0) != 0
            if (sum(tf) < length(group_name)){
              stop("Please input right group name(s)! You should choose one or
                   more in the following names.","\n",
                   paste0("\"",all_group_names,collapse = "\","),"\"")
            }

            if (is.null(index)){
              group_select_info <- do.call("rbind",group_select_info)
              net_select <- group_select_info
            }else{
              if (index_type == "slice"){
                if(length(group_name) != length(index))
                  stop('When the index_type is slice, the length of index must
                       be equal to the selected group numbers')
                if(!is.list(index))
                  stop('When the index_type is slice, index must be a list with
                       the same length with group_name')
                max_slice <- vapply(index,max,0)
                max_group_select <- group(object)$size[group_name]
                group_if <- max_slice > max_group_select
                if (sum(group_if) > 0){
                  stop("You input oversize slice!\n",
                       "The max pathway number is in the following!\n",
                       paste0(group_name,","),"\n",
                       paste0(max_group_select,","))
                }
                net_select <- lapply(seq_len(length(group_name)),
                                     function(i)group_select_info[[i]][index[[i]],])
                net_select <- do.call('rbind',net_select)
              }else{
                group_select_info <- do.call("rbind",group_select_info)
                if (index_type == "pathway_id"){
                  match_tf <- match(unlist(index),
                                    group_select_info$id,
                                    nomatch = 0)
                  if (length(match_tf[match_tf==0])>0){
                    print("The following pathways can't be found!")
                    print(setdiff(unlist(index),
                                  unlist(group_select_info[match_tf[match_tf!=0],
                                                           "id"])))
                  }
                }else if (index_type == "pathway_name"){
                  match_tf <- match(unlist(index),group_select_info$name,
                                    nomatch = 0)
                  if (length(match_tf[match_tf==0])>0){
                    print("The following pathways can't be found!")
                    print(setdiff(unlist(index),
                                  unlist(group_select_info[match_tf[match_tf!=0],
                                                           "name"])))
                  }
                }
                if (length(match_tf[match_tf==0])>0){
                  print("The following pathways can't be found!")
                  print(group_select_info[match_tf[match_tf==0],])
                }
                match_tf <-match_tf[match_tf!=0]
                net_select <- group_select_info[match_tf,]
              }
            }
            return(new(Class = "PFPRefnet",
                       network=object@network[as.vector(net_select$id)],
                       net_info=net_select))
          }
)


#' Show an Object
#'
#' show method short for PFPRefnet object
#'
#'@exportMethod show_net
#'@param object, \code{PFPRefnet} object
#'@docType methods
#'@rdname show_net-methods
#'@aliases show_net show_net-methods
#'@seealso \code{\link{PFPRefnet-class}}
#'@return show the network
#'@examples
#'data(PFPRefnet_hsa)
#'show_net(PFPRefnet_hsa)
setGeneric("show_net",
           function(object){standardGeneric("show_net")})
#' @rdname show_net-methods
#' @aliases show_net show_net-methods

setMethod("show_net", "PFPRefnet",
          function(object){
            group_name <- unique(object@net_info$group)
            group_size <- vapply(group_name,
                                 function(x)sum(x==object@net_info$group),0)
            print(group_size)
          }
)
