

.check.PFP <- function(object){
  if(!is(object, "PFP")) stop("object has to be of class \"PFP\" ")
  errors <- character()
  if(!is.list(object@pathways_score))
    errors <- c(errors, "pathways_score must be a list object")
  if (!is.numeric(object@pathways_score[["PFP_score"]]))
    errors <- c(errors, "PFP_score must be a numeric")
  if (!is.numeric(object@pathways_score[["p_value"]]))
    errors <- c(errors, "p_value must be a a numeric")
  if (!is.data.frame(object@pathways_score[["random_score"]]))
    errors <- c(errors, "random_score must be a dataframe")
  if (!is.list(object@pathways_score[["genes_score"]]))
    errors <- c(errors, "genes_score must be a list")
  if (!is.list(object@net_info))
    errors <- c(errors, "net_info must be a list")
  
  if(length(errors) == 0)
    TRUE
  else
    errors
}

#'\code{PFP-class}
#'
#'An S4 object for storing network fingerprint similarity score information.
#'
#'@slot raw_score, a numeric vector, network fingerprint based on reference
#'networks before standardization.
#'@slot randomized_score, a data frame, the permulated similarity score.
#'@slot standardized_score, a numeric vector, the final standardized network fingerprint.
#'@slot cluster, an \emph{APResult} list, more details see package **apcluster**,
#'each element provides a cluster information of a
#'biological network based on one reference networks.
#' #'@section method:
#'    \itemize{
#'      \item{perm_score, \code{signature(object = "PFP")}:
#'        extract the randomized similarity score}
#'      \item{cluster_info, \code{signature(object = "PFP")}:
#'        extract the cluster information}
#'      \item{sub_PFP, \code{signature(object = "PFP")}:
#'        subset of PFP object}
#'      \item{plot, \code{signature(object, type = "character", p_size = "numeric", l_size = 'numeric')}:
#'        plot PFP results}
#'      \item{show, \code{signature(object = "PFP")}:
#'        display methods for S4 classes PFP, see also
#'        \code{\link[methods]{show}}}
#'    }
#'
#' @name PFP-class
#' @rdname PFP-class
#' @exportClass PFP
#' @seealso \code{\link{show-methods}},
#' \code{\link{plot-methods}}, \code{\link{perm_score-methods}},
#' \code{\link{cluster_info-methods}}, \code{\link{sub_PFP-methods}}
#'
setClass("PFP", slot = list(pathways_score = "list", net_info = "data.frame"),
         prototype = list(pathways_score = NULL, net_info = NULL),
                          validity = .check.PFP)



setGeneric("pathways_score",
           function(object){standardGeneric("pathways_score")})
#' @rdname pathways_score-methods
#' @aliases pathways_score pathways_score-methods
setMethod("pathways_score",signature="PFP",
          function(object){
            object@pathways_score
          }
)


setGeneric("refnet_info",
           function(object){standardGeneric("refnet_info")})
#' @rdname refnet_info-methods
#' @aliases refnet_info refnet_info-methods
setMethod("refnet_info",signature="PFP",
          function(object){
            object@net_info
          }
)


setGeneric("PFP_score",
           function(object){standardGeneric("PFP_score")})
#' @rdname PFP_score-methods
#' @aliases PFP_score PFP_score-methods
setMethod("PFP_score",signature="PFP",
          function(object){
            object@pathways_score[["PFP_score"]]
          }
)

setGeneric("p_value",
           function(object){standardGeneric("p_value")})
#' @rdname p_value-methods
#' @aliases p_value p_value-methods
setMethod("p_value",signature="PFP",
          function(object){
            object@pathways_score[["p_value"]]
          }
)


setGeneric("random_score",
           function(object){standardGeneric("random_score")})
#' @rdname random_score-methods
#' @aliases random_score random_score-methods
setMethod("random_score",signature="PFP",
          function(object){
            object@pathways_score[["random_score"]]
          }
)


setGeneric("genes_score",
           function(object,index=NULL,index_type = c("pathway_id","pathway_name","slice")){standardGeneric("genes_score")})
#' @rdname genes_score-methods
#' @aliases genes_score genes_score-methods
setMethod("genes_score",signature="PFP",
          function(object,index=NULL,index_type = c("pathway_id","pathway_name","slice")){
            index_type <- match.arg(index_type, c("pathway_id","pathway_name","slice"))
            net_info <- object@net_info
            if (is.null(index)){
              net_select <- net_info
            }else{
              if (index_type == "slice"){
                if (max(index) > nrow(object@net_info)){
                  stop("You input oversize slice!\n",
                       "The max pathway number is ",nrow(object@net_info))
                }
                net_select <- net_info[index,]
              }else{
                if (index_type == "pathway_id"){
                  match_tf <- match(index,net_info$id,nomatch = 0)
                }else if (index_type == "pathway_name"){
                  match_tf <- match(index,net_info$name,nomatch = 0)
                }
                match_tf <- match_tf[match_tf!=0]
                net_select <- net_info[match_tf,]
              }
            }
            return(object@pathways_score[["genes_score"]][net_select$id])
          }
)


setGeneric("refnet_names",
           function(object){standardGeneric("refnet_names")})
#' @rdname refnet_name-methods
#' @aliases refnet_name refnet_name-methods
setMethod("refnet_names",signature="PFP",
          function(object){
            object@net_info[c("id","name","group")]
          }
)


#' subset of PFP object
#'
#' This function extract the subsets of PFP-class.
#'
#'@exportMethod sub_PFP
#'@rdname sub_PFP-methods
#'@name sub_PFP-methods
#'@param object, \code{PFP} class
#'@param i, numeric or character indicating the index or the names of the
#'reference network
#'@aliases sub_PFP sub_PFP-methods
#'@docType methods
#'@seealso \code{\link{PFP-class}}
#'@return an similar PFP object contain just the selected elements.

setGeneric("sub_PFP",
           function(object, group_name = NULL, index = NULL, index_type = c("slice","pathway_id","pathway_name")){standardGeneric("sub_PFP")})
#' @rdname sub_PFP-methods
#' @aliases sub_PFP sub_PFP-methods

setMethod("sub_PFP",signature="PFP",
          function(object, group_name = NULL,index = NULL, index_type = c("slice","pathway_id","pathway_name")){
            index_type <- match.arg(index_type, c("slice","pathway_id","pathway_name"))
            if (is.null(group_name)){
              group_name <- unique(object@net_info$group)
            }
            
            net_info <- object@net_info
            group_vec <- as.vector(object@net_info$group)
            group_select_info <- lapply(X = group_name,function(x)net_info[x==group_vec,])
            
            all_group_names <- unique(group_vec)
            tf <- match(group_name,all_group_names,nomatch = 0) != 0
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
                group_size <- vapply(all_group_names,function(x)sum(x==net_info$group),0)
                max_slice <- vapply(index,max,0)
                max_group_select <- group_size[group_name]
                group_if <- max_slice > max_group_select
                if (sum(group_if) > 0){
                  stop("You input oversize slice!\n",
                       "The max pathway number is in the following!\n",
                       paste0(group_name,","),"\n",
                       paste0(max_group_select,","))
                }
                net_select <- lapply(seq_len(length(group_name)),function(i)group_select_info[[i]][index[[i]],])
                net_select <- do.call('rbind',net_select)
              }else{
                group_select_info <- do.call("rbind",group_select_info)
                if (index_type == "pathway_id"){
                  match_tf <- match(unlist(index),group_select_info$id,nomatch = 0)
                  if (length(match_tf[match_tf==0])>0){
                    print("The following pathways can't be found!")
                    print(setdiff(unlist(index),unlist(group_select_info[match_tf[match_tf!=0],"id"])))
                  }
                }else if (index_type == "pathway_name"){
                  match_tf <- match(unlist(index),group_select_info$name,nomatch = 0)
                  if (length(match_tf[match_tf==0])>0){
                    print("The following pathways can't be found!")
                    print(setdiff(unlist(index),unlist(group_select_info[match_tf[match_tf!=0],"name"])))
                  }
                }
                match_tf <-match_tf[match_tf!=0]
                net_select <- group_select_info[match_tf,]
              }
            }
            pathway_select_ids <- as.vector(net_select$id)
            PFP_score <- object@pathways_score[["PFP_score"]][pathway_select_ids]
            p_value <- object@pathways_score[["p_value"]][pathway_select_ids]
            random_score <- object@pathways_score[["random_score"]][,pathway_select_ids]
            genes_score <- object@pathways_score[["genes_score"]][pathway_select_ids]
            return(new(Class = "PFP",
                       pathways_score=list(PFP_score=PFP_score,p_value=p_value,random_score=random_score,genes_score=genes_score),
                       net_info=net_select))
          }
)



setGeneric("show_PFP",
           function(object){standardGeneric("show_PFP")})
#' @rdname show_PFP-methods
#' @aliases show_PFP show_PFP-methods
#' The show generic function
#'
#' Show a shor summary for PFP object, see \code{\link[methods]{show}}.
#'
#'@exportMethod show
#'@param object, \code{PFP} object
#'@docType methods
#'@rdname show-methods
#'@aliases show show-methods
setMethod("show_PFP", "PFP",
          function(object){
            group_name <- unique(object@net_info$group)
            group_size <- vapply(group_name,function(x)sum(x==object@net_info$group),0)
            print(paste0("The PFP object has the following ",length(group_name)," group(s)."))
            print(group_name)
            print("The pathway numbers in the group(s) are displayed in the following.")
            print(group_size)
            # print("\n",)
            print(paste0("The total number of pathways in the PFP object is ",nrow(object@net_info)))
            # print("\n")
            print("The details of PFP scores are displayed in the following.")
            print(object@pathways_score[["PFP_score"]])
          }
)




#' Plot PFP results
#'
#' Function for visualization PFP results.
#'
#'@exportMethod plot_PFP
#'@rdname plot_PFP-methods
#'@name plot_PFP-methods
#'@param object, \code{PFP} class
#'@param type, types of the visaulization of \emph{PFP} object, point or line.
#'Default is point.
#'@param p_size, point size of plot, default is 2.
#'@param l_size, line size of plot, default is 0.5.
#'#'@aliases plot_PFP plot_PFP-methods
#'@docType methods
#'@seealso \code{\link{PFP-class}}

setGeneric("plot_PFP",
           function(object, type = c('matchstick', 'line','point'), p_size = 1, l_size = 0.5)
           {standardGeneric("plot_PFP")})
#' @rdname plot_PFP-methods
#' @aliases plot_PFP plot_PFP-methods

setMethod("plot_PFP",'PFP',
          function(object, type = c('matchstick', 'line','point'), p_size = 1, l_size = 0.5){
            type <- match.arg(type, c('matchstick', 'line','point'))
            if (class(object) == 'PFP'){
              PFP_score <- object@pathways_score[["PFP_score"]]
              PFP_refnet_group <- as.vector(object@net_info$group)
              sim_df <- data.frame(PFP_score = PFP_score, group = PFP_refnet_group,
                                   refnet_index = 1:nrow(object@net_info))
              network_num <- length(PFP_score)
              if(all(!is.na(PFP_score))){ # skip plot if sim is NA
                p <- ggplot(sim_df,aes(x = refnet_index, y = PFP_score))
                if(type == "point")
                  print(p + geom_point(size = p_size, aes(color = group)))
                if(type == "line")
                  print(p + geom_line(size = l_size, aes(color = group, group = 1)))
                if (type == 'matchstick')
                  print(p + geom_point(size = p_size, aes(color = group)) +
                          geom_segment(aes(xend = refnet_index, yend = 0, color = group),
                                       size = l_size))
              }
              else
                stop('PFP score must be NA')
            }
          }
)


setGeneric("rank_PFP",
           function(object,decreasing=TRUE){standardGeneric("rank_PFP")})
#' @rdname rank_PFP-methods
#' @aliases rank_PFP rank_PFP-methods

setMethod("rank_PFP",signature="PFP",
          function(object,decreasing=TRUE){
            net_info <- object@net_info
            
            net_info[["PFP_score"]] <- data.frame(object@pathways_score[["PFP_score"]])
            if (length(net_info[["p_value"]])==nrow(object@net_info)){
              net_info[["p_value"]] <- data.frame(object@pathways_score[["p_value"]])
              net_info <- net_info[order(net_info[,"group"],net_info[,"PFP_score"],-net_info[,"p_value"],decreasing = decreasing),]
            }else{
              net_info <- net_info[order(net_info[,"group"],net_info[,"PFP_score"],decreasing = decreasing),]
            }
            net_info <- net_info[c("index","id","name","group","species")]
          
            match_id <- as.vector(net_info$id)
            PFP_score <- object@pathways_score[["PFP_score"]][match_id]
            genes_score <- object@pathways_score[["genes_score"]][match_id]
            if (length(net_info[["p_value"]])==nrow(object@net_info)){
              p_value <- object@pathways_score[["p_value"]][match_id]
              random_score <- object@pathways_score[["random_score"]][match_id]
            }else{
              p_value <- numeric()
              random_score <- data.frame()
            }
            return(new(Class = "PFP",
                       pathways_score=list(PFP_score=PFP_score,p_value=p_value,random_score=random_score,genes_score=genes_score),
                       net_info=net_info))
          }
)





# library("ggplot2")
# load("~/文档/PFP/RData/PFP.RData")
# 
# sub0 <- sub_PFP(object = PFP,group_name = unique(net_info(PFP)$group)[c(2,4)])
# sub0 <- sub_PFP(object = PFP,group_name = unique(net_info(PFP)$group)[c(2,4)],index = list(1:5,1:5))
# sub0 <- sub_PFP(object = PFP,group_name = unique(net_info(PFP)$group)[c(2,4)],index = list(c("mmu01230","mmu01521"),c("mmu04152")),index_type = "pathway_id")
# sub0 <- sub_PFP(object = PFP,group_name = unique(net_info(PFP)$group)[c(2,4)],index = list(c("MAPK signaling pathway","Phosphatidylinositol signaling system"),c("Staphylococcus aureus infection","Leishmaniasis")),index_type = "pathway_name")
# 
# refnet_names(sub0)
# show(sub0)
# plot_PFP(sub0)
# ran0 <- random_score(sub0)
# gs <- genes_score(sub0)
# plot_PFP(PFP)
# rank_PFP <- rank_PFP(PFP,decreasing=F)
# plot_PFP(rank_PFP)






