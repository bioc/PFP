# PFP-class
.check.PFP <- function(object){
  if(!is(object, "PFP")) stop("object has to be of class \"PFP\" ")
  errors <- character()
  if(!is.list(object@pathways_score))
    errors <- c(errors, "pathways_score must be a list object")
  if(!is.numeric(object@pathways_score[["PFP_score"]]))
    errors <- c(errors, "PFP_score must be a numeric")
  if(!is.numeric(object@pathways_score[["p_value"]]))
    errors <- c(errors, "p_value must be a a numeric")
  if(!is.data.frame(object@pathways_score[["random_score"]]))
    errors <- c(errors, "random_score must be a dataframe")
  if(!is.list(object@pathways_score[["genes_score"]]))
    errors <- c(errors, "genes_score must be a list")
  if(!is.data.frame(object@ref_net_info))
    errors <- c(errors, "ref_net_info must be a dataframe")
  if(length(errors) == 0){
    TRUE
  }else{
    errors
  }
}


#'\code{PFP-class}
#'
#'An S4 object for storing pathway fingerprint scores information.
#'
#' @slot pathways_score, a list contains PFP_score, p_value, random_score, genes_score.
#' PFP_score is a numeric score indicating the performance of a gene_list in some pathways.
#' p_value is a statistic test for the PFP_score.
#' random_score is used for the statistic test.
#' genes_score is the detail scores of every gene in the gene_list.
#' @slot ref_net_info, a data.frame, which contains the specific information of pathway networks.
#' Just be the same as \code{\link{net_info}} in FPFRefnet, including the index, id, name, group and species.
#' #'@section method:
#'    \itemize{
#'      \item{pathways_score, \code{signature(object = "PFP")}:
#'        extract the pathways score}
#'      \item{ref_net_info, \code{signature(object = "PFP")}:
#'        extract the pathway networks information}
#'      \item{PFP_score, \code{signature(object = "PFP")}:
#'        extract the PFP score}
#'      \item{p_value, \code{signature(object = "PFP")}:
#'        extract P value}
#'      \item{random_score, \code{signature(object = "PFP")}:
#'        extract the random score}
#'      \item{genes_score, \code{signature(object = "PFP", index=NULL,
#'      index_type = c("pathway_id","pathway_name","slice"))}:
#'        extract the genes score}
#'      \item{refnet_names, \code{signature(object = "PFP")}:
#'        extract the refnet names}
#'      \item{sub_PFP, \code{signature(object = "PFP", group_name = NULL,
#'      index = NULL, index_type = c("slice","pathway_id","pathway_name"))}
#'        subset of PFP object}
#'      \item{plot_PFP, \code{signature(object, type = "character",
#'      p_size = "numeric", l_size = 'numeric')}:
#'        plot the Pathway Fingerprint.}
#'      \item{rank_PFP \code{signature(object = "PFP", decreasing=TRUE)}
#'        sort the PFP score.}
#'      \item{show_PFP, \code{signature(object = "PFP")}:
#'        display methods for S4 classes PFP, see also}
#'        }
#'
#' @name PFP-class
#' @rdname PFP-class
#' @exportClass PFP
#' @seealso \code{\link{pathways_score-methods}},
#' \code{\link{PFP_score-methods}},
#' \code{\link{p_value-methods}}, \code{\link{random_score-methods}},
#' \code{\link{genes_score-methods}}, \code{\link{refnet_names-methods}},
#' \code{\link{sub_PFP-methods}}, \code{\link{plot_PFP-methods}},
#' \code{\link{rank_PFP-methods}}, \code{\link{show_PFP-methods}},
#'
setClass("PFP", slot = list(pathways_score = "list", ref_net_info = "data.frame"),
         prototype = list(pathways_score = NULL, ref_net_info = NULL),
                          validity = .check.PFP)

#' Basic pathway networks scores of \emph{PFP} class
#' This function can extract the details in pathway fingerprint scores.
#'@exportMethod pathways_score
#'@rdname pathways_score-methods
#'@name pathways_score-methods
#'@param object, \code{PFP} class
#'@aliases pathways_score pathways_score-method
#'@docType methods
#'@seealso \code{\link{PFP-class}}
#'@return details in pathway fingerprint scores.
setGeneric("pathways_score",
           function(object){standardGeneric("pathways_score")})
#' @rdname pathways_score-methods
#' @aliases pathways_score pathways_score-methods
setMethod("pathways_score",signature="PFP",
          function(object){
            object@pathways_score
          }
)


#' Basic network information of \emph{PFP} class
#' This function extract the detail information of reference pathway networks.
#'
#'@exportMethod ref_net_info
#'@rdname ref_net_info-methods
#'@name ref_net_info-methods
#'@param object, \code{PFP} class
#'@aliases ref_net_info ref_net_info-methods
#'@docType methods
#'@seealso \code{\link{PFP-class}}
#'@return detail information of reference pathway networks
setGeneric("ref_net_info",
           function(object){standardGeneric("ref_net_info")})
#' @rdname ref_net_info-methods
#' @aliases ref_net_info ref_net_info-methods
setMethod("ref_net_info",signature="PFP",
          function(object){
            object@ref_net_info
          }
)


#' The score of \emph{PFP}
#' This function can extract the PFP_score of PFP.
#'
#'@exportMethod PFP_score
#'@rdname PFP_score-methods
#'@name PFP_score-methods
#'@param object, \code{PFP} class
#'@aliases PFP_score PFP_score-methods
#'@docType methods
#'@seealso \code{\link{PFP-class}}
#'@return the PFP_score
setGeneric("PFP_score",
           function(object){standardGeneric("PFP_score")})
#' @rdname PFP_score-methods
#' @aliases PFP_score PFP_score-methods
setMethod("PFP_score",signature="PFP",
          function(object){
            object@pathways_score[["PFP_score"]]
          }
)


#' The P value of \emph{PFP}
#' This function can extract the result of statistical analysis
#'@exportMethod p_value
#'@rdname p_value-methods
#'@name p_value-methods
#'@param object, \code{PFP} class
#'@aliases p_value p_value-methods
#'@docType methods
#'@seealso \code{\link{PFP-class}}
#'@return Statistical test result of each pathway score
setGeneric("p_value",
           function(object){standardGeneric("p_value")})
#' @rdname p_value-methods
#' @aliases p_value p_value-methods
setMethod("p_value",signature="PFP",
          function(object){
            object@pathways_score[["p_value"]]
          }
)


#' Extract the randomized similarity score of \emph{PFP}
#' This function extract the randomized similarity score for statistical test.
#'@exportMethod random_score
#'@rdname random_score-methods
#'@name random_score-methods
#'@param object, \code{PFP} class
#'@aliases random_score random_score-methods
#'@docType methods
#'@seealso \code{\link{PFP}}
#'@return a data.frame, each row represents one random test result with the same number of gene_list.
setGeneric("random_score",
           function(object){standardGeneric("random_score")})
#' @rdname random_score-methods
#' @aliases random_score random_score-methods
setMethod("random_score",signature="PFP",
          function(object){
            object@pathways_score[["random_score"]]
          }
)


#' The score of genes in \emph{PFP} class
#' This function extract the detail scores of every gene in the gene_list by specific condition.
#'@exportMethod genes_score
#'@rdname genes_score-methods
#'@name genes_score-methods
#'@param object, \code{PFP} class
#'@param index, character, indicating the groups to subset.
#'@param index_type, "pathway_id","pathway_name","slice"
#'@aliases genes_score genes_score-methods
#'@docType methods
#'@seealso \code{\link{PFP-class}}
#'@return a named vector of numeric scores
setGeneric("genes_score",
           function(object,index=NULL,
                    index_type = c("pathway_id","pathway_name","slice"))
             {standardGeneric("genes_score")})
#' @rdname genes_score-methods
#' @aliases genes_score genes_score-methods
setMethod("genes_score",signature="PFP",
          function(object,index=NULL,index_type = c("pathway_id","pathway_name","slice")){
            index_type <- match.arg(index_type, c("pathway_id","pathway_name","slice"))
            ref_net_info <- object@ref_net_info
            if (is.null(index)){
              net_select <- ref_net_info
            }else{
              if (index_type == "slice"){
                if (max(index) > nrow(object@ref_net_info)){
                  stop("You input oversize slice!\n",
                       "The max pathway number is ",nrow(object@ref_net_info))
                }
                net_select <- ref_net_info[index,]
              }else{
                if (index_type == "pathway_id"){
                  match_tf <- match(index,ref_net_info$id,nomatch = 0)
                }else if (index_type == "pathway_name"){
                  match_tf <- match(index,ref_net_info$name,nomatch = 0)
                }
                match_tf <- match_tf[match_tf!=0]
                net_select <- ref_net_info[match_tf,]
              }
            }
            return(object@pathways_score[["genes_score"]][net_select$id])
          }
)



#' Names of basic networks
#' This function extract the reference pathway network names of PFP.
#'@exportMethod refnet_names
#'@rdname refnet_names-methods
#'@name refnet_names-methods
#'@param object, \code{PFPRefnet} class
#'@aliases refnet_names refnet_names-methods
#'@docType methods
#'@return a vector contains pathway names
setGeneric("refnet_names",
           function(object){standardGeneric("refnet_names")})
#' @rdname refnet_names-methods
#' @aliases refnet_names refnet_names-methods
setMethod("refnet_names",signature="PFP",
          function(object){
            object@ref_net_info$name
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
#'@param group_name, the group name in kegg
#'@param index, the index of pathway, NULL or a list contains slice/numeric, character, specifying elements to extract.
#'This parameter' length must be the same as \code{group_name}.
#'Default is \emph{NULL}, indicating extract all the networks of a group. See
#'\emph{details} for more information.
#'
#'@details This function help users to extract the specific networks PFPscores for
#'customized analysis, which could be of entire group PFP or some part of
#'a specific group PFP.
#'
#'Note, the \code{index} argument is only worked while the group_name argument is
#'consideration, which means group_name is not \emph{NULL}. And the length must be
#'the same as \code{group_name}. Default is \emph{NULL}, indicating extract the entire PFP.
#'
#'@param index_type, the index type,such as "slice","pathway_id","pathway_name"
#'@aliases sub_PFP sub_PFP-methods
#'@docType methods
#'@seealso \code{\link{PFP-class}}
#'@return a PFP object contains just the selected elements.
setGeneric("sub_PFP",
           function(object, group_name = NULL, index = NULL, index_type =
                      c("slice","pathway_id","pathway_name"))
             {standardGeneric("sub_PFP")})
#' @rdname sub_PFP-methods
#' @aliases sub_PFP sub_PFP-methods
setMethod("sub_PFP",signature="PFP",
          function(object, group_name = NULL,index = NULL, index_type = c("slice","pathway_id","pathway_name")){
            index_type <- match.arg(index_type, c("slice","pathway_id","pathway_name"))
            if (is.null(group_name)){
              group_name <- unique(object@ref_net_info$group)
            }

            ref_net_info <- object@ref_net_info
            group_vec <- as.vector(object@ref_net_info$group)
            group_select_info <- lapply(X = group_name,function(x)ref_net_info[x==group_vec,])

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
                group_size <- vapply(all_group_names,function(x)sum(x==ref_net_info$group),0)
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
                       ref_net_info=net_select))
          }
)


#' The show_PFP generic function
#'
#' Show a short summary for PFP object.
#'
#'@exportMethod show_PFP
#'@param object, \code{PFP} object
#'@docType methods
#'@rdname show_PFP-methods
#'@aliases show_PFP show_PFP-methods
setGeneric("show_PFP",
           function(object){standardGeneric("show_PFP")})
#' @rdname show_PFP-methods
#' @aliases show_PFP show_PFP-methods
setMethod("show_PFP", "PFP",
          function(object){
            group_name <- unique(object@ref_net_info$group)
            group_size <- vapply(group_name,function(x)sum(x==object@ref_net_info$group),0)
            print(paste0("The PFP object has the following ",length(group_name)," group(s)."))
            print(group_name)
            print("The pathway numbers in the group(s) are displayed in the following.")
            print(group_size)
            # print("\n",)
            print(paste0("The total number of pathways in the PFP object is ",nrow(object@ref_net_info)))
            # print("\n")
            print("The details of PFP scores are displayed in the following.")
            print(object@pathways_score[["PFP_score"]])
          }
)


#' Plot PFP results
#' Function for visualization PFP results.
#'@exportMethod plot_PFP
#'@rdname plot_PFP-methods
#'@name plot_PFP-methods
#'@param object, \code{PFP} class
#'@param type, types of the visaulization of \emph{PFP} object, 'matchstick', 'line','point'.
#'Default is 'matchstick'.
#'@param p_size, point size of plot, default is 1.
#'@param l_size, line size of plot, default is 0.5.
#'@aliases plot_PFP plot_PFP-methods
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
              PFP_refnet_group <- as.vector(object@ref_net_info$group)
              sim_df <- data.frame(PFP_score = PFP_score, group = PFP_refnet_group,
                                   refnet_index = 1:nrow(object@ref_net_info))
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



#' rank of the PFP object by the value of PFP_score.
#' Function for show the rank of PFP results.
#'@exportMethod rank_PFP
#'@rdname rank_PFP-methods
#'@name rank_PFP-methods
#'@param object, \code{PFP} class
#'@param decreasing, Sorting method, the default is deceasing
#'@aliases rank_PFP rank_PFP-methods
#'@docType methods
#'@seealso \code{\link{PFP-class}}
#'@return a ranked PFP object.
setGeneric("rank_PFP",
           function(object,decreasing=TRUE){standardGeneric("rank_PFP")})
#' @rdname rank_PFP-methods
#' @aliases rank_PFP rank_PFP-methods
setMethod("rank_PFP",signature="PFP",
          function(object,decreasing=TRUE){
            ref_net_info <- object@ref_net_info

            ref_net_info[["PFP_score"]] <- data.frame(object@pathways_score[["PFP_score"]])
            if (length(ref_net_info[["p_value"]])==nrow(object@ref_net_info)){
              ref_net_info[["p_value"]] <- data.frame(object@pathways_score[["p_value"]])
              ref_net_info <- ref_net_info[order(ref_net_info[,"group"],ref_net_info[,"PFP_score"],-ref_net_info[,"p_value"],decreasing = decreasing),]
            }else{
              ref_net_info <- ref_net_info[order(ref_net_info[,"group"],ref_net_info[,"PFP_score"],decreasing = decreasing),]
            }
            ref_net_info <- ref_net_info[c("index","id","name","group","species")]

            match_id <- as.vector(ref_net_info$id)
            PFP_score <- object@pathways_score[["PFP_score"]][match_id]
            genes_score <- object@pathways_score[["genes_score"]][match_id]
            if (length(ref_net_info[["p_value"]])==nrow(object@ref_net_info)){
              p_value <- object@pathways_score[["p_value"]][match_id]
              random_score <- object@pathways_score[["random_score"]][match_id]
            }else{
              p_value <- numeric()
              random_score <- data.frame()
            }
            return(new(Class = "PFP",
                       pathways_score=list(PFP_score=PFP_score,p_value=p_value,random_score=random_score,genes_score=genes_score),
                       ref_net_info=ref_net_info))
          }
)

