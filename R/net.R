# get_exp_cor_edges
#' @title get co-expression genes
#' @description compute the correlation coefficient of gene expression data,return the most related genes
#' @param gene_list, a vector of characters
#' @param data_std, a matrix of data, such as gene expression data, whose rownames are gene
#' names or ids and colnames are sample names
#' @param method, a chareater, which method to compare the correlation of gene expression data
#' it could be "pearson", "kendall", "spearman", "spearman" is default
#' @param num, an integer, the top number of co-expressed genes to choose, 5 is default
#' @param cor_threshold, a numeric, the threshold of the correlation coefficient to choose, default is \emph{NULL}
#' @details, This function computes the correlation coefficient of gene expression data between \code{gene_list}
#' and \code{data_std}, it will return a data.frame which can be translated a graph or network.
#' In the data.frame, \code{source} refers to the genes in \code{gene_list}, \code{target} refers to the top coexpressed genes,
#' \code{weight} refers to the correlated coefficient of genes in \code{source} and \code{target},
#' \code{pathway} is "uncertain" and \code{edge_type} is "coexp".
#'
#' Note, when choosing the top co-expressed genes, we will use the \code{num} param if
#' the \code{cor_threshold} param is \emph{NULL}. If not, we will choose the \code{cor_threshold} param.
#' @return the coexp of edges.
#' @examples
#' # Load the data
#' data(data_std)
#' data(gene_list_hsa)
#' # Get the coexp of edges
#' #edges_coexp <- get_exp_cor_edges(gene_list_hsa,data_std)
#' @export
get_exp_cor_edges <- function(gene_list,data_std,method="spearman",num=5,cor_threshold=NULL){
  bg_genelist <- rownames(data_std)
  match_list <- match(x = gene_list,table = bg_genelist,nomatch = 0)
  if (sum(match_list==0)){
    stop("gene_list must have common elements with rownames of data_std, maybe the gene id types you choose in the two data set are not consistent.")
  }else if(sum(match_list!=0)<length(match_list)){
    print("Some genes in gene_list can't be found in data_std. These genes will be removed.")
    print(gene_list[match_list==0])
    gene_list <- gene_list[match_list!=0]
  }

  get_gene_cor_num <- function(gene,data_std,method,num){
    cor_list <- vapply(X =rownames(data_std),FUN = function(name0)cor(x = unlist(data_std[gene,]),y = unlist(data_std[name0,]),method = method),0)
    cor_list <- cor_list[order(abs(cor_list),decreasing = TRUE)]
    cor_data_frame <- data.frame(source=rep(x = gene,num),target=names(cor_list)[2:(num+1)],weight=cor_list[2:(num+1)])
  }
  get_gene_cor_thresh <- function(gene,data_std,method,cor_threshold){
    cor_list <- vapply(X =rownames(data_std),FUN = function(name0)cor(x = unlist(data_std[gene,]),y = unlist(data_std[name0,]),method = method),0)
    target <- as.character(names(cor_list)[cor_list>cor_threshold])
    weight <- cor_list[cor_list>cor_threshold]
    cor_data_frame <- data.frame(source=rep(x = gene,length(target)),target=target,weight=weight)
  }
  if (is.null(cor_threshold)){
    edges_coeff_list <- lapply(X = gene_list,FUN = get_gene_cor_num, data_std,method,num)
  }else if(is.numeric(cor_threshold)){
    edges_coeff_list <- lapply(X = gene_list,FUN = get_gene_cor_thresh, data_std,method,cor_threshold)
  }else{
    stop("cor_threshold must be a numeric!")
  }
  edges_coeff <- do.call(rbind,edges_coeff_list)
  edges_coeff2 <- data.frame(source=edges_coeff$target,target=edges_coeff$source,weight=edges_coeff$weight)

  delete_row <- lapply(X = seq_len(nrow(edges_coeff)),
                       FUN = function(x)seq_len(nrow(edges_coeff))[edges_coeff[x,1]==edges_coeff2[,1] & edges_coeff[x,2]==edges_coeff2[,2]])
  delete_row <- lapply(X = seq_len(length(delete_row)),function(x)delete_row[[x]]>x)
  names(delete_row) <- seq_len(nrow(edges_coeff))
  delete_row <- unlist(delete_row)
  delete_row <- delete_row[delete_row==TRUE]
  delete_row <- as.numeric(names(delete_row))
  if (length(delete_row)!=0){
    edges_coeff <- edges_coeff[-delete_row,]
  }
  rownames(edges_coeff) <- 1:nrow(edges_coeff)
  edges_coeff[["pathway"]] <- rep("uncertain",nrow(edges_coeff))
  edges_coeff[["edge_type"]] <- rep("coexp",nrow(edges_coeff))
  return(edges_coeff)
}



# get_bg_related_kegg
#' @title get_bg_related_kegg
#' @description This function will select all genes in all kegg pathways which are directly connected with the genes in \code{gene_list}
#' @param gene_list, a vector of characters, refers to genes ids
#' @param PFPRefnet, an object of PFPRefnet class, it contains all kegg pathways.
#' @param rm_duplicated, a logical, whether to remove the duplicated kegg edges in different pathways.
#' Defalut is \emph{FALSE}
#' @details, It will return a data.frame which can be translated a graph or network.
#' In the data.frame, \code{source} refers to the genes in \code{gene_list}, \code{target} refers to the directly
#' connected genes in kegg, \code{weight} is 0.5, no real means, \code{pathway}
#' refers to the pathway which the edge emerge and \code{edge_type} is "kegg".
#' Note, if \code{rm_duplicated} is \emph{FALSE}, it may return many duplicated edges,
#' which will be complex when plotting a network. If \code{rm_duplicated} is \emph{TRUE},
#' it will retain the first pathway which contains the duplicated edge.
#' @return the related kegg network.
#' @examples
#' # Load the PFPRefnet of human
#' data(PFPRefnet_hsa)
#' # Load the list of diff genes
#' data(gene_list_hsa)
#' # Get the related kegg network
#' #edges_kegg <- get_bg_related_kegg(gene_list_hsa,PFPRefnet_hsa)
#' @export
get_bg_related_kegg <- function(gene_list,PFPRefnet,rm_duplicated = FALSE){
  if (sum(is.na(as.numeric(gene_list))) > 0)
    stop("You should translate all your gene ids into ENTREZID!")
  kegg_edges <- lapply(X = network(PFPRefnet),FUN = getEdgeList)
  data_tf <- llply(kegg_edges,function(kegg_edge)apply(X=kegg_edge,MARGIN=1, FUN=function(x)((x[1] %in% gene_list) | (x[2] %in% gene_list))))
  kegg_edges1 <- lapply(X = names(kegg_edges), function(x)kegg_edges[[x]][data_tf[[x]],])
  kegg_edges1 <-  lapply(X = seq_len(length(kegg_edges1)), function(x)data.frame(source=kegg_edges1[[x]][,1],
                                                                                 target=kegg_edges1[[x]][,2],
                                                                                 weight=rep(0.5,nrow(kegg_edges1[[x]])),
                                                                                 pathway=rep(names(kegg_edges)[x],nrow(kegg_edges1[[x]])),
                                                                                 edge_type=rep("kegg",nrow(kegg_edges1[[x]]))))
  edges_kegg <- do.call(rbind,kegg_edges1)
  if (rm_duplicated == TRUE){
    edges_kegg <- edges_kegg[!duplicated(edges_kegg[,c("source","target")]),]
  }
  return(edges_kegg)
}



# trans_edges_id
#' @title trans_edges_id
#' @description translate the id name in edges_data
#' @param edges_data, the edges_data to translate, it can be the data.frame got
#' from \code{\link{get_exp_cor_edges}} or \code{\link{get_asso_net}}, or a data.frame
#' contains the same colnames with them.
#' @param from_type, a character,the type of gene ID, "ENSEMBL","GO","SYMBOL" and so on.
#' @param to_type, a character,the type of gene ID, "ENSEMBL","GO","SYMBOL" and so on.
#' @param gene_info_db, a gene
#' @details, Translate the id name in edges_data.
#' Note, the \code{from_type} must be consistent with the genes id type in \code{edges_data}.
#' The \code{gene_info_db} must be consistent with the species in \code{edges_data}
#' @return the id of the edges.
#' @examples
#' # Library the datebase of org.Hs.eg.db
#' library(org.Hs.eg.db)
#' # Load the data of human's PFPRefnet
#' data(PFPRefnet_hsa)
#' # Load the list of gene
#' data(gene_list_hsa)
#' # Get the related kegg network
#' #edges_kegg <- get_bg_related_kegg(gene_list_hsa,PFPRefnet_hsa)
#' # Trans the id of edges
#' #edges_kegg <- trans_edges_id(edges_kegg,gene_info_db=org.Hs.eg.db)
#' @export
trans_edges_id <- function(edges_data,from_type="ENTREZID",to_type="SYMBOL",gene_info_db=NULL){
  source_exp <- bitr(geneID = edges_data$source,fromType = from_type,toType = to_type,OrgDb = gene_info_db)
  target_exp <- bitr(geneID = edges_data$target,fromType = from_type,toType = to_type,OrgDb = gene_info_db)
  colnames(source_exp) <- c("source","source_SYMBOL")
  colnames(target_exp) <- c("target","target_SYMBOL")
  edges_data <- merge(x = edges_data,y=source_exp,by="source",all.x=TRUE)
  edges_data <- merge(x = edges_data,y=target_exp,by="target",all.x=TRUE)
  edges_data <- edges_data[!is.na(edges_data$source_SYMBOL)&!is.na(edges_data$target_SYMBOL),]
  edges_data <- edges_data[c("source_SYMBOL","target_SYMBOL","weight","pathway","edge_type")]
  colnames(edges_data) <- c("source","target","weight","pathway","edge_type")
  return(edges_data)
}




# get_asso_net
#' @title merge the edges_coexp and edges_kegg
#' @description This function will remove the co-expressed edges in edges_coexp which also emerge in edges_kegg.
#' @param edges_coexp, a data.frame whose colnames is "source","target","weight","pathway","edge_type".
#' @param edges_kegg, a data.frame whose colnames is "source","target","weight","pathway","edge_type".
#' @param if_symbol, a logical,whether to translate the gene id type. Default is TRUE.
#' @param trans_fun, a function, when \code{if_symbol} is \emph{TRUE},it will use the \code{trans_fun} function
#' to translate the gene ids. Default is \code{trans_edges_id}.
#' @param from_type, a character,the parameter in \code{trans_fun}. It is the type of gene ID, "ENSEMBL","GO","SYMBOL" and so on.
#' @param to_type, a character,the parameter in \code{trans_fun}. It is the type of gene ID, "ENSEMBL","GO","SYMBOL" and so on.
#' @param gene_info_db, an AnnotationDb-object for gene annotation, such as "org.Hs.eg.db".
#' @details, This function will remove the co-expressed edges in edges_coexp which also emerge in edges_kegg.
#' It will return a list contains two data.frames. One is the merged data. Another is
#' the nodes information of the edges.
#' @return the nodes information of the edges.
#' @examples
#' # Load the depends
#' library(org.Hs.eg.db)
#' # Load the data
#' data(PFPRefnet_hsa)
#' data(gene_list_hsa)
#' PFP_s10 <- calc_PFP_score(genes = gene_list_hsa,PFPRefnet = PFPRefnet_hsa,coeff1 = 1,coeff2 = 0.1)
#' rank1 <- rank_PFP(object = PFP_s10,total_rank = TRUE)
#' #pathway_select <- rank1@ref_net_info[1,"id"]
#' #gene_test <- rank1@pathways_score$genes_score[[pathway_select]]$ENTREZID
#' # Get the correlation of edges
#' #edges_coexp <- get_exp_cor_edges(gene_test,data_std)
#' # Get the related kegg
#' #edges_kegg <- get_bg_related_kegg(gene_list_hsa,PFPRefnet_hsa)
#' # Get the related kegg
#' #asso_net <- get_asso_net(edges_coexp,edges_kegg)
#' @export
get_asso_net <- function(edges_coexp,edges_kegg,file_dir=NULL,if_symbol=TRUE,trans_fun = trans_edges_id,from_type="ENTREZID",to_type="SYMBOL",gene_info_db=NULL){
    # trans_id
    if (if_symbol==TRUE){
        if (is.null(gene_info_db)){
            stop("Please input a genome wide annotation packade, for example: org.Hs.eg.db.")
        }else{
          edges_kegg <- trans_edges_id(edges_data = edges_kegg,from_type = from_type,to_type = to_type,gene_info_db = gene_info_db)
          edges_coexp <- trans_edges_id(edges_data = edges_coexp,from_type = from_type,to_type = to_type,gene_info_db = gene_info_db)
        }
    }

    # remove duplicated and combine
    genes_select <- unique(edges_coexp$source)
    delete_row <- vapply(X = seq_len(nrow(edges_coexp)),
                         FUN = function(x)(sum((edges_coexp[x,1]==edges_kegg[,1])&(edges_coexp[x,2]==edges_kegg[,2])|(edges_coexp[x,1]==edges_kegg[,2])&(edges_coexp[x,2]==edges_kegg[,1]))>0)*x,
                         0)
    edges_coexp_new <- edges_coexp[delete_row==0,]
    genes_coexp_new <- unique(c(edges_coexp_new$source,edges_coexp_new$target))
    genes_coexp <- setdiff(genes_coexp_new,genes_select)
    asso_net_edges <- rbind(edges_coexp_new,edges_kegg)
    genes_kegg <- setdiff(unique(c(asso_net_edges$source,asso_net_edges$target)),union(genes_coexp,genes_select))

    asso_net <- list(nodes=data.frame(nodes=c(genes_select,genes_coexp,genes_kegg), group=c(rep("select",length(genes_select)),rep("coexp",length(genes_coexp)),rep("kegg",length(genes_kegg)))),edges=asso_net_edges)

    if (!is.null(file_dir)){
        if (substr(file_dir,nchar(file_dir),nchar(file_dir)) == "/"){
            file_dir <- substr(file_dir,1,(nchar(file_dir)-1))
        }
        write.csv(x = asso_net_edges,file = paste0(file_dir,"/","asso_net_edges.csv"))
        write.csv(x = asso_net[["nodes"]],file=paste0(file_dir,"/","asso_net_nodes.csv"))
    }
    # save(list = c("asso_net"),file = "/home/zx/文档/test/asso_net.RData")
    return(asso_net)
}


#
# library("ggplot2")
#
# load("/home/zx/文档/kangqichuang/PFP/data/PFPRefnet_hsa.RData")
# load("/home/zx/文档/kangqichuang/PFP/data/gene_list_hsa.RData")
# load("/home/zx/文档/PFP/RData/data_std.RData")
# PFP_s10 <- calc_PFP_score(genes = gene_list,PFPRefnet = PFPRefnet_hsa,coeff1 = 1,coeff2 = 0.1)
# rank1 <- rank_PFP(object = PFP_s10,total_rank = TRUE)
# pathway_select <- rank1@ref_net_info[1,"id"]
# gene_test <- rank1@pathways_score$genes_score[[pathway_select]]$ENTREZID
# edges_coexp <- get_exp_cor_edges(gene_test,data_std)
# gene_list2 <- unique(c(edges_coexp$source,edges_coexp$target))
# edges_kegg <- get_bg_related_kegg(gene_list2,PFPRefnet=PFPRefnet_hsa,rm_duplicated = TRUE)
#
# net_test <- get_asso_net(edges_coexp = edges_coexp,edges_kegg = edges_kegg,if_symbol = T,gene_info_db = org.Hs.eg.db)
#
# write.csv(x = net_test$nodes,"/home/zx/文档/nodes.csv")
# write.csv(x = net_test$edges,"/home/zx/文档/edges2.csv")

