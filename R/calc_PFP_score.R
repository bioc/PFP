#' @title Get the pathway fingerprint of a gene_list
#' @description It can evaluate the performance of a gene list in the pathway networks.
#' @param genes, a vector of characters
#' @param PFPRefnet, A PFPRefnet class
#' @param coeff1, a numeric, the coefficient 1 in PFP model
#' @param coeff2, a numeric, the coefficient 2 in PFP model
#' @param statistic, a logical,whether to do the statistical test
#' @param bg_genelist, a vector of characters, background gene set for the statistical test
#' @param adjust_method, statistic test method for adjust the p_value.
#' It could be "holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none".
#' @details The main part of pathway fingerprint. PFP is used to evaluate the performance of
#' a gene_list in some pathway networks by considering the genes' topological location in
#' a pathway. Then we can get every gene's score and the pathway score is caculated by sum
#' all genes' score. All pathways' scores combine the pathway fingerprint.
#' @return The main part of pathway fingerprint.
#' @examples
#' # Get the example data.
#' data(gene_list_hsa)
#' data("PFPRefnet_hsa")
#' # Calculate the PFP score.
#' PFP <- calc_PFP_score(gene_list_hsa,PFPRefnet)
#' @export
calc_PFP_score <- function(genes,PFPRefnet,coeff1=1,coeff2=0.1,statistic = FALSE,
                           bg_genelist=NULL,adjust_method="BH"){
  # get graph score
  get_graph_score <- function(graph0,genes,coeff1,coeff2){

    kegg_edges0 <- strsplit(x = names(edgeData(graph0)),split = "\\|")
    kegg_edges1 <- data.frame(t(data.frame(kegg_edges0)))
    kegg_nodes0 <- nodes(graph0)
    con_rate <- length(kegg_edges0)/length(kegg_nodes0)
    genes_inter <- as.vector(intersect(genes,kegg_nodes0))
    graph_score <- data.frame(ENTREZID=genes_inter,
                              score0=rep(1,length(genes_inter)))

    # score1
    score1_tmp <- vapply(X = kegg_edges0,
                         FUN = function(x)sum(match(x = x,
                                                    table = genes_inter,
                                                    nomatch = 0) == 0)==0,
                         FUN.VALUE = TRUE)
    if (sum(score1_tmp) == 0){
      score1 <- data.frame(matrix(nrow = 0,ncol=2))
      colnames(score1) <- c("ENTREZID","score1")
    }else{
      score1 <- data.frame(table(unlist(kegg_edges0[score1_tmp])))
      colnames(score1) <- c("ENTREZID","score1")
      score1[["score1"]] <- score1[,"score1"]*coeff1
    }

    # score2
    left_con <- lapply(genes_inter,function(x)kegg_edges1[kegg_edges1[,1]==x,2])
    names(left_con) <- genes_inter
    right_con <- lapply(genes_inter,function(x)kegg_edges1[kegg_edges1[,2]==x,1])
    names(right_con) <- genes_inter
    fun_intersect <- function(set,set_list){
      vapply(X = set_list, FUN = function(x)(length(intersect(x,set))),0)
    }
    score2_tmp <- data.frame(lapply(X = left_con,FUN = fun_intersect,right_con))
    score2 <- vapply(X = seq_len(length(genes_inter)),
                     function(x)sum(score2_tmp[,x])+sum(score2_tmp[x,])-
                       score2_tmp[x,x],0)
    score2 <- data.frame(ENTREZID=genes_inter,score2=score2)
    score2[["score2"]] <- score2[,"score2"]*coeff2

    graph_score <- merge(x = graph_score, y = score1, by = "ENTREZID",
                         all.x = TRUE)
    graph_score <- merge(x = graph_score, y = score2, by = "ENTREZID",
                         all.x = TRUE)
    graph_score[["score"]] <- rowSums(graph_score[,c("score1","score2")],
                                      na.rm = TRUE)/con_rate+graph_score[,c("score0")]
    graph_score[is.na(graph_score)] <- 0
    return(graph_score)
  }

  test_pavlue <- function(genes,bg_genelist,PFPRefnet,adjust_method){
    genes <- intersect(genes,bg_genelist)
    sep_nodes <- lapply(X = names(network(PFPRefnet)),FUN = function(x)intersect(nodes(network(PFPRefnet)[[x]]),bg_genelist))
    names(sep_nodes) <- names(network(PFPRefnet))
    sep_nodes_hit <- lapply(X = names(sep_nodes),FUN = function(x)intersect(sep_nodes[[x]],genes))
    names(sep_nodes_hit) <- names(network(PFPRefnet))
    p_values <- vapply(X = names(sep_nodes),FUN = function(x)phyper(q = length(sep_nodes_hit[[x]])-1,
                                                                    m = length(sep_nodes[[x]]),
                                                                    n = length(bg_genelist)-length(sep_nodes[[x]]),
                                                                    k = length(genes),
                                                                    lower.tail = FALSE),FUN.VALUE = 1)
    p_adjust <- p.adjust(p = p_values, method = adjust_method)
    return(list(p_value = p_values,p_adj_value = p_adjust))
  }


  if (sum(is.na(as.numeric(genes))) > 0)
    stop("You should translate all your gene ids into ENTREZID!")

  genes_score <- lapply(network(PFPRefnet),get_graph_score,genes=genes,
                        coeff1=coeff1,coeff2=coeff2)
  PFP_score <- vapply(genes_score,function(x)sum(x$score),0)

  if (statistic==TRUE){
    if(is.null(bg_genelist)){
      bg_genelist <- unique(unlist(lapply(X = network(PFPRefnet),nodes)))
    }else{
      if (sum(is.na(as.numeric(bg_genelist))) > 0){
        stop("You should translate all your background genes ids into ENTREZID!")
      }
    }
    random_tests <- test_pavlue(genes,bg_genelist,PFPRefnet,adjust_method = adjust_method)
    p_value <- random_tests[["p_value"]]
    p_adj_value <- random_tests[["p_adj_value"]]
  }else{
    p_value <- numeric()
    p_adj_value <- numeric()
  }

  PFP_res <- new(Class = "PFP",
                 pathways_score = list(PFP_score=PFP_score,
                                       stats_test=data.frame(p_value=p_value,p_adj_value=p_adj_value),
                                       genes_score=genes_score),
                 refnet_info = net_info(PFPRefnet))
  return(PFP_res)
}


