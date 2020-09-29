
library("clusterProfiler") # translate id
library("ggplot2")
library("graph")


# translate gene id to ENTREZID
trans_id <- function(gene_list,gene_id_type,gene_info_db){
  gene_t <- gene_list
  if (gene_id_type != "ENTREZID"){
    gene_n <- bitr(gene_t,fromType = gene_id_type,
                   toType = "ENTREZID",#"ENSEMBL"
                   OrgDb = gene_info_db,drop = F)
    gene_n <- gene_n[!duplicated(gene_n$ENSEMBL),]
  }else{
    gene_n <- data.frame(ENTREZID=gene_t)
  }
  return (gene_n)
}


get_PFP_score <- function(genes,PFPRefnet,coeff1=1,coeff2=0.1,statistic = FALSE,bg_genelist=NULL,num_loop=10){
  
  # get graph score
  get_graph_score <- function(graph,genes,coeff1,coeff2){

    kegg_edges0 <- strsplit(x = names(edgeData(graph0)),split = "\\|")
    kegg_edges1 <- data.frame(t(data.frame(kegg_edges0)))
    kegg_nodes0 <- nodes(graph0)
    con_rate <- length(kegg_edges0)/length(kegg_nodes0)
    genes_inter <- as.vector(intersect(genes,kegg_nodes0))
    graph_score <- data.frame(ENTREZID=genes_inter,score0=rep(1,length(genes_inter)))
    
    # score1
    score1_tmp <- vapply(X = kegg_edges0,FUN = function(x)sum(match(x = x,table = genes_inter,nomatch = 0) == 0)==0,FUN.VALUE = TRUE)
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
    score2 <- vapply(X = seq_len(length(genes_inter)),function(x)sum(score2_tmp[,x])+sum(score2_tmp[x,])-score2_tmp[x,x],0)
    score2 <- data.frame(ENTREZID=genes_inter,score2=score2)
    score2[["score2"]] <- score2[,"score2"]*coeff2
    
    graph_score <- merge(x = graph_score, y = score1, by = "ENTREZID", all.x = TRUE)
    graph_score <- merge(x = graph_score, y = score2, by = "ENTREZID", all.x = TRUE)
    
    graph_score[["score"]] <- rowSums(graph_score[,c("score1","score2")],na.rm = T)/con_rate+graph_score[,c("score0")]
    graph_score[is.na(graph_score)] <- 0
    return(graph_score)
  }
  
  
  # test p.value
  test_pavlue <- function(num_genes,bg_genelist,PFP_score,num_loop,PFPRefnet,coeff1,coeff2,get_graph_score){
    
    genes_list <- replicate(num_loop, sample(x = bg_genelist,size = num_genes,replace = F))
    genes_list <- split(genes_list,col(genes_list))
    
    test_once <- function(genes,PFPRefnet,coeff1,coeff2){
      genes_score_t <- lapply(network(PFPRefnet),get_graph_score,genes=genes,coeff1=coeff1,coeff2=coeff2)
      PFP_score_t <- vapply(genes_score_t,function(x)sum(x$score),0) 
    }
    
    test_total <- data.frame(t(data.frame(lapply(X = genes_list,test_once,PFPRefnet,coeff1,coeff2))))
    # test normal distribution
    test_once_pvalue <- function(name,mu_s,vec_s){
      if (length(unique(vec_s[,name]))==1){
        tt = 1
      }else{
        tt <- shapiro.test(vec_s[,name])$p.value
      }
      if (tt < 0.05){
        # student t test
        pvalue0 <- t.test(x = vec_s[,name], mu=mu_s[name])$p.value
      }else{
        # wilcox.test
        pvalue0 <- wilcox.test(x = vec_s[,name], mu = mu_s[name])$p.value
        if (is.na(pvalue0)){
          pvalue0 <- 1
        }
      }
      return(pvalue0)
    }
    p_values_t <- vapply(X = names(PFP_score),test_once_pvalue,PFP_score,test_total,FUN.VALUE = 0)
    return(list(random_score=test_total,p_value = p_values_t))
  }
  
  
  if (sum(is.na(as.numeric(genes))) > 0)
    stop("You should translate all your gene ids into ENTREZID!")
  
  genes_score <- lapply(network(PFPRefnet),get_graph_score,genes=genes,coeff1=coeff1,coeff2=coeff2)
  PFP_score <- vapply(genes_score,function(x)sum(x$score),0)
  
  
  if (statistic==TRUE){
    if(is.null(bg_genelist)){
      bg_genelist <- unique(unlist(lapply(X = network(PFPRefnet),nodes)))
      num_genes <- length(unique(intersect(genes,bg_genelist)))
    }else{
      if (sum(is.na(as.numeric(bg_genelist))) > 0){
        stop("You should translate all your background genes ids into ENTREZID!")
      num_genes <- length(genes)
      }
    }
    random_tests <- test_pavlue(num_genes,bg_genelist,PFP_score,num_loop,PFPRefnet,coeff1,coeff2,get_graph_score)
    random_score <- random_tests[["random_score"]]
    p_value <- random_tests[["p_value"]]
  }else{
    random_score <- data.frame()
    p_value <- numeric()
  }
  
  PFP_res <- new(Class = "PFP",
                 pathways_score = list(PFP_score=PFP_score,p_value=p_value,random_score=random_score,genes_score=genes_score),
                 net_info = net_info(PFPRefnet))
  
}



# test
gene_list0 = as.vector(read.csv(file = paste0("/home/zx/文档/drug_database/result/",1,"/diff/diff_sig.csv"),header = T)$gene)
gene_id_type0 = "ENSEMBL"
con_coeff10 = 1 
con_coeff20 = 0.1 

library("org.Mm.eg.db")
load(file = "~/文档/PFP/RData/PFPRefnet.RData")
gene_n <- trans_id(gene_list = gene_list0,gene_id_type = gene_id_type0,gene_info_db = org.Mm.eg.db)
genes <- unique(gene_n$ENTREZID[!is.na(gene_n$ENTREZID)])
res <- get_PFP_score(genes,PFPRefnet,coeff1=con_coeff10,coeff2=con_coeff20,statistic = TRUE)
