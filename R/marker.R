setwd("/Users/xiaochang/Desktop/drug_database") # 设置文件夹

# library("graphite") # get pathway
library("clusterProfiler") # 基因ID转换
library("ggplot2")
library("xlsx")



# 参数
parameters <- commandArgs(T)

print(parameters[1])
print(parameters[2])
print(parameters[3])
print(parameters[4])
print(parameters[5])
print(parameters[6])


# gene_list0 = unlist(strsplit(x = parameters[1] ,split = "\\|")) # "ENSMUSG00000076492|ENSMUSG00000076604|ENSMUSG00000093346|ENSMUSG00000091159|ENSMUSG00000065612" #parameter1 基因列表
# gene_list0 = parameters[1] # "ENSMUSG00000076492|ENSMUSG00000076604|ENSMUSG00000093346|ENSMUSG00000091159|ENSMUSG00000065612" #parameter1 基因列表
# gene_id_type0 = parameters[2] # "ENSEMBL"/"ENTREZID"/"SYMBOL" #parameter2 基因ID类型
# species0 = parameters[3] # "hsa" "mmu" #parameter3 物种名称（依据kegg对物种的表示方法）
# con_coeff10 = as.numeric(parameters[4]) # 1 #parameter4 通路指纹模型相邻节点得分参数
# con_coeff20 = as.numeric(parameters[5]) # 0.1 #parameter5 通路指纹模型相隔一个节点得分参数
# task_id0 = parameters[6] # 3 # 任务编号,系统生成


gene_list0 = as.vector(read.csv(file = paste0("result/",1,"/diff/diff_sig.csv"),header = T)$gene)
species0 = "mmu" # "hsa" "mmu" #parameter3 物种名称（依据kegg对物种的表示方法）
gene_id_type0 = "ENSEMBL"
con_coeff10 = 1 # 1 #parameter4 通路指纹模型相邻节点得分参数
con_coeff20 = 0.1 # 0.1 #parameter5 通路指纹模型相隔一个节点得分参数
task_id0 = 1
num_loop0 =200


if (!is.na(con_coeff10)){
  con_coeff10 <- 1
}
if (!is.na(con_coeff20)){
  con_coeff20 <- 0.1
}

pathways_id_dir <- "RData/pathways_id.xlsx" # 含有通路组别信息


# 获取kegg通路信息
if (species0=="hsa"){
  load(file = paste0("RData/hsa_kegg_info.RData"))
  library("org.Hs.eg.db")
  gene_info_db <- org.Hs.eg.db
}else if (species0=="mmu"){
  load(file = paste0("RData/mmu_kegg_info.RData")) 
  library("org.Mm.eg.db")
  gene_info_db <- org.Mm.eg.db
}

# 创建结果文件夹
if (!dir.exists(paste0("./result/",task_id0))){
  dir.create(paste0("./result/",task_id0))
}
if (!dir.exists(paste0("./result/",task_id0,"/PFP"))){
  dir.create(paste0("./result/",task_id0,"/PFP"))
}
if (!dir.exists(paste0("./result/",task_id0,"/enrich"))){
  dir.create(paste0("./result/",task_id0,"/enrich"))
}

# 对基因列表进行id 转换
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


get_PFP_score <- function(gene_n,pathways_id_dir,kegg_edges,kegg_nodes,coeff1=1,coeff2=0.1,task_id){
  
  # 通路指纹得分
  get_edges_score <- function(x,genes,kegg_edges,kegg_nodes,coeff1,coeff2){
    con_rate <- nrow(kegg_edges[[x]])/length(kegg_nodes[[x]])
    genes_inter <- as.vector(intersect(genes,kegg_nodes[[x]]))
    edges_score <- data.frame(genes_inter)
    colnames(edges_score) <- c("ENTREZID")
    score1 <- c()
    for (j in genes_inter){
      score1[j] = 0
    }
    for (j in nrow(kegg_edges[[x]])){
      source0 <- as.vector(kegg_edges[[x]][j,1])
      target0 <- as.vector(kegg_edges[[x]][j,2])
      if ((source0 %in% genes_inter) && (target0 %in% genes_inter)){
        score1[source0] <- score1[source0]+coeff1
        score1[target0] <- score1[target0]+coeff1
      }
    } 
    for (gene0 in genes_inter){
      gene1 <- kegg_edges[[x]][kegg_edges[[x]][,1] == gene0,2]
      for (gene10 in gene1){
        gene2 <- kegg_edges[[x]][kegg_edges[[x]][,1] == gene10,2]
        for (gene20 in gene2 ){
          if (gene20 %in% genes_inter){
            score1[gene0] <- score1[gene0]+coeff2
            score1[gene20] <- score1[gene20]+coeff2
          }
        }
      }
    }
    edges_score[["score"]] <- score1/con_rate
    return(edges_score)
  }
  
  # PFP绘图
  get_pfp_single <- function(result_PFP,task_id){
    title <- "PFP_score"
    save_root <- paste0("result/",task_id,"/PFP")
    sim_df <- data.frame(refnet_index=1:length(result_PFP$id),PFP_score= as.vector(result_PFP$PFP_score),group = as.vector(result_PFP$group))
    
    pic <- ggplot(sim_df,aes(x = refnet_index, y = PFP_score)) +
      geom_point(size = 1, aes(color = group)) +
      geom_segment(aes(xend = refnet_index, yend = 0, color = group),size = 0.5) +
      labs(title = title)+
      theme(plot.title = element_text(hjust = 0.5))
    ggsave(filename = paste0(save_root,"/","plot_PFP.pdf"),plot = pic,device = "pdf",width = 12,height = 4)
  }
  
  # 统计检验
  test_pavlue <- function(PFP_score,x,genelist,kegg_edges,kegg_nodes,coeff1,coeff2,gene_id_type,gene_info_db,trans_id,get_edges_score,num_loop){
    verify  <- c()
    bg_genelist <- keys(x = gene_info_db,keytype = "ENSEMBL")
    for(i in 1:num_loop){
      genelist1 <- sample(x = bg_genelist,size = length(gene_list0),replace = F)
      gene_n <- trans_id(gene_list = genelist1,gene_id_type = gene_id_type,gene_info_db = gene_info_db)
      gene <- unique(gene_n[!is.na(gene_n$ENTREZID),"ENTREZID"])
      res <-get_edges_score(x = x,genes = gene,kegg_edges,kegg_nodes,coeff1,coeff2)
      res0 <- sum(res$score)+nrow(res)
      verify <- c(verify,res0)
    }
    # 检验数据是否服从正态分布
    tt <- shapiro.test(verify)$p.value
    if (tt < 0.05){
      # 参数检验
      pvalue0 <- t.test(x = verify,mu=PFP_score)$p.value
    }else{
      # 非参数检验
      pvalue0 <- wilcox.test(x = verify,mu = PFP_score)$p.value
    }
    return(pvalue0)
  }
  
  
  # id translate get ENTREZID
  gene <- unique(gene_n[!is.na(gene_n$ENTREZID),"ENTREZID"])
  
  # get PFP_score
  kegg_pathway_scores <- c()
  kegg_pathway_pvalues <- c()
  diff_sig_PFP <- gene_n
  for (x in names(kegg_edges)){
    kegg_edges_score <- get_edges_score(x=x,genes = gene,kegg_edges = kegg_edges,kegg_nodes = kegg_nodes,coeff1 = coeff1,coeff2 = coeff2)
    PFP_score0 <- sum(kegg_edges_score$score)+length(kegg_edges_score$score)
    test_pavlue0 <- test_pavlue(PFP_score = PFP_score0, x=x,genelist = gene_list0,kegg_edges=kegg_edges,kegg_nodes=kegg_nodes,coeff1=con_coeff10,coeff2=con_coeff20,gene_id_type=gene_id_type0,gene_info_db=gene_info_db,trans_id,get_edges_score,num_loop=num_loop0)
    kegg_pathway_scores <- c(kegg_pathway_scores,PFP_score0)
    kegg_pathway_pvalues <- c(kegg_pathway_pvalues,test_pavlue0)
    diff_sig_PFP <- merge(diff_sig_PFP,kegg_edges_score,by = "ENTREZID",all.x = T)
  }
  colnames(diff_sig_PFP) <- c("ENTREZID","ENSEMBL",names(kegg_edges))
  write.csv(x = diff_sig_PFP,file = paste0("result/",task_id,"/PFP/diff_sig_PFP.csv"),row.names = F)
  
  #pathway_net.csv
  pathway_ids <- substr(x = names(kegg_edges),start = 4,stop = 8)
  result_PFP <- data.frame(id=pathway_ids,PFP_score=kegg_pathway_scores,pvalue = kegg_pathway_pvalues)
  group_name_index <- read.xlsx(file = pathways_id_dir,sheetIndex = 1,header = T)
  result_PFP <- merge(x = result_PFP,y = group_name_index,by="id",all.x=T)
  result_PFP <- result_PFP[order(result_PFP[,"group"],decreasing = F),]
  rownames(result_PFP) <- 1:nrow(result_PFP)
  write.xlsx(x = result_PFP[,c("id","name","PFP_score","pvalue","group")],file = paste0("result/",task_id0,"/PFP/PFP_score.xlsx"))
  # plot PFP
  get_pfp_single(result_PFP = result_PFP,task_id = task_id)
}


# KEGG富集分析
kegg_enrich <- function(gene_n,spec,gene_info_db,task_id){
  
  gene <- gene_n$ENTREZID[!is.na(as.vector(gene_n$ENTREZID))]
  bg_genes <- keys(x = gene_info_db,keytype  = "ENTREZID")
  # KEGG富集分析
  
  kk <- enrichKEGG(gene=gene, universe=bg_genes, organism=spec)
  write.csv(kk, paste0("result/",task_id0,"/enrich/kegg.enrich.csv"), quote=F, row.names=F)
  
  # KEGG条形图
  pdf(paste0("result/",task_id,"/enrich/", "kegg.bar.pdf"))#, width=5000, height=5000, res=300)
  p <- barplot(kk, title="Enrichment KEGG")
  print(p)
  dev.off()
  # # KEGG气泡图
  # pdf(paste0("result/",task_id,"/enrich/", "kegg.dot.pdf"))#, width=5000, height=5000, res=300)
  # p2 <- dotplot(kk, title="Enrichment KEGG")
  # print(p2)
  # dev.off()
  
}


# 函数执行流程
gene_n <- trans_id(gene_list = gene_list0,gene_id_type = gene_id_type0,gene_info_db = gene_info_db)
kegg_enrich(gene_n = gene_n,spec = species0,gene_info_db = gene_info_db,task_id = task_id0)
get_PFP_score(gene_n,pathways_id_dir,kegg_edges,kegg_nodes,coeff1=con_coeff10,coeff2=con_coeff20,task_id = task_id0)
