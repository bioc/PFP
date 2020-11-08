data("gene_list_hsa")
data("PFPRefnet_hsa")
data("PFP_test1")
data("data_std")

test_that('function calc_PFP_score of PFP are as expected',{
  PFP_test <- calc_PFP_score(genes = gene_list_hsa,PFPRefnet = PFPRefnet_hsa)
  expect_that(PFP_test, is_a('PFP'))
})

test_that('functions in net.R of PFP are as expected',{
  rank1 <- rank_PFP(object = PFP_test1,total_rank = T)
  pathway_select <- refnet_info(rank1)[1,"id"]
  gene_test <- pathways_score(rank1)$genes_score[[pathway_select]]$ENTREZID
  edges_coexp <- get_exp_cor_edges(gene_test,data_std)
  gene_list2 <- unique(c(edges_coexp$source,edges_coexp$target))
  edges_kegg <- get_bg_related_kegg(gene_list2,PFPRefnet=PFPRefnet_hsa,rm_duplicated = TRUE)
  net_test <- get_asso_net(edges_coexp = edges_coexp,edges_kegg = edges_kegg,if_symbol = T,gene_info_db = org.Hs.eg.db)
  expect_that(net_test, is_a('list'))
  expect_that(net_test$nodes, is_a('data.frame'))
  expect_that(net_test$edges, is_a('data.frame'))
})

test_that('function plot_PFPlist of PFP are as expected',{
  PFPlist <- list(PFP_test1,PFP_test1)
  plot_PFPlist(PFPlist)
})


