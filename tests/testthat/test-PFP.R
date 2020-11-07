data("PFP_test1")
data("PFP_test2")

test_that('method pathways_score of PFP are as expected',{
  pathways.score <- pathways_score(PFP_test1)
  expect_that(pathways.score, is_a('list'))
  expect_that(length(pathways.score), is_equivalent_to(4))
})

test_that('method refnet_info of PFP are as expected',{
  refnet.info <- refnet_info(PFP_test1)
  expect_that(refnet.info, is_a('data.frame'))
  expect_that(ncol(refnet.info), is_equivalent_to(5))
})

test_that('method PFP_score of PFP are as expected',{
  PFP.score <- PFP_score(PFP_test1)
  expect_that(PFP.score, is_a('numeric'))
})

test_that('method stats_test of PFP are as expected',{
  stats.test1 <- stats_test(PFP_test1)
  stats.test2 <- stats_test(PFP_test2)
  expect_that(stats.test1, is_a('data.frame'))
  expect_that(stats.test2, is_a('data.frame'))
})

#test_that('method random_score of PFP are as expected',{
#  random.score1 <- random_score(PFP_test1)
#  random.score2 <- random_score(PFP_test2)
#  expect_that(random.score1, is_a('data.frame'))
#  expect_that(random.score2, is_a('data.frame'))
#  expect_that(nrow(random.score2), is_equivalent_to(0))
#})

test_that('method genes_score of PFP are as expected',{
  genes.score <- genes_score(PFP_test1)
  expect_that(genes.score, is_a('list'))
})

test_that('method sub_PFP of PFP are as expected',{
  sub1_PFP1 <- sub_PFP(object = PFP_test1,group_name = c("Human Diseases","Metabolism"))
  sub2_PFP1 <- sub_PFP(object = PFP_test1,index = rep(list(1:5),6))
  sub3_PFP1 <- sub_PFP(object = PFP_test1,group_name = c("Human Diseases","Metabolism"),
                       index = list(1:5,1:5))
  sub4_PFP1 <- sub_PFP(object = PFP_test1,group_name = c("Human Diseases","Metabolism"),
                       index = list(c("hsa05202","hsa05205"),c("hsa00040","hsa00030")),
                       index_type = "pathway_id")
  sub5_PFP1 <- sub_PFP(object = PFP_test1,group_name = c("Human Diseases","Metabolism"),
                       index = list(c("MicroRNAs in cancer","Chemical carcinogenesis"),
                                    c("Citrate cycle (TCA cycle)","Fructose and mannose metabolism")),
                       index_type = "pathway_name")
  sub5_PFP2 <- sub_PFP(object = PFP_test2,group_name = c("Human Diseases","Metabolism"),
                       index = list(c("MicroRNAs in cancer","Chemical carcinogenesis"),
                                    c("Citrate cycle (TCA cycle)","Fructose and mannose metabolism")),
                       index_type = "pathway_name")
  expect_that(sub1_PFP1, is_a('PFP'))
  expect_that(sub2_PFP1, is_a('PFP'))
  expect_that(sub3_PFP1, is_a('PFP'))
  expect_that(sub4_PFP1, is_a('PFP'))
  expect_that(sub5_PFP1, is_a('PFP'))
  expect_that(sub5_PFP2, is_a('PFP'))
  expect_that(nrow(refnet_info(sub2_PFP1)), is_equivalent_to(30))
  expect_that(length(PFP_score(sub2_PFP1)), is_equivalent_to(30))
})

test_that('method rank_PFP of PFP are as expected',{
  rank1_PFP1 <- rank_PFP(object = PFP_test1,total_rank = TRUE)
  rank2_PFP1 <- rank_PFP(object = PFP_test1,total_rank = FALSE)
  rank1_PFP2 <- rank_PFP(object = PFP_test2,total_rank = TRUE)
  rank2_PFP2 <- rank_PFP(object = PFP_test2,total_rank = FALSE)
  expect_that(rank1_PFP1,is_a('PFP'))
  expect_that(rank2_PFP1,is_a('PFP'))
  expect_that(rank1_PFP2,is_a('PFP'))
  expect_that(rank2_PFP2,is_a('PFP'))
})
