data(PFP_test)

## PFP object
test_that('class of PFP and its components are as expected',{
  expect_that(PFP_test, is_a('PFP'))
  expect_that(PFP_test@pathways_score, is_a('list'))
  expect_that(PFP_test@ref_net_info, is_a('data.frame'))
})

test_that('method refnet_names of PFP are as expected',{
  refnet.names <- refnet_names(PFP_test)
  expect_that(refnet.names, is_a('character'))
  expect_that(length(refnet.names), is_equivalent_to(300))

})

test_that('method p value of PFP are as expected',{
  p_value.PFP <- p_value(PFP_test)
  expect_that(p_value.PFP, is_a('numeric'))
  expect_that(length(p_value.PFP), is_identical_to(300L))

})

test_that('genes scores of PFP are as expected',{
  genes_score.PFP <- genes_score(PFP_test)
  expect_that(genes_score.PFP, is_a('list'))
  expect_that(length(genes_score.PFP), is_identical_to(300L))
  expect_that(length(genes_score.PFP[1]), is_identical_to(1L))

})

test_that('subnet of PFP are as expected',{
  group <- unique(ref_net_info(PFP_test)$group)[c(2,4)]
  sub_PFP<- sub_PFP(PFP_test,group)
  expect_that(sub_PFP,is_a('PFP'))
  expect_that(length(sub_PFP),is_equivalent_to(1))
})
