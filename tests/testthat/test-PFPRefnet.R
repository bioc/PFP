data("PFPRefnet_hsa")

test_that('method network of PFPRefnet are as expected',{
  refnet.networks <- network(PFPRefnet_hsa)
  expect_that(refnet.networks, is_a('list'))
})

test_that('method net_info of PFPRefnet are as expected',{
  refnet.info <- net_info(PFPRefnet_hsa)
  expect_that(refnet.info, is_a('data.frame'))
})

test_that('method group of PFPRefnet are as expected',{
  refnet.group <- group(PFPRefnet_hsa)
  expect_that(refnet.group, is_a('list'))
  expect_that(length(refnet.group), is_equivalent_to(3))
})

test_that('method net_names of PFPRefnet are as expected',{
  refnet.netnames <- net_names(PFPRefnet_hsa)
  expect_that(refnet.netnames, is_a('character'))
})

test_that('method subnet of PFPRefnet are as expected',{
  group_name <- unique(net_info(PFPRefnet_hsa)$group)[c(2,4)]
  sub_net1 <- subnet(PFPRefnet_hsa, group_name, index = list(1:5,1:5))
  sub_net2 <- subnet(PFPRefnet_hsa,index = rep(list(1:5),6))
  sub_net3 <- subnet(object = PFPRefnet_hsa,group_name = c("Human Diseases","Metabolism"),
                       index = list(c("hsa05202","hsa05205"),c("hsa00040","hsa00030")),
                       index_type = "pathway_id")
  sub_net4 <- subnet(object = PFPRefnet_hsa,group_name = c("Human Diseases","Metabolism"),
                       index = list(c("MicroRNAs in cancer","Chemical carcinogenesis"),
                                    c("Citrate cycle (TCA cycle)","Fructose and mannose metabolism")),
                       index_type = "pathway_name")
  expect_that(sub_net1,is_a('PFPRefnet'))
  expect_that(sub_net2,is_a('PFPRefnet'))
  expect_that(sub_net3,is_a('PFPRefnet'))
  expect_that(sub_net4,is_a('PFPRefnet'))
  expect_that(nrow(net_info(sub_net2)), is_equivalent_to(30))
  expect_that(nrow(net_info(sub_net3)), is_equivalent_to(4))
})
