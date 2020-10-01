data(PFPRefnet)
PFPRefnet <- new(Class = "PFPRefnet",network = network(PFPRefnet),
                 net_info = net_info(PFPRefnet))

test_that('class of PFPRefnet and its components are as expected',{
  expect_that(PFPRefnet, is_a('PFPRefnet'))
  expect_that(PFPRefnet@network, is_a('list'))
  expect_that(PFPRefnet@net_info, is_a('data.frame'))
})

test_that('method group of PFPRefnet are as expected',{
  refnet.group <- group(PFPRefnet)
  expect_that(refnet.group, is_a('list'))
  expect_that(length(refnet.group), is_equivalent_to(3))

})

test_that('method netnames of PFPRefnet are as expected',{
  refnet.netnames <- net_names(PFPRefnet)
  expect_that(refnet.netnames, is_a('character'))
  expect_that(length(refnet.netnames), is_equivalent_to(300))
})

test_that('subnet of PFPRefnet are as expected',{
  group_name <- unique(net_info(PFPRefnet)$group)[c(2,4)]
  sub_PFPRefnet <- subnet(PFPRefnet, group_name, index = list(1:5,1:5))
  expect_that(sub_PFPRefnet,is_a('PFPRefnet'))
  expect_that(length(sub_PFPRefnet),is_equivalent_to(1))
})
