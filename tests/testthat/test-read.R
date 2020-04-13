skip_if(as.logical(Sys.getenv("SKIP_NP_SERVER_TESTS")))

test_that("hemibrain_read_neurons works", {
  kcg2ids=neuprint_search('type:KCg')[1:2,]
  expect_is(kcg2 <- hemibrain_read_neurons(kcg2ids, savedir = F, local = F, clean = T),
            'neuronlist')
  expect_false(any(duplicated(colnames(kcg2))))

})
