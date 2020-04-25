skip_if(as.logical(Sys.getenv("SKIP_NP_SERVER_TESTS")))

test_that("grouped_adjacency works", {

  expect_is(da2pnkc <- grouped_adjacency("/.*DA2.*PN.*", 'KC', ingroup = NULL),
            'matrix')
  expect_true(any(grepl("KC", colnames(da2pnkc))))
})

