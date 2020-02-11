
test_that("ids", {
  fakeids=c(123, 456, 789)
  expect_equal(neuprint_ids(fakeids), as.character(fakeids))
  expect_error(neuprint_ids(-1))

})

skip_if(as.logical(Sys.getenv("SKIP_NP_SERVER_TESTS")))

test_that("ids w search", {
  expect_error(neuprint_ids("-1"))
  expect_is(neuprint_ids("MBON"), "character")
})
