skip_if(as.logical(Sys.getenv("SKIP_NP_SERVER_TESTS")))

test_that("hemibrain_meshes works", {
  expect_is(cal <- hemibrainr::hemibrain_roi_meshes('CA(R)', microns = TRUE), 'neuronlist')
  expect_is(ca <- cal[[1]], 'mesh3d')
  # less than 100 microns across
  expect_lt(max(diff(nat::boundingbox(ca))), 100)
})
