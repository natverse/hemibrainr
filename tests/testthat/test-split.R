skip_if(as.logical(Sys.getenv("SKIP_NP_SERVER_TESTS")))

test_that("hemibrain_skeleton_check works", {
  ids = c("203253072")

  neurons = neuprint_read_neurons(ids)

  expect_is(neurons.flow <- hemibrain_flow_centrality(neurons),
            'neuronlist')
  expect_is(neurons.flow[[1]], 'neuron')
  expect_length(neurons.flow[[1]][['primary.branch.point']], 1L)
})
