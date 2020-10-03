skip_if(as.logical(Sys.getenv("SKIP_NP_SERVER_TESTS")))

test_that("hemibrain_clean_skeleton works", {
  ids = c("203253072")
  neurons = neuprint_read_neurons(ids)
  neurons.flow = hemibrain_flow_centrality(neurons)
  neurons.cleaned = hemibrain_clean_skeleton(neurons.flow)

  expect_is(neurons.cleaned,
            'neuronlist')
  expect_is(neurons.cleaned[[1]], 'neuron')
  expect_length(neurons.cleaned[[1]][['primary.branch.point']], 1L)
})
