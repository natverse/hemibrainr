skip_if(as.logical(Sys.getenv("SKIP_NP_SERVER_TESTS")))

test_that("hemibrain_skeleton_check works", {
  ids = c("5813056323", "579912201", "5813015982", "973765182", "885788485",
          "915451074", "5813032740", "1006854683", "5813013913", "5813020138",
          "853726809", "916828438", "5813078494", "420956527", "486116439",
          "573329873", "5813010494", "5813040095", "514396940", "665747387",
          "793702856", "451644891", "482002701", "391631218", "390948259",
          "390948580", "452677169", "511262901", "422311625", "451987038"
  )

  # Read in these neurons
  expect_is(neurons <- neuprintr::neuprint_read_neurons(ids[1:2]), 'neuronlist')
  expect_is(neurons[[1]], 'neuron')

  # Re-root
  expect_is(neuron.checked <- hemibrain_skeleton_check(neurons, OmitFailures = TRUE),
            'neuronlist')
  expect_is(neuron.checked[[1]], 'neuprintneuron')
})

