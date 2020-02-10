library(testthat)
library(hemibrainr)

conn=try(neuprint_login())
if(isTRUE(conn$server=="https://neuprint.janelia.org")) {
  Sys.setenv(SKIP_NP_SERVER_TESTS="FALSE")
} else {
  Sys.setenv(SKIP_NP_SERVER_TESTS="TRUE")
}

test_check("hemibrainr")
