## code to prepare datasets
hemibrain_splitpoints_distance <- read.csv("data-raw/hemibrain_all_neurons_splitpoints.csv")
hemibrain_metrics <- read.csv("data-raw/hemibrain_all_neurons_metrics.csv")
rownames(hemibrain_metrics) <- hemibrain_metrics$bodyid
usethis::use_data(hemibrain_splitpoints_distance, overwrite = TRUE)
usethis::use_data(hemibrain_metrics, overwrite = TRUE)
