### Code to prepare datasets ###

### Split points
hemibrain_splitpoints_polypre_centrifugal_distance <- read.csv("data-raw/hemibrain_all_neurons_splitpoints_polypre_centrifugal_distance.csv")
hemibrain_splitpoints_pre_centrifugal_distance <- read.csv("data-raw/hemibrain_all_neurons_splitpoints_pre_centrifugal_distance.csv")
hemibrain_splitpoints_polypre_centrifugal_presynapses <- read.csv("data-raw/hemibrain_all_neurons_splitpoints_polypre_centrifugal_presynapses.csv")
usethis::use_data(hemibrain_splitpoints_polypre_centrifugal_distance, overwrite = TRUE)
usethis::use_data(hemibrain_splitpoints_polypre_centrifugal_presynapses, overwrite = TRUE)

### Metrics
hemibrain_metrics <- read.csv("data-raw/hemibrain_all_neurons_metrics.csv")
rownames(hemibrain_metrics) <- hemibrain_metrics$bodyid
usethis::use_data(hemibrain_metrics, overwrite = TRUE)

### Consensus

### Mesh
hemibrain.surf = readobj::read.obj("data-raw/hemibrain_raw.obj", convert.rgl = TRUE)
hemibrain.surf = nat::as.hxsurf(hemibrain.surf[[1]])
hemibrain_microns.surf = hemibrain.surf*(8/1000)
nat.templatebrains::regtemplate(hemibrain_microns.surf) = "JRCFIB2018F"
nat.templatebrains::regtemplate(hemibrain.surf) = "JRCFIB2018Fraw"
usethis::use_data(hemibrain_microns.surf, overwrite = TRUE)
usethis::use_data(hemibrain.surf, overwrite = TRUE)
