###### Script to split all neurons in the hemibrain into putative axon / dendrite #######
######## Originally run from JData5, the Jefferis lab drive on Hex at the MRC LMB #######

# source /public/gcc/gcc7_1_0.csh
# In order to then update Rvcg in R if needed or remotes::install_github('natverse/nat.h5reg')
# setenv PATH ${PATH}:/public/flybrain/hdf5/lib/
# setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:/public/flybrain/hdf5/lib/
# At the LMB, this code must be run on Hex/Hal not Max
# setenv LD_LIBRARY_PATH /public/flybrain/R/R-3.5.1/lib64/R/library/rJava/libs/
# source /public/flybrain/flybrain.csh

# # Transform example for debugging
# library(nat.jrcbrains)
# library(rJava)
# pts=structure(c(20107.5, 20075.3, 20114.8, 20122.4, 20200.8, 12214.9,
#                 12279.3, 12318.2, 12325, 12394.7, 12058.7, 12198.2, 12244.8,
#                 12252.4, 12330.8), .Dim = c(5L, 3L), .Dimnames = list(NULL, c("X",
#                                                                               "Y", "Z")))
# pts = matrix(rep(pts,100000), byrow = FALSE, ncol = 3)
# # With java
# start.time <- Sys.time()
# t = xform_brain(pts,
#             sample = "JRCFIB2018F",
#             reference = "FAFB14",
#             .parallel = FALSE,
#             Verbose = TRUE,
#             OmitFailures = TRUE,
#             progress.rjava=TRUE,
#             method = "rjava")
# end.time <- Sys.time()
# time.taken <- end.time - start.time
# time.taken
# ## Time difference of 1.188741 mins
# # Without
# start.time <- Sys.time()
# t = xform_brain(pts,
#             sample = "JRCFIB2018F",
#             reference = "FAFB14",
#             .parallel = FALSE,
#             Verbose = TRUE,
#             OmitFailures = TRUE,
#             method = "auto")
# end.time <- Sys.time()
# time.taken <- end.time - start.time
# time.taken
# ## Time difference of 1.048975 mins
# library(nat.jrcbrains)
# library(hemibrainr)
# neurons = neuprint_read_neurons(hemibrainr::pn.ids[1:10])
# testx=function(method='rjava'){
#   xform_brain(neurons,
#                           sample = "JRCFIB2018F",
#                           reference = "FAFB14",
#                           .parallel = FALSE,
#                           Verbose = TRUE,
#                           OmitFailures = TRUE,
#                           method = method)
#   TRUE
# }
# bench::mark(java=testx('java'), rjava=testx('rjava'))
# #


# Code to transform flywire neurons under construction by the Cambridge FlyConnectome group
library(hemibrainr)
library(nat.jrcbrains)
library(doMC)
library(rJava)

# Parallelise
numCores <- parallel::detectCores()
numCores <- numCores/2
doMC::registerDoMC(numCores)
message("Using ", numCores, " cores")

# Get flywire neurons from google drive
# fw.neurons = hemibrainr::flywire_neurons()
load(file = paste0("/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/flywire/flywire_neurons.rda"))

# Strip meshes for transformed neurons
fw.neurons.skels = hemibrainr:::strip_meshes(fw.neurons)
good = sapply(fw.neurons.skels,nat::is.neuron)
fw.neurons.skels = fw.neurons.skels[good]

# Mirror neurons
fw.neurons.fafb.m = mirror_fafb(fw.neurons.skels, .parallel = FALSE, OmitFailures = TRUE)
hemibrainr:::googledrive_upload_neuronlistfh(fw.neurons.fafb.m,
                                             hemibrainr_team_drive(),
                                             file_name = "flywire_neurons_FlyWire_mirrored.rds",
                                             folder = "flywire_neurons",
                                             subfolder = "FlyWire")
1

# Make the transforms that rely on nat.jrcbrains first
fw.neurons.JRCFIB2018F = hemibrainr:::xform_brain_parallel(fw.neurons.skels, reference = "JRCFIB2018F", sample = "FAFB14", .parallel = FALSE, verbose = TRUE, OmitFailures = TRUE, progress.rjava=TRUE, numCores = numCores)
save(fw.neurons.JRCFIB2018F, file = paste0("/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/flywire/flywire_neurons_JRCFIB2018F.rda"))
hemibrainr:::googledrive_upload_neuronlistfh(fw.neurons.JRCFIB2018F,
                                             hemibrainr_team_drive(),
                                             file_name = "flywire_neurons_JRCFIB2018F.rds",
                                             folder = "flywire_neurons",
                                             subfolder = "JRCFIB2018F")

# Transform into hemibrain space, mirrored
fw.neurons.JRCFIB2018F.m = hemibrainr:::xform_brain_parallel(fw.neurons.fafb.m, reference = "JRCFIB2018F", sample = "FAFB14", .parallel = FALSE, verbose = TRUE, OmitFailures = TRUE, progress.rjava=TRUE, numCores = numCores)
save(fw.neurons.JRCFIB2018F.m, file = paste0("/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/flywire/flywire_neurons_JRCFIB2018F_mirrored.rda"))
hemibrainr:::googledrive_upload_neuronlistfh(fw.neurons.JRCFIB2018F.m,
                                             hemibrainr_team_drive(),
                                             file_name = "flywire_neurons_JRCFIB2018F_mirrored.rds",
                                             folder = "flywire_neurons",
                                             subfolder = "JRCFIB2018F")

# Transform into JFRC2 space
fw.neurons.JFRC2 = xform_brain(fw.neurons.skels, reference = "JFRC2", sample = "FAFB14", .parallel = FALSE, verbose = TRUE, OmitFailures = TRUE)
hemibrainr:::googledrive_upload_neuronlistfh(fw.neurons.JFRC2,
                                hemibrainr_team_drive(),
                                file_name = "flywire_neurons_JFRC2.rds",
                                folder = "flywire_neurons",
                                subfolder = "JFRC2")
fw.neurons.JFRC2.m = nat.templatebrains::mirror_brain(fw.neurons.JFRC2, brain = JFRC2, .parallel = FALSE, OmitFailures = TRUE)
hemibrainr:::googledrive_upload_neuronlistfh(fw.neurons.JFRC2.m,
                                hemibrainr_team_drive(),
                                file_name = "flywire_neurons_JFRC2_mirrored.rds",
                                folder = "flywire_neurons",
                                subfolder = "JFRC2")

# Transform into JRC2018F
fw.neurons.JRC2018F = hemibrainr:::xform_brain_parallel(fw.neurons.skels, reference = "JRC2018F", sample = "FAFB14", .parallel = FALSE, verbose = TRUE, OmitFailures = TRUE, progress.rjava=TRUE, numCores = numCores)
hemibrainr:::googledrive_upload_neuronlistfh(fw.neurons.JRC2018F,
                                             hemibrainr_team_drive(),
                                             file_name = "flywire_neurons_JRC2018F.rds",
                                             folder = "flywire_neurons",
                                             subfolder = "JRC2018F")
# fw.neurons.JRC2018F.m = mirror_brain(fw.neurons.JRC2018F, brain = JRC2018F, .parallel = FALSE, OmitFailures = TRUE)
fw.neurons.JRC2018F.m = xform_brain(fw.neurons.JFRC2.m, reference = "JRC2018F", sample = "JFRC2", .parallel = FALSE, verbose = TRUE, OmitFailures = TRUE, progress.rjava=TRUE)
hemibrainr:::googledrive_upload_neuronlistfh(fw.neurons.JRC2018F.m,
                                             hemibrainr_team_drive(),
                                             file_name = "flywire_neurons_JRC2018F_mirrored.rds",
                                             folder = "flywire_neurons",
                                             subfolder = "JRC2018F")

# Transform into FCWB space
fw.neurons.FCWB = xform_brain(fw.neurons.skels, reference = "FCWB", sample = "FAFB14", .parallel = FALSE, verbose = TRUE, OmitFailures = TRUE)
hemibrainr:::googledrive_upload_neuronlistfh(fw.neurons.FCWB,
                                             hemibrainr_team_drive(),
                                             file_name = "flywire_neurons_FCWB.rds",
                                             folder = "flywire_neurons",
                                             subfolder = "FCWB")
fw.neurons.FCWB.m =  nat.templatebrains::mirror_brain(fw.neurons.FCWB, brain = FCWB, .parallel = FALSE, OmitFailures = TRUE)
hemibrainr:::googledrive_upload_neuronlistfh(fw.neurons.FCWB.m,
                                             hemibrainr_team_drive(),
                                             file_name = "flywire_neurons_FCWB_mirrored.rds",
                                             folder = "flywire_neurons",
                                             subfolder = "FCWB")
