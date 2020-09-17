###### Script to split all neurons in the hemibrain into putative axon / dendrite #######
######## Originally run from JData5, the Jefferis lab drive on Hex at the MRC LMB #######

# source /public/gcc/gcc7_1_0.csh
# In order to then update Rvcg in R if needed

# At the LMB, this code must be run on Hex/Hal not Max

# Code to transform flywire neurons under construction by the Cambridge FlyConnectome group
library(hemibrainr)
library(nat.jrcbrains)
library(foreach)
library(doParallel)

# Parallelise
numCores <- parallel::detectCores()
numCores <- numCores/2
doMC::registerDoMC(numCores)
message("Using ", numCores, " cores")

# Get flywire neurons from google drive
fw.neurons = hemibrainr::flywire_neurons()
load(file = paste0("/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/flywire/flywire_neurons.rda"))

# Strip meshes for transformed neurons
fw.neurons.skels = hemibrainr:::strip_meshes(fw.neurons)

# Make the transforms that rely in nat.jrcbrains first
fw.neurons.JRCFIB2018F = xform_brain(fw.neurons.skels, reference = "JRCFIB2018F", sample = "FAFB14", .parallel = TRUE, verbose = TRUE)
hemibrainr:::googledrive_upload_neuronlistfh(fw.neurons.JRCFIB2018F,
                                             team_drive = "hemibrain",
                                             file_name = "flywire_neurons_JRCFIB2018F.rds",
                                             folder = "flywire_neurons",
                                             subfolder = "JRCFIB2018F")
1
fw.neurons.JRC2018F = xform_brain(fw.neurons.skels, reference = "JRC2018F", sample = "FAFB14", .parallel = TRUE, verbose = TRUE)
hemibrainr:::googledrive_upload_neuronlistfh(fw.neurons.JRC2018F,
                                             team_drive = "hemibrain",
                                             file_name = "flywire_neurons_JRC2018F.rds",
                                             folder = "flywire_neurons",
                                             subfolder = "JRC2018F")
save(fw.neurons.JRCFIB2018F, file = paste0("/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/flywire/flywire_neurons_JRCFIB2018F.rda"))

# Transform into JFRC2 space
fw.neurons.JFRC2 = xform_brain(fw.neurons.skels, reference = "JFRC2", sample = "FAFB14", .parallel = TRUE, verbose = TRUE)
hemibrainr:::googledrive_upload_neuronlistfh(fw.neurons.JFRC2,
                                team_drive = "hemibrain",
                                file_name = "flywire_neurons_JFRC2.rds",
                                folder = "flywire_neurons",
                                subfolder = "JFRC2")
fw.neurons.JFRC2.m = mirror_brain(fw.neurons.JFRC2, brain = JFRC2, .parallel = TRUE)
hemibrainr:::googledrive_upload_neuronlistfh(fw.neurons.JFRC2.m,
                                team_drive = "hemibrain",
                                file_name = "flywire_neurons_JFRC2_mirrored.rds",
                                folder = "flywire_neurons",
                                subfolder = "JFRC2")
fw.neurons.fafb.m = xform_brain(fw.neurons.JFRC2.m, reference = "FAFB14", sample = "JFRC2", .parallel = TRUE, verbose = TRUE)
hemibrainr:::googledrive_upload_neuronlistfh(fw.neurons.fafb.m,
                                team_drive = "hemibrain",
                                file_name = "flywire_neurons_FlyWire_mirrored.rds",
                                folder = "flywire_neurons",
                                subfolder = "FlyWire")

# Transform into FCWB space
fw.neurons.FCWB = xform_brain(fw.neurons.skels, reference = "FCWB", sample = "FAFB14", .parallel = TRUE, verbose = TRUE)
hemibrainr:::googledrive_upload_neuronlistfh(fw.neurons.FCWB,
                                team_drive = "hemibrain",
                                file_name = "flywire_neurons_FCWB.rds",
                                folder = "flywire_neurons",
                                subfolder = "FCWB")
fw.neurons.FCWB.m = xform_brain(fw.neurons.JFRC2.m, reference = "FCWB", sample = "JFRC2", .parallel = TRUE, verbose = TRUE)
hemibrainr:::googledrive_upload_neuronlistfh(fw.neurons.FCWB.m,
                                team_drive = "hemibrain",
                                file_name = "flywire_neurons_FCWB_mirrored.rds",
                                folder = "flywire_neurons",
                                subfolder = "FCWB")

# Transform into hemibrain space
fw.neurons.JRCFIB2018F.m = xform_brain(fw.neurons.JFRC2.m, reference = "JRCFIB2018F", sample = "JFRC2", .parallel = TRUE, verbose = TRUE)
hemibrainr:::googledrive_upload_neuronlistfh(fw.neurons.JRCFIB2018F.m,
                                team_drive = "hemibrain",
                                file_name = "flywire_neurons_JRCFIB2018F_mirrored.rds",
                                folder = "flywire_neurons",
                                subfolder = "JRCFIB2018F")
save(fw.neurons.JRCFIB2018F.m, file = paste0("/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/flywire/flywire_neurons_JRCFIB2018F_mirrored.rda"))

# Transform into JRC2018F space
fw.neurons.JRC2018F.m = xform_brain(fw.neurons.JFRC2.m, reference = "JRC2018F", sample = "JFRC2", .parallel = TRUE, verbose = TRUE)
hemibrainr:::googledrive_upload_neuronlistfh(fw.neurons.JRC2018F.m,
                                             team_drive = "hemibrain",
                                             file_name = "flywire_neurons_JRC2018F_mirrored.rds",
                                             folder = "flywire_neurons",
                                             subfolder = "JRC2018F")
