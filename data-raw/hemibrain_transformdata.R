###### Script to split all neurons in the hemibrain into putative axon / dendrite #######
######## Originally run from JData5, the Jefferis lab drive on Max at the MRC LMB #######

# Parameters
polypre = TRUE
split = "synapses"
mode = "centrifugal"
identifier = paste(ifelse(polypre,"polypre","pre"),mode,split,sep="_")
message("Splitting with paramaters: ", identifier)

# Set things ups
message("Start")
library(nat.jrcbrains)
library(catmaid)
library(hemibrainr)
library(fafbseg)
message("packages loaded")

# Get functions
source("/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/R/startup/functions.R")
message("functions loaded")

# Parallelise
numCores <- detectCores()
doMC::registerDoMC(numCores/2)
message("Using ", numCores/2, " cores")

# Load pre-processed neurons
load(paste0("/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/hemibrain_all_neurons_flow_",identifier,".rda"))
message("Checked ", length(all.neurons.flow), " neurons")

### Checked neurons ###

# Move neurons into different spaces
dir.create("/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/FAFB/data/",recursive = TRUE)
dir.create("/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/JRC2018F/data/",recursive = TRUE)
dir.create("/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/JFRC2/data/",recursive = TRUE)
dir.create("/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/FCWB/data/",recursive = TRUE)
dir.create("/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/JRCFIB2018F/data/",recursive = TRUE)
dir.create("/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/flywire/data/",recursive = TRUE)
all.neurons.flow.microns = scale_neurons(all.neurons.flow)
all.neurons.fafb = xform_brain(all.neurons.flow.microns, reference = "FAFB14", sample = "JRCFIB2018F", .parallel = TRUE, verbose = TRUE, OmitFailures = TRUE)
all.neurons.jrc2018f = xform_brain(all.neurons.flow.microns, reference = "JRC2018F", sample = "JRCFIB2018F", .parallel = FALSE, verbose = TRUE, OmitFailures = TRUE)
all.neurons.jfrc2 = xform_brain(all.neurons.fafb, reference = "JFRC2", sample = "FAFB14", .parallel = FALSE, verbose = TRUE, OmitFailures = TRUE)
all.neurons.fcwb = xform_brain(all.neurons.fafb, reference = "FCWB", sample = "FAFB14", .parallel = FALSE, verbose = TRUE, OmitFailures = TRUE)
all.neurons.flywire = xform_brain(all.neurons.fafb, reference = "FlyWire", sample = "FAFB14", .parallel = TRUE, verbose = TRUE, OmitFailures = TRUE)

# Mirror
all.neurons.fcwb.m = mirror_brain(all.neurons.fcwb, brain = FCWB)
all.neurons.jfrc2.m = mirror_brain(all.neurons.jfrc2, brain = JFRC2)
all.neurons.jrc2018f.m = xform_brain(all.neurons.jfrc2.m, reference = "JRC2018F", sample = "JFRC2", .parallel = FALSE, verbose = TRUE)
all.neurons.fafb.m = xform_brain(all.neurons.jfrc2.m, reference = "FAFB14", sample = "JFRC2", .parallel = FALSE, verbose = TRUE)
all.neurons.flow.microns.m = xform_brain(all.neurons.jfrc2.m, reference = "JRCFIB2018F", sample = "JFRC2", .parallel = FALSE, verbose = TRUE)
all.neurons.flow.m = scale_neurons(all.neurons.flow, scaling = (1000/8))
all.neurons.flywire = xform_brain(all.neurons.fafb.m, reference = "FlyWire", sample = "FAFB14", .parallel = TRUE, verbose = TRUE)

# Save as neuronlistfh objects
all.neurons.fafb = as.neuronlistfh(all.neurons.fafb, dbdir= "/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/FAFB/data", WriteObjects = "yes")
write.neuronlistfh(all.neurons.fafb, file="/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/FAFB/hemibrain_all_neurons_checked_FAFB.rds", overwrite=TRUE)
all.neurons.jrc2018f = as.neuronlistfh(all.neurons.jrc2018f, dbdir= "/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/JRC2018F/data", WriteObjects = "yes")
write.neuronlistfh(all.neurons.jrc2018f, file="/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/JRC2018F/hemibrain_all_neurons_checked_JRC2018F.rds", overwrite=TRUE)
all.neurons.jfrc2 = as.neuronlistfh(all.neurons.jfrc2, dbdir= "/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/JFRC2/data", WriteObjects = "yes")
write.neuronlistfh(all.neurons.jfrc2, file="/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/JFRC2/hemibrain_all_neurons_checked_JFRC2.rds", overwrite=TRUE)
all.neurons.fcwb = as.neuronlistfh(all.neurons.fcwb, dbdir= "/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/FCWB/data", WriteObjects = "yes")
write.neuronlistfh(all.neurons.fcwb, file="/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/FCWB/hemibrain_all_neurons_checked_FCWB.rds", overwrite=TRUE)
all.neurons.flow.microns = as.neuronlistfh(all.neurons.flow.microns, dbdir= "/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/JRCFIB2018F/data", WriteObjects = "yes")
write.neuronlistfh(all.neurons.flow.microns, file="/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/JRCFIB2018F/hemibrain_all_neurons_checked_microns.rds", overwrite=TRUE)
all.neurons.fafb = as.neuronlistfh(all.neurons.flywire, dbdir= "/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/flywire/data", WriteObjects = "yes")
write.neuronlistfh(all.neurons.flywire, file="/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/flywire/hemibrain_all_neurons_checked_FAFB.rds", overwrite=TRUE)

# And mirrored versions
all.neurons.fafb.m = as.neuronlistfh(all.neurons.fafb.m, dbdir= "/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/FAFB/data", WriteObjects = "yes")
write.neuronlistfh(all.neurons.fafb.m, file="/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/FAFB/hemibrain_all_neurons_checked_FAFB_mirrored.rds", overwrite=TRUE)
all.neurons.jrc2018f.m = as.neuronlistfh(all.neurons.jrc2018f.m, dbdir= "/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/JRC2018F/data", WriteObjects = "yes")
write.neuronlistfh(all.neurons.jrc2018f.m, file="/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/JRC2018F/hemibrain_all_neurons_checked_JRC2018F_mirrored.rds", overwrite=TRUE)
all.neurons.jfrc2.m = as.neuronlistfh(all.neurons.jfrc2.m, dbdir= "/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/JFRC2/data", WriteObjects = "yes")
write.neuronlistfh(all.neurons.jfrc2.m, file="/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/JFRC2/hemibrain_all_neurons_checked_JFRC2_mirrored.rds", overwrite=TRUE)
all.neurons.fcwb.m = as.neuronlistfh(all.neurons.fcwb.m, dbdir= "/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/FCWB/data", WriteObjects = "yes")
write.neuronlistfh(all.neurons.fcwb.m, file="/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/FCWB/hemibrain_all_neurons_checked_FCWB_mirrored.rds", overwrite=TRUE)
all.neurons.flow.microns.m = as.neuronlistfh(all.neurons.flow.microns.m, dbdir= "/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/JRCFIB2018F/data", WriteObjects = "yes")
write.neuronlistfh(all.neurons.flow.microns.m, file="/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/JRCFIB2018F/hemibrain_all_neurons_checked_microns_mirrored.rds", overwrite=TRUE)
all.neurons.flow.m = as.neuronlistfh(all.neurons.flow.m, dbdir= "/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/data", WriteObjects = "yes")
write.neuronlistfh(all.neurons.flow.m, file="/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/hemibrain_all_neurons_checked_microns.rds", overwrite=TRUE)
