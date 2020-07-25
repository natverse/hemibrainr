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
source("/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/R/startup/packages.R")
message("packages loaded")
library()

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
dir.create("FAFB/data/",recursive = TRUE)
dir.create("JRC2018F/data/",recursive = TRUE)
dir.create("JFRC2/data/",recursive = TRUE)
dir.create("FCWB/data/",recursive = TRUE)
dir.create("mirror/data/",recursive = TRUE)
all.neurons.fafb = xform_brain(all.neurons.flow, reference = "FAFB14", sample = "JRCFIB2018Fraw", .parallel = TRUE)
all.neurons.jrc2018f = xform_brain(all.neurons.flow, reference = "JRC2018F", sample = "JRCFIB2018Fraw", .parallel = TRUE)
all.neurons.jfrc2 = xform_brain(all.neurons.flow, reference = "JFRC2", sample = "JRCFIB2018Fraw", .parallel = TRUE)
all.neurons.fcwb = xform_brain(all.neurons.flow, reference = "FCWB", sample = "JRCFIB2018Fraw", .parallel = TRUE)

# Mirror
all.neurons.jrc2018f.m = mirror_brain(all.neurons.jrc2018f, brain = JRC2018F)
all.neurons.fcwb.m = mirror_brain(all.neurons.jrc2018f, brain = FCWB)
all.neurons.jfrc2.m = mirror_brain(all.neurons.jrc2018f, brain = JFRC2)
all.neurons.fafb.m = xform_brain(all.neurons.flow, reference = "FAFB14", sample = "JRC2018F", .parallel = TRUE)
all.neurons.flow.m = xform_brain(all.neurons.jrc2018f.m, reference = "JRCFIB2018Fraw", sample = "JRC2018F", .parallel = TRUE)

# Save as neuronlistfh objects
all.neurons.flow.fafb = as.neuronlistfh(all.neurons.flow.fafb, dbdir= "/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/FAFB/data", WriteObjects = "yes")
write.neuronlistfh(all.neurons.flow.fafb, file="/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/FAFB/hemibrain_all_neurons_checked_FAFB.rds", overwrite=TRUE)
all.neurons.jrc2018f = as.neuronlistfh(all.neurons.jrc2018f, dbdir= "/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/JRC2018F/data", WriteObjects = "yes")
write.neuronlistfh(all.neurons.jrc2018f, file="/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/JRC2018F/hemibrain_all_neurons_checked_JRC2018F.rds", overwrite=TRUE)
all.neurons.jfrc2 = as.neuronlistfh(all.neurons.jfrc2, dbdir= "/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/JFRC2/data", WriteObjects = "yes")
write.neuronlistfh(all.neurons.jfrc2, file="/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/JFRC2/hemibrain_all_neurons_checked_JFRC2.rds", overwrite=TRUE)
all.neurons.fcwb = as.neuronlistfh(all.neurons.fcwb, dbdir= "/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/FCWB/data", WriteObjects = "yes")
write.neuronlistfh(all.neurons.fcwb, file="/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/FCWB/hemibrain_all_neurons_checked_FCWB.rds", overwrite=TRUE)

# And mirrored versions
all.neurons.flow.fafb.m = as.neuronlistfh(all.neurons.flow.fafb.m, dbdir= "/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/FAFB/data", WriteObjects = "yes")
write.neuronlistfh(all.neurons.flow.fafb.m, file="/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/FAFB/hemibrain_all_neurons_checked_FAFB_mirrored.rds", overwrite=TRUE)
all.neurons.jrc2018f.m = as.neuronlistfh(all.neurons.jrc2018f.m, dbdir= "/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/JRC2018F/data", WriteObjects = "yes")
write.neuronlistfh(all.neurons.jrc2018f.m, file="/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/JRC2018F/hemibrain_all_neurons_checked_JRC2018F_mirrored.rds", overwrite=TRUE)
all.neurons.jfrc2.m = as.neuronlistfh(all.neurons.jfrc2.m, dbdir= "/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/JFRC2/data", WriteObjects = "yes")
write.neuronlistfh(all.neurons.jfrc2.m, file="/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/JFRC2/hemibrain_all_neurons_checked_JFRC2_mirrored.rds", overwrite=TRUE)
all.neurons.fcwb.m = as.neuronlistfh(all.neurons.fcwb.m, dbdir= "/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/FCWB/data", WriteObjects = "yes")
write.neuronlistfh(all.neurons.fcwb.m, file="/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/FCWB/hemibrain_all_neurons_checked_FCWB_mirrored.rds", overwrite=TRUE)
all.neurons.flow.m = as.neuronlistfh(all.neurons.flow.m, dbdir= "/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/mirror/data", WriteObjects = "yes")
write.neuronlistfh(all.neurons.flow.m, file="/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/mirror/hemibrain_all_neurons_checked_FCWB_mirrored.rds", overwrite=TRUE)
