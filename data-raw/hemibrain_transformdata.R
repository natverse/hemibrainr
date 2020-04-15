###### Script to split all neurons in the hemibrain into putative axon / dendrite #######
######## Originally run from JData5, the Jefferis lab drive on Max at the MRC LMB #######

# Parameters
# polypre = TRUE
# split = "synapses"
# mode = "centrifugal"
# identifier = paste(ifelse(polypre,"polypre","pre"),mode,split,sep="_")
# message("Splitting with paramaters: ", identifier)

# Set things ups
message("Start")
source("/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/R/startup/packages.R")
message("packages loaded")

# Get functions
source("/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/R/startup/functions.R")
message("functions loaded")

# Parallelise
numCores <- detectCores()
doMC::registerDoMC(numCores/2)
message("Using ", numCores/2, " cores")

# Load pre-processed neurons
load("/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/hemibrain_all_neurons_flow_checked.rda")
message("Checked ", length(all.neurons.checked), " neurons")

### Checked neurons ###

# Move neurons into different spaces
dir.create("FAFB/data/",recursive = TRUE)
dir.create("JRC2018F/data/",recursive = TRUE)
dir.create("JFRC2/data/",recursive = TRUE)
dir.create("FCWB/data/",recursive = TRUE)
all.neurons.fafb = xform_brain(all.neurons.checked, reference = "FAFB14", sample = "JRCFIB2018Fraw", .parallel = TRUE)
all.neurons.jrc2018f = xform_brain(all.neurons.checked, reference = "JRC2018F", sample = "JRCFIB2018Fraw", .parallel = TRUE)
all.neurons.jfrc2 = xform_brain(all.neurons.checked, reference = "JFRC2", sample = "JRCFIB2018Fraw", .parallel = TRUE)
all.neurons.fcwb = xform_brain(all.neurons.checked, reference = "FCWB", sample = "JRCFIB2018Fraw", .parallel = TRUE)

# Save as neuronlistfh objects
all.neurons.flow.fafb = as.neuronlistfh(all.neurons.flow.fafb, dbdir= "/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/FAFB/data", WriteObjects = "yes")
write.neuronlistfh(all.neurons.flow.fafb, file="/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/FAFB/hemibrain_all_neurons_checked_FAFB.rds", overwrite=TRUE)
all.neurons.jrc2018f = as.neuronlistfh(all.neurons.jrc2018f, dbdir= "/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/JRC2018F/data", WriteObjects = "yes")
write.neuronlistfh(all.neurons.jrc2018f, file="/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/JRC2018F/hemibrain_all_neurons_checked_JRC2018F.rds", overwrite=TRUE)
all.neurons.jfrc2 = as.neuronlistfh(all.neurons.jfrc2, dbdir= "/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/JFRC2/data", WriteObjects = "yes")
write.neuronlistfh(all.neurons.jfrc2, file="/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/JFRC2/hemibrain_all_neurons_checked_JFRC2.rds", overwrite=TRUE)
all.neurons.fcwb = as.neuronlistfh(all.neurons.fcwb, dbdir= "/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/FCWB/data", WriteObjects = "yes")
write.neuronlistfh(all.neurons.fcwb, file="/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/FCWB/hemibrain_all_neurons_checked_FCWB.rds", overwrite=TRUE)

