###### Script to split all neurons in the hemibrain into putative axon / dendrite #######
######## Originally run from JData5, the Jefferis lab drive on Max at the MRC LMB #######

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

# Get the neuronlist
message("Searching neuprint ...")
all.neurons.meta = neuprint_search("Traced",field="status")
all.neurons.meta = subset(all.neurons.meta, statusLabel!="Leaves")
all.bodyids = all.neurons.meta$bodyid
all.bodyids = c(all.bodyids,
                rn.ids,
                orn.ids,
                hrn.ids,
                pn.ids,
                upn.ids,
                mpn.ids,
                vppn.ids,
                dan.ids,
                mbon.ids,
                alln.ids,
                ton.ids
                )
all.bodyids = unique(all.bodyids)
write.csv(all.bodyids, file = "/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/hemibrain_all_neuron_bodyids.csv", row.names=FALSE)
message("Neuron bodyids: ", length(all.bodyids))

# Read all neurons
all.neurons = neuprint_read_neurons(all.bodyids,.parallel = TRUE, all_segments = TRUE)
save(all.neurons, file = "/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/hemibrain_all_neurons.rda")
message("Neurons read: ", length(all.neurons))

# Load meshes
hemibrain.rois = hemibrain_roi_meshes()
message("meshes loaded: ", length(hemibrain.rois))

# Preprocess neurons
all.neurons.checked = hemibrain_skeleton_check(all.neurons, meshes = hemibrain.rois, OmitFailures = TRUE, .parallel = TRUE)
message("Checked ", length(all.neurons.checked), " neurons")

# Save
save(all.neurons.checked, file = "/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/hemibrain_all_neurons_flow_checked.rda")
message("Checked ", length(all.neurons.checked), " neurons")
