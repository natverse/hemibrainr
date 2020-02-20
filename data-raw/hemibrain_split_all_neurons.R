###### Script to split all neurons in the hemibrain into putative axon / dendrite #######
######## Originally run from JData5, the Jefferis lab drive on Max at the MRC LMB #######

# Set things ups
message("Start")
library(hemibrainr)
library(parallel)
library(doMC)
message("packages loaded")

# Parallelise
numCores <- detectCores()
doMC::registerDoMC(numCores/2)
message("Using ", numCores/2, " cores")

# Get the neuronlist
message("Searching neuprint ...")
all.neurons.meta = neuprint_search("Traced",field="status")
all.neurons.meta = subset(all.neurons.meta, statusLabel!="Leaves")
all.bodyids = all.neurons.meta$bodyid
message("Neuron bodyids: ", length(all.bodyids))

# Read all neurons
all.neurons = neuprint_read_neurons(all.bodyids,.parallel = TRUE)
save(all.neurons, file = "/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/hemibrain_all_neurons.rda")
message("Neurons read: ", length(all.neurons))

# Load meshes
hemibrain.rois = hemibrain_roi_meshes()
message("meshes loaded: ", length(hemibrain.rois))

# Preprocess neurons
all.neurons.checked = hemibrain_skeleton_check(all.neurons, meshes = hemibrain.rois, OmitFailures = TRUE, .parallel = TRUE)
message("Checked ", length(all.neurons.checked), " neurons")

# Save
save(all.neurons.checked, file = "/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/hemibrain_all_neurons_flow_checked_3.rda")

# Split skeletons
all.neurons.flow = flow_centrality(all.neurons.checked, mode = "centrifugal", polypre = FALSE, split = "distance", .parallel = TRUE, OmitFailures = TRUE)
message("Flow centrality calculated for ", length(all.neurons.flow), " neurons")

# Save
save(all.neurons.flow, file = "/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/hemibrain_all_neurons_flow_3.rda")

# Save key points
pns.points = hemibrain_splitpoints(all.neurons.flow)
write.csv(pns.points, file = "/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/hemibrain_all_neurons_splitpoints_3.csv")
message("Split points calculated for neurons: ", nrow(pns.points))

# Calculate and save interesting metrics on this split
mets = hemibrain_compartment_metrics(all.neurons.flow, OmitFailures = TRUE, .parallel = TRUE)
write.csv(mets, file = "/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/hemibrain_all_neurons_metrics_3.csv")
message("Metrics calculated for neurons: ", nrow(mets))


