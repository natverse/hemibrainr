###### Script to split all neurons in the hemibrain into putative axon / dendrite #######
######## Originally run from JData5, the Jefferis lab drive on Max at the MRC LMB #######

# Parameters
polypre = TRUE
split = "distance"
mode = "centrifugal"
identifier = paste(ifelse(polypre,"polypre","pre"),mode,split,sep="_")
message("Splitting with paramaters: ", identifier)

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
#message("Searching neuprint ...")
#all.neurons.meta = neuprint_search("Traced",field="status")
#all.neurons.meta = subset(all.neurons.meta, statusLabel!="Leaves")
#all.bodyids = all.neurons.meta$bodyid
#message("Neuron bodyids: ", length(all.bodyids))

# Read all neurons
#all.neurons = neuprint_read_neurons(all.bodyids,.parallel = TRUE)
#save(all.neurons, file = "/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/hemibrain_all_neurons.rda")
#load("/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/hemibrain_all_neurons.rda")
#message("Neurons read: ", length(all.neurons))

# Load meshes
#hemibrain.rois = hemibrain_roi_meshes()
#message("meshes loaded: ", length(hemibrain.rois))

# Preprocess neurons
#all.neurons.checked = hemibrain_skeleton_check(all.neurons, meshes = hemibrain.rois, OmitFailures = TRUE, .parallel = TRUE)
#message("Checked ", length(all.neurons.checked), " neurons")

# Save
#save(all.neurons.checked, file = "/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/hemibrain_all_neurons_flow_checked_3.rda")
load("/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/hemibrain_all_neurons_flow_checked.rda")
message("Checked ", length(all.neurons.checked), " neurons")

# Split skeletons
all.neurons.flow = flow_centrality(all.neurons.checked, mode = mode, polypre = polypre, split = split, .parallel = TRUE, OmitFailures = TRUE)
message("Flow centrality calculated for ", length(all.neurons.flow), " neurons")

# Save
save(all.neurons.flow, file = paste0("/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/hemibrain_all_neurons_flow_",identifier,".rda"))

# Save key points
all.neurons.points = hemibrain_splitpoints(all.neurons.flow)
write.csv(all.neurons.points, file = paste0("/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/hemibrain_all_neurons_splitpoints_",identifier,".csv"))
message("Split points calculated for neurons: ", nrow(all.neurons.points))

# Get Synapses
all.neurons.flow.syns = hemibrain_extract_synapses(all.neurons.flow)
write.csv(all.neurons.flow.syns, file = paste0("/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/hemibrain_all_neurons_synapses_",identifier,".csv"))
message("Saved synapses: ", nrow(all.neurons.flow.syns))

# Get connections
all.neurons.flow.conns = hemibrain_extract_connections(all.neurons.flow)
write.csv(all.neurons.flow.conns, file = paste0("/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/hemibrain_all_neurons_connections_",identifier,".csv"))
message("Saved connections: ", nrow(all.neurons.flow.conns))

# Calculate and save interesting metrics on this split
all.neurons.flow = all.neurons.flow*(8/1000) # convert to microns
mets = hemibrain_compartment_metrics(all.neurons.flow, OmitFailures = TRUE, .parallel = TRUE, delta = 5, resample = NULL)
write.csv(mets, file = paste0("/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/hemibrain_all_neurons_metrics_",identifier,".csv"))
message("Metrics calculated for neurons: ", nrow(mets))

