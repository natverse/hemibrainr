###### Script to split all neurons in the hemibrain into putative axon / dendrite #######
######## Originally run from JData5, the Jefferis lab drive on Max at the MRC LMB #######

# setenv LD_LIBRARY_PATH /public/flybrain/R/R-3.5.1/lib64/R/library/rJava/libs/
# source /public/flybrain/flybrain.csh
# source /public/gcc/gcc7_1_0.csh

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
                kc.ids,
                apl.ids,
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
                lhn.ids,
                ton.ids,
                lc.ids
                )
all.bodyids = unique(all.bodyids)
write.csv(all.bodyids, file = "/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/hemibrain_all_neuron_bodyids.csv", row.names=FALSE)
message("Neuron bodyids: ", length(all.bodyids))

# Read all neurons
all.neurons = neuprint_read_neurons(all.bodyids,.parallel = TRUE, all_segments = TRUE, OmitFailures = TRUE)
missing = setdiff(all.bodyids, names(all.neurons))
if(length(missing)){
  missing.neurons = neuprint_read_neurons(missing,.parallel = TRUE, all_segments = TRUE, heal = 10000, OmitFailures = TRUE)
  all.neurons = c(all.neurons, missing.neurons)
}
missing = setdiff(all.bodyids, names(all.neurons))
if(length(missing)){
  missing.neurons = neuprint_read_neurons(missing,.parallel = TRUE, all_segments = TRUE, heal = FALSE, OmitFailures = TRUE)
  all.neurons = c(all.neurons, missing.neurons)
}
good = sapply(all.neurons, function(x) nrow(x$d)>2)
all.neurons = all.neurons[good]
all.neurons = add_field_seq(all.neurons,all.neurons[,"bodyid"],field="bodyid")
save(all.neurons, file = "/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/hemibrain_all_neurons.rda")
message("Neurons read: ", length(all.neurons))

# Load meshes
hemibrain.rois = hemibrain_roi_meshes()
message("meshes loaded: ", length(hemibrain.rois))

# Preprocess neurons
all.neurons.checked = hemibrain_skeleton_check(all.neurons, meshes = hemibrain.rois, OmitFailures = TRUE, .parallel = TRUE, googlesheet = TRUE)
1
message("Checked ", length(all.neurons.checked), " neurons")

# Get the right meta data
meta = hemibrain_get_meta(x = names(all.neurons.checked))
meta = meta[,]
all.neurons.checked[,] = meta[names(all.neurons.checked),]
all.neurons.checked = hemibrainr:::hemibrain_neuron_class(all.neurons.checked)

# Save
save(all.neurons.checked, file = "/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/hemibrain_all_neurons_flow_checked.rda")
message("Checked ", length(all.neurons.checked), " neurons")

