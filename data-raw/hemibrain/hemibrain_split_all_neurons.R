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

# Get functions
source("/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/R/startup/functions.R")
message("functions loaded")

# Parallelise
numCores <- parallel::detectCores()
doMC::registerDoMC(numCores/2)
message("Using ", numCores/2, " cores")

# Load pre-processed neurons
load("/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/hemibrain_all_neurons_flow_checked.rda")
message("Checked ", length(all.neurons.checked), " neurons")

# Split skeletons
all.neurons.flow = flow_centrality(all.neurons.checked, mode = mode, polypre = polypre, split = split, .parallel = TRUE, OmitFailures = TRUE)
missed = setdiff(names(all.neurons.checked),names(all.neurons.flow))
if(length(missed)){
  missed.flow = hemibrainr::add_Label(all.neurons.checked[as.character(missed)], Label = 3, internal.assignments = TRUE)
  all.neurons.flow = all.neurons.flow[setdiff(names(all.neurons.flow),names(missed.flow))]
  all.neurons.flow = c(all.neurons.flow,missed.flow)
}
message("Flow centrality calculated for ", length(all.neurons.flow), " neurons")
sensories = unique(c(rn.ids,orn.ids,hrn.ids))
rn.flow = hemibrainr::add_Label(all.neurons.flow[intersect(as.character(sensories),names(all.neurons.flow))], Label = 2, internal.assignments = TRUE)
all.neurons.flow[names(rn.flow)] = rn.flow
good = sapply(all.neurons.flow,nat::is.neuron)
all.neurons.flow = all.neurons.flow[good]
rm("all.neurons.checked")

# Make sure connectors also have labels
all.neurons.flow = nat::nlapply(all.neurons.flow, hemibrainr:::carryover_labels, .parallel = TRUE)

# Save
save(all.neurons.flow, file = paste0("/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/hemibrain_all_neurons_flow_",identifier,".rda"))

# Save as neuronlistfh
savedirfh = paste0("/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/hemibrain_all_neurons_flow_",identifier,"/")
unlink(paste0(savedirfh,"data/"))
dir.create = dir.create(paste0(savedirfh,"data/"),recursive = TRUE)
all.neurons.flow.fh = as.neuronlistfh(all.neurons.flow, dbdir= "/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/JRCFIB2018Fraw/data/", WriteObjects = "yes")
write.neuronlistfh(all.neurons.flow.fh, file="/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/JRCFIB2018Fraw/hemibrain_all_neurons_flow_JRCFIB2018Fraw.rds", overwrite=TRUE)
hemibrainr:::remove_unused_filehash("/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/JRCFIB2018Fraw/")

# Save key points
all.neurons.points = hemibrain_splitpoints(all.neurons.flow)
write.csv(all.neurons.points, file = paste0("/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/hemibrain_all_neurons_splitpoints_",identifier,".csv"), row.names=FALSE)
message("Split points: ", nrow(all.neurons.points))

# Save compressed version
all.neurons.points.comp = subset(all.neurons.points, !grepl("null",all.neurons.points$point))
write.csv(all.neurons.points.comp, file = paste0("/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/hemibrain_all_neurons_splitpoints_",identifier,"_compressed.csv"), row.names=FALSE)
message("Compressed split points: ", nrow(all.neurons.points.comp))

# Get Synapses
all.neurons.flow.syns = hemibrain_extract_synapses(all.neurons.flow, .parallel = TRUE, OmitFailures = TRUE)
write.csv(all.neurons.flow.syns, file = paste0("/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/hemibrain_all_neurons_synapses_",identifier,".csv"), row.names=FALSE)
message("Saved synapses: ", nrow(all.neurons.flow.syns))

# Get connections
all.neurons.flow.conns = hemibrain_extract_connections(all.neurons.flow, .parallel = TRUE, OmitFailures = TRUE)
write.csv(all.neurons.flow.conns, file = paste0("/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/hemibrain_all_neurons_connections_",identifier,".csv"), row.names=FALSE)
message("Saved connections: ", nrow(all.neurons.flow.conns))

# Get edgelist
all.neurons.flow.elist = hemibrain_extract_compartment_edgelist(all.neurons.flow, .parallel = TRUE, OmitFailures = TRUE)
write.csv(all.neurons.flow.elist, file = paste0("/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/hemibrain_all_neurons_edgelist_",identifier,".csv"), row.names=FALSE)
message("Saved connections: ", nrow(all.neurons.flow.elist))

# Calculate and save interesting metrics on this split
all.neurons.flow.microns = all.neurons.flow*(8/1000) # convert to microns
mets = hemibrain_compartment_metrics(all.neurons.flow.microns, OmitFailures = TRUE, .parallel = TRUE, delta = 5, resample = NULL)
mets2 = lapply(1:ncol(mets), function(x) unlist(mets[,x]))
mets2 = as.data.frame(t(do.call(rbind, mets2)))
dimnames(mets2) = dimnames(mets)
all.neurons.flow[,] = mets2[names(all.neurons.flow),]
all.neurons.flow.microns[,] = mets2[names(all.neurons.flow),]
write.csv(mets2, file = paste0("/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/hemibrain_all_neurons_metrics_",identifier,".csv"), row.names=FALSE)
message("Metrics calculated for neurons: ", nrow(mets))
save(all.neurons.flow, file = paste0("/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/hemibrain_all_neurons_flow_",identifier,".rda"))

# Make a dotprops object
all.neurons.flow.dps = nat::dotprops(all.neurons.flow.microns, Labels = TRUE, resample = 1, OmitFailures = TRUE, .parallel = TRUE)
save(all.neurons.flow.dps, file = paste0("/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/hemibrain_all_neurons_flow_",identifier,"_dps.rda"))
all.neurons.flow.dps.fh = as.neuronlistfh(all.neurons.flow.dps, dbdir= "/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/JRCFIB2018F/data/", WriteObjects = "yes")
write.neuronlistfh(all.neurons.flow.dps.fh, file="/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/JRCFIB2018F/hemibrain_all_neurons_flow_JRCFIB2018F_dotprops.rds", overwrite=TRUE)
hemibrainr:::remove_unused_filehash("/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/JRCFIB2018F/")



