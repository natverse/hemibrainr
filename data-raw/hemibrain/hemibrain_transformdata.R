###### Script to split all neurons in the hemibrain into putative axon / dendrite #######
######## Originally run from JData5, the Jefferis lab drive on Hal at the MRC LMB #######
# make sure java is in your path
# java --version
# tell R about your current version java
# R CMD javareconf
# setenv PATH ${PATH}:/public/flybrain/hdf5/lib
# setenv LD_LIBRARY_PATH /public/flybrain/hdf5/lib
# source /public/flybrain/flybrain.csh

# Parameters
polypre = TRUE
split = "synapses"
mode = "centrifugal"
identifier = paste(ifelse(polypre,"polypre","pre"),mode,split,sep="_")
message("Splitting with paramaters: ", identifier)

# Set things ups
message("Start")
library(hemibrainr)
library(nat.jrcbrains)
library(doMC)
message("packages loaded")

# Parallelise
numCores <- parallel::detectCores()
doMC::registerDoMC(numCores/2)
message("Using ", numCores/2, " cores")

# Load pre-processed neurons
## "/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/hemibrain_all_neurons_flow_polypre_centrifugal_synapses.rda"
load(paste0("/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/hemibrain_all_neurons_flow_",identifier,".rda"))
message("Checked ", length(all.neurons.flow), " neurons")
all.neurons.flow.microns = hemibrainr::scale_neurons(all.neurons.flow)
good = sapply(all.neurons.flow.microns,nat::is.neuron)
all.neurons.flow.microns = all.neurons.flow.microns[good]

### Checked neurons ###

# Move neurons into different spaces
dir.create("/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/FAFB14/data/",recursive = TRUE)
dir.create("/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/JRC2018F/data/",recursive = TRUE)
dir.create("/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/JFRC2/data/",recursive = TRUE)
dir.create("/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/FCWB/data/",recursive = TRUE)
dir.create("/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/JRCFIB2018F/data/",recursive = TRUE)
dir.create("/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/JRCFIB2018Fraw/data/",recursive = TRUE)
dir.create("/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/FlyWire/data/",recursive = TRUE)

# FAFB14
all.neurons.fafb = hemibrainr:::java_xform_brain(all.neurons.flow.microns, reference = "FAFB14", sample = "JRCFIB2018F", .parallel = FALSE, verbose = TRUE, OmitFailures = TRUE, method = "rjava", progress.rjava=TRUE)
if(length(all.neurons.fafb)<length(all.neurons.flow.microns)){
  all.neurons.fafb = xform_brain(all.neurons.flow.microns, reference = "FAFB14", sample = "JRCFIB2018F", .parallel = FALSE, verbose = TRUE, OmitFailures = TRUE, method = "rjava", progress.rjava=TRUE)
}
all.neurons.fafb = as.neuronlistfh(all.neurons.fafb, dbdir= "/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/FAFB14/data", WriteObjects = "yes")
write.neuronlistfh(all.neurons.fafb, file="/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/FAFB14/hemibrain_all_neurons_flow_FAFB14.rds", overwrite=TRUE)

# Mirror in FAFB14
all.neurons.fafb.m = mirror_fafb(all.neurons.fafb, .parallel = FALSE, OmitFailures = TRUE)
all.neurons.fafb.m = as.neuronlistfh(all.neurons.fafb.m, dbdir= "/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/FAFB14/data", WriteObjects = "yes")
write.neuronlistfh(all.neurons.fafb.m, file="/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/FAFB14/hemibrain_all_neurons_flow_FAFB14_mirrored.rds", overwrite=TRUE)

# JRC2018F
all.neurons.jrc2018f = hemibrainr:::java_xform_brain(all.neurons.flow.microns, reference = "JRC2018F", sample = "JRCFIB2018F", .parallel = FALSE, verbose = TRUE, OmitFailures = TRUE, progress.rjava=TRUE)
all.neurons.jrc2018f = as.neuronlistfh(all.neurons.jrc2018f, dbdir= "/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/JRC2018F/data", WriteObjects = "yes")
write.neuronlistfh(all.neurons.jrc2018f, file="/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/JRC2018F/hemibrain_all_neurons_flow_JRC2018F.rds", overwrite=TRUE)

# JFRC2
all.neurons.jfrc2 = xform_brain(all.neurons.fafb, reference = "JFRC2", sample = "FAFB14", .parallel = FALSE, verbose = TRUE, OmitFailures = TRUE)
all.neurons.jfrc2 = as.neuronlistfh(all.neurons.jfrc2, dbdir= "/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/JFRC2/data", WriteObjects = "yes")
write.neuronlistfh(all.neurons.jfrc2, file="/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/JFRC2/hemibrain_all_neurons_flow_JFRC2.rds", overwrite=TRUE)

# FCWB
all.neurons.fcwb = xform_brain(all.neurons.fafb, reference = "FCWB", sample = "FAFB14", .parallel = FALSE, verbose = TRUE, OmitFailures = TRUE)
all.neurons.fcwb = as.neuronlistfh(all.neurons.fcwb, dbdir= "/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/FCWB/data", WriteObjects = "yes")
write.neuronlistfh(all.neurons.fcwb, file="/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/FCWB/hemibrain_all_neurons_flow_FCWB.rds", overwrite=TRUE)

# Mirror
all.neurons.jfrc2.m = mirror_brain(all.neurons.jfrc2, brain = JFRC2, .parallel = FALSE, verbose = TRUE, OmitFailures = TRUE)
all.neurons.jrc2018f.m = xform_brain(all.neurons.jfrc2.m, reference = "JRC2018F", sample = "JFRC2", .parallel = FALSE, verbose = TRUE, OmitFailures = TRUE)
all.neurons.flow.microns.m = xform_brain(all.neurons.jfrc2.m, reference = "JRCFIB2018F", sample = "JFRC2", .parallel = FALSE, verbose = TRUE, OmitFailures = TRUE)
all.neurons.flow.m = scale_neurons(all.neurons.flow, scaling = (1000/8))
all.neurons.fcwb.m = mirror_brain(all.neurons.fcwb, brain = FCWB, .parallel = FALSE, verbose = TRUE, OmitFailures = TRUE)
all.neurons.flow.microns.m.dps = dotprops(all.neurons.flow.microns)

# Hemibrain, microns
all.neurons.flow.microns.m.dps = nat::dotprops(all.neurons.flow.microns.m, Labels = TRUE, resample = 1, OmitFailures = TRUE, .parallel = FALSE)
all.neurons.flow.microns = as.neuronlistfh(all.neurons.flow.microns, dbdir= "/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/JRCFIB2018F/data", WriteObjects = "yes")
write.neuronlistfh(all.neurons.flow.microns, file="/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/JRCFIB2018F/hemibrain_all_neurons_flow_microns.rds", overwrite=TRUE)
all.neurons.flow.microns.m.dps = as.neuronlistfh(all.neurons.flow.microns.m.dps, dbdir= "/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/JRC2018F/data/", WriteObjects = "yes")
write.neuronlistfh(all.neurons.flow.microns.m.dps, file="/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/JRC2018F/hemibrain_all_neurons_flow_dotprops_mirrored.rds", overwrite=TRUE)

# And mirrored versions
all.neurons.jrc2018f.m = as.neuronlistfh(all.neurons.jrc2018f.m, dbdir= "/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/JRC2018F/data", WriteObjects = "yes")
write.neuronlistfh(all.neurons.jrc2018f.m, file="/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/JRC2018F/hemibrain_all_neurons_flow_JRC2018F_mirrored.rds", overwrite=TRUE)
all.neurons.jfrc2.m = as.neuronlistfh(all.neurons.jfrc2.m, dbdir= "/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/JFRC2/data", WriteObjects = "yes")
write.neuronlistfh(all.neurons.jfrc2.m, file="/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/JFRC2/hemibrain_all_neurons_flow_JFRC2_mirrored.rds", overwrite=TRUE)
all.neurons.fcwb.m = as.neuronlistfh(all.neurons.fcwb.m, dbdir= "/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/FCWB/data", WriteObjects = "yes")
write.neuronlistfh(all.neurons.fcwb.m, file="/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/FCWB/hemibrain_all_neurons_flow_FCWB_mirrored.rds", overwrite=TRUE)
all.neurons.flow.m = as.neuronlistfh(all.neurons.flow.m, dbdir= "/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/JRC2018raw/data", WriteObjects = "yes")
write.neuronlistfh(all.neurons.flow.m, file="/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/JRC2018raw/hemibrain_all_neurons_flow_microns.rds", overwrite=TRUE)







