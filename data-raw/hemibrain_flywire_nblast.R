###### Script to split all neurons in the hemibrain into putative axon / dendrite #######
######## Originally run from JData5, the Jefferis lab drive on Max at the MRC LMB #######

## use the python here: /lmb/home/abates/miniconda3/bin/python

### blender
# scp Downloads/Blender.zip abates@atg.mrc-lmb.cam.ac.uk:/home/abates/Blender.zip
# ssh abates@atg.mrc-lmb.cam.ac.uk
# scp Blender.zip abates@max:/public/flybrain/Blender.zip
### blender

# Set things ups
message("Start")
nblastfolder = paste0("/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/NBLAST/",identifier,"/")
dir.create(nblastfolder, recursive = TRUE)
source("/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/R/startup/packages.R")
library(hemibrainr)
library(fafbseg)
library(nat.jrcbrains)
library(foreach)
library(doParallel)
message("packages loaded")

# Parallelise
numCores <- detectCores()
numCores <- numCores/2
registerDoParallel(numCores)  # use multicore, set to the number of our cores
doMC::registerDoMC(numCores)
registerDoParallel(numCores)  # use multicore, set to the number of our cores
message("Cores recruited: ", foreach::getDoParWorkers())

# Get flywire neurons from google drive
fw.neurons = flywire_neurons(brain= "JRCFIB2018F")
fw.neurons.m = flywire_neurons(brain= "JRCFIB2018F", mirror = TRUE)
fw.neurons.dps = dotprops(fw.neurons, .parallel = TRUE, OmitFailures = TRUE)
fw.neurons.m.dps = dotprops(fw.neurons.m, .parallel = TRUE, OmitFailures = TRUE)

# Load split neurons
load(paste0("/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/hemibrain_all_neurons_flow_",identifier,".rda"))
message("Loaded: ", length(all.neurons.flow), " neurons")

# Scale
fib.twigs5=nlapply(all.neurons.flow, prune_twigs, twig_length=5000, .parallel = TRUE, OmitFailures = TRUE)
all.neurons.flow.microns = hemibrainr:::scale_neurons.neuronlist(fib.twigs5, .parallel = TRUE, OmitFailures = TRUE)
all.neurons.flow.microns.dps = dotprops(all.neurons.flow.microns, .parallel = TRUE, OmitFailures = TRUE)
rm("all.neurons.flow")
rm("fib.twigs5")
message("Neurons scaled!")

### NBLAST native
hemibrain.flywire.mean.1 = nat.nblast::nblast(query = fw.neurons.dps,
                                                 target = all.neurons.flow.microns.dps,
                                                 .parallel=TRUE,
                                                 normalised = TRUE)
hemibrain.flywire.mean.2 = nat.nblast::nblast(query = all.neurons.flow.microns.dps,
                                                 target = fw.neurons.dps,
                                                 .parallel=TRUE,
                                                 normalised = TRUE)
hemibrain.flywire.mean.native = (hemibrain.flywire.mean.1+t(hemibrain.flywire.mean.2))/2

### NBLAST mirrored
hemibrain.flywire.mean.3 = nat.nblast::nblast(query = fw.neurons.m.dps,
                                              target = all.neurons.flow.microns.dps,
                                              .parallel=TRUE,
                                              normalised = TRUE)
hemibrain.flywire.mean.4 = nat.nblast::nblast(query = all.neurons.flow.microns.dps,
                                              target = fw.neurons.m.dps,
                                              .parallel=TRUE,
                                              normalised = TRUE)
hemibrain.flywire.mean.m = (hemibrain.flywire.mean.3+t(hemibrain.flywire.mean.4))/2

## Average native and mirrored
hemibrain.flywire.mean = rbind(hemibrain.flywire.mean.native, hemibrain.flywire.mean.m)
hemibrain.flywire.mean = hemibrainr:::collapse_matrix_by_names(hemibrain.flywire.mean)

# Save NBLAST
hemibrainr:::save_compressed_nblast_mat(hemibrain.flywire.mean,
                                        file = nblastfolder)

# Move files to googledrive



