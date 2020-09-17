###### Script to split all neurons in the hemibrain into putative axon / dendrite #######
######## Originally run from JData5, the Jefferis lab drive on Max at the MRC LMB #######

## use the python here: /lmb/home/abates/miniconda3/bin/python

### blender
# scp Downloads/Blender.zip abates@atg.mrc-lmb.cam.ac.uk:/home/abates/Blender.zip
# ssh abates@atg.mrc-lmb.cam.ac.uk
# scp Blender.zip abates@max:/public/flybrain/Blender.zip
# unzip /public/flybrain/Blender.zip
# setenv PATH ${PATH}:/public/flybrain/blender-2.90.0-linux64/blender
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
numCores <- parallel::detectCores()
numCores <- numCores/2
doMC::registerDoMC(numCores)
message("Using ", numCores, " cores")

# Get flywire neurons from google drive
fw.neurons = flywire_neurons(brain= "JRCFIB2018F")
fw.neurons.m = flywire_neurons(brain= "JRCFIB2018F", mirror = TRUE)

# Prune
fw.neurons=nat::nlapply(fw.neurons, nat::prune_twigs, twig_length=.5, .parallel = TRUE, OmitFailures = TRUE)
fw.neurons.m=nat::nlapply(fw.neurons.m, nat::prune_twigs, twig_length=.5, .parallel = TRUE, OmitFailures = TRUE)

# Create vector clouds
fw.neurons.dps = dotprops(fw.neurons, .parallel = TRUE, OmitFailures = TRUE)
fw.neurons.m.dps = dotprops(fw.neurons.m, .parallel = TRUE, OmitFailures = TRUE)

# Load split neurons
load(paste0("/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/hemibrain_all_neurons_flow_",identifier,".rda"))
message("Loaded: ", length(all.neurons.flow), " neurons")

# Scale
fib.twigs5=nlapply(all.neurons.flow, prune_twigs, twig_length=500, .parallel = TRUE, OmitFailures = TRUE)
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
hemibrain.flywire.mean = hemibrainr:::collapse_matrix_by_names(hemibrain.flywire.mean, FUN = max)

# Save NBLAST on hemibrain googledrive
hemibrainr:::googledrive_upload_nblast(hemibrain.flywire.mean)

# Just Flywire left-right
hemibrain.flywire.mirror.mean.1 = nat.nblast::nblast(query = fw.neurons.dps,
                                              target = fw.neurons.m.dps,
                                              .parallel=TRUE,
                                              normalised = TRUE)
hemibrain.flywire.mirror.mean.2 = nat.nblast::nblast(query = fw.neurons.m.dps,
                                              target = fw.neurons.dps,
                                              .parallel=TRUE,
                                              normalised = TRUE)
hemibrain.flywire.mirror.mean = (hemibrain.flywire.mirror.mean.1+t(hemibrain.flywire.mirror.mean.2))/2
hemibrainr:::googledrive_upload_nblast(hemibrain.flywire.mirror.mean)


# Just Flywire left-right
hemibrain.flywire.mean = nat.nblast::nblast_allbyall(query = fw.neurons.dps,
                                                              .progress='text',
                                                              .parallel=TRUE,
                                                              normalisation='mean')
hemibrainr:::googledrive_upload_nblast(hemibrain.flywire.mean)
