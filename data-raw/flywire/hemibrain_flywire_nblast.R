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
nblastfolder = paste0("/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/NBLAST/")
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
numCores <- numCores/4
doMC::registerDoMC(numCores)
message("Using ", numCores, " cores")

# Get flywire neurons from Google drive
# fw.neurons.JRCFIB2018F = flywire_neurons(brain= "JRCFIB2018F")
# fw.neurons.JRCFIB2018F.m = flywire_neurons(brain= "JRCFIB2018F", mirror = TRUE)
load("/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/flywire/flywire_neurons_JRCFIB2018F.rda")
load("/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/flywire/flywire_neurons_JRCFIB2018F_mirrored.rda")

# Prune
fw.neurons.JRCFIB2018F=nat::nlapply(fw.neurons.JRCFIB2018F, nat::prune_twigs, twig_length=.5, .parallel = TRUE, OmitFailures = TRUE)
fw.neurons.JRCFIB2018F.m=nat::nlapply(fw.neurons.JRCFIB2018F.m, nat::prune_twigs, twig_length=.5, .parallel = TRUE, OmitFailures = TRUE)

# Create vector clouds
fw.neurons.JRCFIB2018F.dps = nat::dotprops(fw.neurons.JRCFIB2018F, .parallel = TRUE, OmitFailures = TRUE)
fw.neurons.JRCFIB2018F.m.dps = nat::dotprops(fw.neurons.JRCFIB2018F.m, .parallel = TRUE, OmitFailures = TRUE)

# Cut to fit hemibrain bounding box

# Load split neurons
polypre = TRUE
split = "synapses"
mode = "centrifugal"
identifier = paste(ifelse(polypre,"polypre","pre"),mode,split,sep="_")
message("Splitting with paramaters: ", identifier)
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
hemibrain.flywire.mean.1 = nat.nblast::nblast(query = fw.neurons.JRCFIB2018F.dps,
                                                 target = all.neurons.flow.microns.dps,
                                                 .parallel=TRUE,
                                                 normalised = TRUE)
hemibrain.flywire.mean.2 = nat.nblast::nblast(query = all.neurons.flow.microns.dps,
                                                 target = fw.neurons.JRCFIB2018F.dps,
                                                 .parallel=TRUE,
                                                 normalised = TRUE)
hemibrain.flywire.mean.native = (hemibrain.flywire.mean.1+t(hemibrain.flywire.mean.2))/2

### NBLAST mirrored
hemibrain.flywire.mean.3 = nat.nblast::nblast(query = fw.neurons.JRCFIB2018F.m.dps,
                                              target = all.neurons.flow.microns.dps,
                                              .parallel=TRUE,
                                              normalised = TRUE)
hemibrain.flywire.mean.4 = nat.nblast::nblast(query = all.neurons.flow.microns.dps,
                                              target = fw.neurons.JRCFIB2018F.m.dps,
                                              .parallel=TRUE,
                                              normalised = TRUE)
hemibrain.flywire.mean.m = (hemibrain.flywire.mean.3+t(hemibrain.flywire.mean.4))/2

## Average native and mirrored
hemibrain.flywire.mean = plyr::rbind.fill.matrix(t(hemibrain.flywire.mean.native), t(hemibrain.flywire.mean.m))
hemibrain.flywire.mean[is.na(hemibrain.flywire.mean)] = 0
rownames(hemibrain.flywire.mean) = c(colnames(hemibrain.flywire.mean.native), colnames(hemibrain.flywire.mean.m))
hemibrain.flywire.mean = hemibrainr:::collapse_matrix_by_names(hemibrain.flywire.mean, FUN = max)

# Save NBLAST on hemibrain googledrive
hemibrainr:::googledrive_upload_nblast(hemibrain.flywire.mean)
save(hemibrain.flywire.mean, file = "/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/NBLAST/hemibrain.flywire.mean.rda")

# Just FlyWire left-right
flywire.mirror.mean.1 = nat.nblast::nblast(query = fw.neurons.JRCFIB2018F.dps,
                                              target = fw.neurons.JRCFIB2018F.m.dps,
                                              .parallel=TRUE,
                                              normalised = TRUE)
flywire.mirror.mean.2 = nat.nblast::nblast(query = fw.neurons.JRCFIB2018F.m.dps,
                                              target = fw.neurons.JRCFIB2018F.dps,
                                              .parallel=TRUE,
                                              normalised = TRUE)
flywire.mirror.mean = plyr::rbind.fill.matrix(flywire.mirror.mean.1+t(flywire.mirror.mean.2))/2
flywire.mirror.mean[is.na(flywire.mirror.mean)] = 0
rownames(flywire.mirror.mean) = colnames(flywire.mirror.mean.2)
save(flywire.mirror.mean, file = "/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/NBLAST/flywire.mirror.mean.rda")
hemibrainr:::googledrive_upload_nblast(flywire.mirror.mean)

# Just FlyWire left-right
flywire.mean = nat.nblast::nblast_allbyall(x = fw.neurons.JRCFIB2018F.dps,
                                                              .progress='text',
                                                              .parallel=TRUE,
                                                              normalisation='mean')
save(flywire.mean, file = "/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/NBLAST/flywire.mean.rda")
hemibrainr:::googledrive_upload_nblast(flywire.mean)
