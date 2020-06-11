############## Script to split NBALST all neurons and their compartments  ###############
######## Originally run from JData5, the Jefferis lab drive on Max at the MRC LMB #######

# Parameters
polypre = TRUE
split = "synapses"
mode = "centrifugal"
identifier = paste(ifelse(polypre,"polypre","pre"),mode,split,sep="_")
message("Working with paramaters: ", identifier)

# Set things ups
message("Start")
nblastfolder = paste0("/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/NBLAST/",identifier,"/")
dir.create(nblastfolder, recursive = TRUE)
source("/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/R/startup/packages.R")
message("packages loaded")

# Get functions
source("/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/R/startup/functions.R")
message("functions loaded")

# Parallelise
numCores <- detectCores()
doMC::registerDoMC(numCores/4)
message("Using ", numCores/4, " cores")

# Load split neurons
load(paste0("/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/hemibrain_all_neurons_flow_",identifier,".rda"))
message("Loaded: ", length(all.neurons.flow), " neurons")

# Scale
fib.twigs5=nlapply(all.neurons.flow, prune_twigs, twig_length=5000, .parallel = TRUE, OmitFailures = TRUE)
all.neurons.flow.microns = hemibrainr:::scale_neurons.neuronlist(fib.twigs5, .parallel = TRUE, OmitFailures = TRUE)
rm("all.neurons.flow")
rm("fib.twigs5")
message("Neurons scaled!")

# Extract compartments
all.neurons.pnt = primary_neurite_cable(x = all.neurons.flow.microns, .parallel = TRUE, OmitFailures = TRUE)
all.neurons.pd = primary_dendrite_cable(x = all.neurons.flow.microns, .parallel = TRUE, OmitFailures = TRUE)
all.neurons.axon = axonic_cable(x = all.neurons.flow.microns, .parallel = TRUE, OmitFailures = TRUE)
all.neurons.dendrite = dendritic_cable(x = all.neurons.flow.microns, .parallel = TRUE, OmitFailures = TRUE)
all.neurons.arbour = arbour_cable(x = all.neurons.flow.microns, .parallel = TRUE, OmitFailures = TRUE)
all.neurons.simp = nat::nlapply(X = all.neurons.flow.microns, FUN = nat::simplify_neuron, n = 1, invert = FALSE, .parallel = TRUE, OmitFailures = TRUE)
all.neurons.tract = tract_cable(x = all.neurons.flow.microns, .parallel = TRUE, OmitFailures = TRUE)
all.neurons.spine = nat::nlapply(X = all.neurons.flow.microns, FUN = nat::spine, n = 1, UseStartPoint = TRUE, .parallel = TRUE, OmitFailures = TRUE)

# Make dps objects
hemibrain.dps=dotprops(all.neurons.flow.microns, k=5, resample=1, .parallel=T, OmitFailures = T)
hemibrain.dps.pnt=dotprops(all.neurons.pnt, k=5, resample=1, .parallel=T, OmitFailures = T)
rm("all.neurons.pnt")
hemibrain.dps.pd=dotprops(all.neurons.pd, k=5, resample=1, .parallel=T, OmitFailures = T)
rm("all.neurons.pd")
hemibrain.dps.axon=dotprops(all.neurons.axon, k=5, resample=1, .parallel=T, OmitFailures = T)
rm("all.neurons.axon")
hemibrain.dps.dendrite=dotprops(all.neurons.dendrite, k=5, resample=1, .parallel=T, OmitFailures = T)
rm("all.neurons.dendrite")
hemibrain.dps.arbour=dotprops(all.neurons.arbour, k=5, resample=1, .parallel=T, OmitFailures = T)
rm("all.neurons.arbour")
hemibrain.dps.spine=dotprops(all.neurons.spine, k=5, resample=1, .parallel=T, OmitFailures = T)
rm("all.neurons.spine")
hemibrain.dps.simp=dotprops(all.neurons.simp, k=5, resample=1, .parallel=T, OmitFailures = T)
rm("all.neurons.simp")
hemibrain.dps.tract=dotprops(all.neurons.tract, k=5, resample=1, .parallel=T, OmitFailures = T)
rm("all.neurons.tract")

# Forward FIB->FIB nblasts (normalised)
# hemibrain.aba.mean=nat.nblast::nblast_allbyall(hemibrain.twigs5.dps,
#                                     .progress='text',
#                                     .parallel=TRUE,
#                                     normalisation='mean')
# hemibrainr:::save_compressed_nblast_mat(hemibrain.aba.mean,
#                                         file = nblastfolder)

### NBLAST pnt
hemibrain.pnt.aba.mean=nat.nblast::nblast_allbyall(hemibrain.dps.pnt,
                                    .progress='text',
                                    .parallel=TRUE,
                                    normalisation='mean')
hemibrainr:::save_compressed_nblast_mat(hemibrain.pnt.aba.mean,
                                        file = nblastfolder)

### NBLAST pd
hemibrain.pd.aba.mean=nat.nblast::nblast_allbyall(hemibrain.dps.pd,
                                        .progress='text',
                                        .parallel=TRUE,
                                        normalisation='mean')
hemibrainr:::save_compressed_nblast_mat(hemibrain.pd.aba.mean,
                                        file = nblastfolder)

### NBLAST axon
hemibrain.axon.aba.mean=nat.nblast::nblast_allbyall(hemibrain.dps.axon,
                                        .progress='text',
                                        .parallel=TRUE,
                                        normalisation='mean')
hemibrainr:::save_compressed_nblast_mat(hemibrain.axon.aba.mean,
                                        file = nblastfolder)

### NBLAST dendrite
hemibrain.dendrite.aba.mean=nat.nblast::nblast_allbyall(hemibrain.dps.dendrite,
                                        .progress='text',
                                        .parallel=TRUE,
                                        normalisation='mean')
hemibrainr:::save_compressed_nblast_mat(hemibrain.dendrite.aba.mean,
                                        file = nblastfolder)

### NBLAST arbour
hemibrain.arbour.aba.mean=nat.nblast::nblast_allbyall(hemibrain.dps.arbour,
                                        .progress='text',
                                        .parallel=TRUE,
                                        normalisation='mean')
hemibrainr:::save_compressed_nblast_mat(hemibrain.arbour.aba.mean,
                                        file = nblastfolder)

### NBLAST spine
hemibrain.spine.aba.mean=nat.nblast::nblast_allbyall(hemibrain.dps.spine,
                                                      .progress='text',
                                                      .parallel=TRUE,
                                                      normalisation='mean')
hemibrainr:::save_compressed_nblast_mat(hemibrain.spine.aba.mean,
                                        file = nblastfolder)

### NBLAST simp
hemibrain.simp.aba.mean=nat.nblast::nblast_allbyall(hemibrain.dps.simp,
                                                     .progress='text',
                                                     .parallel=TRUE,
                                                     normalisation='mean')
hemibrainr:::save_compressed_nblast_mat(hemibrain.simp.aba.mean,
                                        file = nblastfolder)

### NBLAST tract
hemibrain.tract.aba.mean=nat.nblast::nblast_allbyall(hemibrain.dps.tract,
                                                     .progress='text',
                                                     .parallel=TRUE,
                                                     normalisation='mean')
hemibrainr:::save_compressed_nblast_mat(hemibrain.tract.aba.mean,
                                        file = nblastfolder)


####################
# Light level data #
####################

# Light level neurons
lhfolder = paste0("/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/light_level/lhns/")
flycircuitfolder = paste0("/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/light_level/flycircuit/")
most.lhns.hemi.dps = hemibrain_lm_lhns(savedir = lhfolder,
                                       folder="",
                                       local = TRUE,
                                       data = "dps",
                                       brainspace = "JRCFIB2018F")
most.lhins.hemi.dps = hemibrainr::hemibrain_lm_lhns(savedir = lhfolder,
                                       folder="",
                                       local = TRUE,
                                       data = "dps",
                                       cable = "lhins",
                                       brainspace = "JRCFIB2018F")
most.lh.hemi.dps = nat::union(most.lhns.hemi.dps,most.lhins.hemi.dps)
check = unlist(sapply(most.lh.hemi.dps, function(e) nrow(e$points)>2))
most.lh.hemi.dps = most.lh.hemi.dps[check]
fcs.hemi.dps = flycircuit_neurons(savedir = flycircuitfolder,
                                  folder="",
                                  local = TRUE,
                                  cable = "all",
                                  data = "dps",
                                  brainspace = "JRCFIB2018F")

### NBLAST with light level
hemibrain.lhns.mean.1 = nat.nblast::nblast(query = most.lh.hemi.dps,
                                         target = hemibrain.dps,
                                         .parallel=TRUE,
                                         normalised = TRUE)
hemibrain.lhns.mean.2 = nat.nblast::nblast(query = hemibrain.dps,
                                         target = most.lh.hemi.dps,
                                         .parallel=TRUE,
                                         normalised = TRUE)
hemibrain.lhns.mean = (hemibrain.lhns.mean.1+t(hemibrain.lhns.mean.2))/2
hemibrainr:::save_compressed_nblast_mat(hemibrain.lhns.mean,
                                        file = nblastfolder)

### NBLAST with FlyCircuit
hemibrain.flycircuit.mean.1 = nat.nblast::nblast(query = fcs.hemi.dps,
                                         target = hemibrain.dps,
                                         .parallel=TRUE,
                                         normalised = TRUE)
hemibrain.flycircuit.mean.2 = nat.nblast::nblast(query = hemibrain.dps,
                                               target = fcs.hemi.dps,
                                               .parallel=TRUE,
                                               normalised = TRUE)
hemibrain.flycircuit.mean = (hemibrain.flycircuit.mean.1+t(hemibrain.flycircuit.mean.2))/2
hemibrainr:::save_compressed_nblast_mat(hemibrain.flycircuit.mean,
                                        file = nblastfolder)



