###### Script to split NBALST all neurons and their compartments  #######
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
doMC::registerDoMC(numCores/2)
message("Using ", numCores/2, " cores")

# Load split neurons
load(paste0("/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/hemibrain_all_neurons_flow_",identifier,".rda"))
message("Loaded: ", length(all.neurons.flow), " neurons")

# Scale
all.neurons.flow.microns.microns = hemibrainr:::scale_neurons.neuronlist(all.neurons.flow, .parallel = TRUE)
message("Neurons scaled!")

# Extract compartments
all.neurons.pnt = primary_neurite_cable(x = all.neurons.flow.microns, .parallel = TRUE, OmitFailures = TRUE)
all.neurons.pd = primary_dendrite_cable(x = all.neurons.flow.microns, .parallel = TRUE, OmitFailures = TRUE)
all.neurons.axon = axonic_cable(x = all.neurons.flow.microns, .parallel = TRUE, OmitFailures = TRUE)
all.neurons.dendrite = dendritic_cable(x = all.neurons.flow.microns, .parallel = TRUE, OmitFailures = TRUE)
all.neurons.arbour = arbour_cable(x = all.neurons.flow.microns, .parallel = TRUE, OmitFailures = TRUE)
all.neurons.spine = nat::nlapply(X = all.neurons.flow.microns, FUN = nat::simplify_neuron, n = 1, invert = FALSE, .parallel = TRUE, OmitFailures = TRUE)

# Make dps objects
## hemibrain.twigs5.dps=dotprops(all.neurons.flow.microns/1e3, k=5, resample=1, .parallel=T, OmitFailures = T)
hemibrain.dps.pnt=dotprops(all.neurons.pnt, k=5, resample=1, .parallel=T, OmitFailures = T)
hemibrain.dps.pd=dotprops(all.neurons.pd, k=5, resample=1, .parallel=T, OmitFailures = T)
hemibrain.dps.axon=dotprops(all.neurons.axon, k=5, resample=1, .parallel=T, OmitFailures = T)
hemibrain.dps.dendrite=dotprops(all.neurons.dendrite, k=5, resample=1, .parallel=T, OmitFailures = T)
hemibrain.dps.arbour=dotprops(all.neurons.arbour, k=5, resample=1, .parallel=T, OmitFailures = T)
hemibrain.dps.spine=dotprops(all.neurons.spine, k=5, resample=1, .parallel=T, OmitFailures = T)

# Forward FIB->FIB nblasts (normalised)
# hemibrain.twigs5.aba.mean=nat.nblast::nblast_allbyall(hemibrain.twigs5.dps,
#                                     .progress='text',
#                                     .parallel=TRUE,
#                                     normalisation='mean')
# hemibrainr:::save_compressed_nblast_mat(hemibrain.twigs5.aba.mean,
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
