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

# Extract compartments
all.neurons.pnt = primary_neurite_cable(x = all.neurons.flow, .parallel = TRUE, OmitFailures = TRUE)
all.neurons.pd = primary_dendrite_cable(x = all.neurons.flow, .parallel = TRUE, OmitFailures = TRUE)
all.neurons.axon = axonic_cable(x = all.neurons.flow, .parallel = TRUE, OmitFailures = TRUE)
all.neurons.dendrite = dendritic_cable(x = all.neurons.flow, .parallel = TRUE, OmitFailures = TRUE)
all.neurons.arbour = arbour_cable(x = all.neurons.flow, .parallel = TRUE, OmitFailures = TRUE)

# Make dps objects
## hemibrain.twigs5.dps=dotprops(all.neurons.flow/1e3, k=5, resample=1, .parallel=T, OmitFailures = T)
hemibrain.dps.pnt=dotprops(all.neurons.pnt/1e3, k=5, resample=1, .parallel=T, OmitFailures = T)
hemibrain.dps.pd=dotprops(all.neurons.pd/1e3, k=5, resample=1, .parallel=T, OmitFailures = T)
hemibrain.dps.axon=dotprops(all.neurons.axon/1e3, k=5, resample=1, .parallel=T, OmitFailures = T)
hemibrain.dps.dendrite=dotprops(all.neurons.dendrite/1e3, k=5, resample=1, .parallel=T, OmitFailures = T)
hemibrain.dps.arbour=dotprops(all.neurons.arbour/1e3, k=5, resample=1, .parallel=T, OmitFailures = T)

# Forward FIB->FIB nblasts (normalised)
# hemibrain.twigs5.aba.mean=nblast_allbyall(hemibrain.twigs5.dps,
#                                     .progress='text',
#                                     .parallel=TRUE,
#                                     normalisation='mean')
# hemibrainr:::save_compressed_nblast_mat(hemibrain.twigs5.aba.mean,
#                                         file = nblastfolder)

### NBLAST pnt
hemibrain.pnt.aba.mean=nblast_allbyall(hemibrain.dps.pnt,
                                    .progress='text',
                                    .parallel=TRUE,
                                    normalisation='mean')
hemibrainr:::save_compressed_nblast_mat(hemibrain.pnt.aba.mean,
                                        file = nblastfolder)

### NBLAST pd
hemibrain.pd.aba.mean=nblast_allbyall(hemibrain.dps.pd,
                                        .progress='text',
                                        .parallel=TRUE,
                                        normalisation='mean')
hemibrainr:::save_compressed_nblast_mat(hemibrain.pd.aba.mean,
                                        file = nblastfolder)

### NBLAST axon
hemibrain.axon.aba.mean=nblast_allbyall(hemibrain.dps.axon,
                                        .progress='text',
                                        .parallel=TRUE,
                                        normalisation='mean')
hemibrainr:::save_compressed_nblast_mat(hemibrain.axon.aba.mean,
                                        file = nblastfolder)

### NBLAST dendrite
hemibrain.dendrite.aba.mean=nblast_allbyall(hemibrain.dps.dendrite,
                                        .progress='text',
                                        .parallel=TRUE,
                                        normalisation='mean')
hemibrainr:::save_compressed_nblast_mat(hemibrain.dendrite.aba.mean,
                                        file = nblastfolder)

### NBLAST arbour
hemibrain.arbour.aba.mean=nblast_allbyall(hemibrain.dps.arbour,
                                        .progress='text',
                                        .parallel=TRUE,
                                        normalisation='mean')
hemibrainr:::save_compressed_nblast_mat(hemibrain.arbour.aba.mean,
                                        file = nblastfolder)


