###### Script to split all neurons in the hemibrain into putative axon / dendrite #######
######## Originally run from JData5, the Jefferis lab drive on Max at the MRC LMB #######

###### Script to split all neurons in the hemibrain into putative axon / dendrite #######
######## Originally run from JData5, the Jefferis lab drive on Max at the MRC LMB #######

## use the python here: /lmb/home/abates/miniconda3/bin/python

### blender
# scp Downloads/blender-2.90.0-linux64.zip abates@atg.mrc-lmb.cam.ac.uk:/home/abates/blender-2.90.0-linux64.zip
# ssh abates@atg.mrc-lmb.cam.ac.uk
# scp blender-2.90.0-linux64.zip abates@max:/public/flybrain/blender-2.90.0-linux64.zip
# ssh max
# unzip /public/flybrain/blender-2.90.0-linux64.zip
# Add /public/flybrain/Blender.app to your PATH:
# setenv PATH ${PATH}:/public/flybrain/blender-2.90.0-linux64/
# setenv LD_LIBRARY_PATH /public/flybrain/R/R-3.5.1/lib64/R/library/Rvcg/libs/
# source /public/flybrain/flybrain.csh
# source /public/gcc/gcc7_1_0.csh
### blender

# Find blender!
reticulate::py_run_string("import os")
# Check what's in PATH
reticulate::py_run_string("print(os.environ['PATH'])")
# This is how trimesh tries to find Blender (on Linux):
reticulate::py_run_string("from distutils.spawn import find_executable")
reticulate::py_run_string("blender_executable = find_executable('blender')")
reticulate::py_run_string("print('Blender executable:', blender_executable)")

# Code to skeletonise flywire neurons under construction by the Cambridge FlyConnectome group
library(hemibrainr)
library(fafbseg)
library(nat.jrcbrains)
library(foreach)
library(doParallel)
fafbseg::choose_segmentation("flywire")

# Parallelise
numCores <- detectCores()
numCores <- numCores/2
registerDoParallel(numCores)  # use multicore, set to the number of our cores
doMC::registerDoMC(numCores)
registerDoParallel(numCores)  # use multicore, set to the number of our cores
message("Cores recruited: ", foreach::getDoParWorkers())

# Neuron positions flagged by the community as of interest
master = hemibrainr:::flywire_ids_update(numCores=numCores)
1
ids = unique(master$flywire.id)

# Load neurons we have previously downloaded
## Only update the ones that need updating
load(file = paste0("/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/flywire/flywire_neurons.rda"))
old.ids = names(fw.neurons)
fw.neurons.old = fw.neurons[intersect(ids,old.ids)]
ids = setdiff(ids, old.ids)

# Read neurons
batches = split(ids, round(seq(from = 1, to = numCores, length.out = length(ids))))
foreach.skeletons <- foreach::foreach (batch = 1:numCores) %dopar% {
  fw.ids = batches[[batch]]
  j = tryCatch(fafbseg::skeletor(fw.ids, clean = FALSE, brain = elmr::FAFB14.surf, mesh3d = FALSE, SL = 10), error = function(e) NULL)
}
isnl = sapply(foreach.skeletons, nat::is.neuronlist)
fw.neurons = do.call(c, foreach.skeletons[isnl])
fw.neurons = union(fw.neurons, fw.neurons.old)

# Add metadata
colnames(fw.neurons[,]) = c("flywire.id")
coords = master[match(names(fw.neurons),master$flywire.id),]
fw.neurons[,] = coords
fw.neurons[,"dataset"] = "flywire"

# Save locally just in case
save(fw.neurons, file = paste0("/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/flywire/flywire_neurons.rda"))

# Save flywire skeletons
hemibrainr:::googledrive_upload_neuronlistfh(fw.neurons,
                                team_drive = hemibrainr_team_drive(),
                                file_name = "flywire_neurons_FlyWire.rds",
                                folder = "flywire_neurons",
                                subfolder = "FlyWire",
                                numCores = numCores)
1
