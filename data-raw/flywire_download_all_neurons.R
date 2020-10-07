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
### blender

# Find blender!
reticulate::py_run_string("import os")
# Check what's in PATH
reticulate::py_run_string("print(os.environ['PATH'])")
# This is how trimesh tries to find Blender (on Linux):
reticulate::py_run_string("from distutils.spawn import find_executable")
reticulate::py_run_string("blender_executable = find_executable('blender')")
reticulate::py_run_string("print('Blender executable:', blender_executable)")
# Update skeletor
reticulate::py_run_string("pip3 install git+git://github.com/schlegelp/skeletor@master")

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
positions.gi = hemibrainr:::gsheet_manipulation(FUN = googlesheets4::read_sheet,
                                        ss = "1rzG1MuZYacM-vbW7100aK8HeA-BY6dWAVXQ7TB6E2cQ",
                                        sheet = "flywire",
                                        guess_max = 1,
                                        return = TRUE)
1
positions.gi$flywire.xyz = apply(positions.gi[,c("fw.x","fw.y",'fw.z')],1,paste,sep=",",collapse=",")

## Read matching google sheet
tabs = googlesheets4::sheet_names("1spGSuhUX6Hhn-8HH0U_ArIWUuPpMBFNjIjeSSh_MFVY")
tabs = tabs[grepl("ItoLee",tabs)]
gs = data.frame()
for(tab in tabs){
  gs.t = hemibrainr:::gsheet_manipulation(FUN = googlesheets4::read_sheet,
                                          ss = "1spGSuhUX6Hhn-8HH0U_ArIWUuPpMBFNjIjeSSh_MFVY",
                                          sheet = tab,
                                          guess_max = 3000,
                                          return = TRUE)
  gs = plyr::rbind.fill(gs,gs.t)
}
gs = gs[!is.na(gs$fw.x),]
gs$flywire.xyz = apply(gs[,c("fw.x","fw.y",'fw.z')],1,paste,sep=",",collapse=",")

# Get xyz positions for FAFBv14 neurons.
fg = hemibrainr:::gsheet_manipulation(FUN = googlesheets4::read_sheet,
                         ss = "1OSlDtnR3B1LiB5cwI5x5Ql6LkZd8JOS5bBr-HTi0pOw",
                         sheet = "FAFB",
                         guess_max = 3000,
                         return = TRUE)
positions.gs = fg[!is.na(fg$flywire.xyz),]
positions.gs = sapply(positions.gs$flywire.xyz,strsplit,",")
positions.gs = positions.gs[sapply(positions.gs,function(p) length(p)==3)]
positions.gs = do.call(rbind, positions.gs)
positions.gs = apply(positions.gs,2,as.numeric)

# Get positions for flywire neurons=
positions = gs[,c("fw.x","fw.y",'fw.z')]
positions = apply(positions,2,as.numeric)
positions = rbind(positions, positions.gs)
positions = rbind(positions, positions.gi)
positions = as.data.frame(positions)
positions$flywire.xyz = apply(positions[,c("fw.x","fw.y",'fw.z')],1,paste,sep=",",collapse=",")
positions = positions[!duplicated(positions$flywire.xyz),]

# Batch positions for grabbing IDs from flywire
batches = split(1:nrow(positions), round(seq(from = 1, to = numCores, length.out = nrow(positions))))

# Find flywire IDs for each specified position
foreach.ids <- foreach::foreach (batch = 1:numCores) %dopar% {
  pos = positions[batches[[batch]],]
  colnames(pos) = c("x","y","z","flywire.xyz")
  j <- tryCatch({i = fafbseg::flywire_xyz2id(pos[,1:3], rawcoords = TRUE)
  names(i) = pos$flywire.xyz
  if(sum(i==0)>0){
    i[i==0] =  fafbseg::flywire_xyz2id(pos[i==0,], rawcoords = FALSE)
  }
  i}, error = function(e) NULL)
}
j[]
# Load neurons we have previously downloaded
## Only update the ones that need updating
save(fw.neurons, file = paste0("/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/flywire/flywire_neurons.rda"))
old.ids = names(fw.neurons)
fw.neurons.old = fw.neurons[intersect(ids,old.ids)]
ids = setdiff(ids, old.ids)

# Read neurons
batches = split(ids, round(seq(from = 1, to = numCores, length.out = length(ids))))
foreach.skeletons <- foreach::foreach (batch = 1:numCores) %dopar% {
  fw.ids = batches[[batch]]
  j = tryCatch(fafbseg::skeletor(fw.ids, clean = FALSE, brain = elmr::FAFB14.surf, mesh3d = FALSE), error = function(e) NULL)
}
isnl = sapply(foreach.skeletons, nat::is.neuronlist)
fw.neurons = do.call(c, foreach.skeletons[isnl])
fw.neurons = union(fw.neurons, fw.neurons.old)

# Add metadata
colnames(fw.neurons[,]) = c("flywire.id")
coords = master[match(names(fw.neurons),master$flywire.id),]
flywire.xyz = apply(coords[,c("fw.x","fw.y",'fw.z')],1,paste,sep=",",collapse=",")
fw.neurons[,"flywire.xyz"] = flywire.xyz
fw.neurons[,"FAFB.xyz"] = coords[,"FAFB.xyz"]
fw.neurons[,"side"] = coords[,"side"]
fw.neurons[,"ItoLee_Hemilineage"] = coords[,"ItoLee_Hemilineage"]
fw.neurons[,"Hartenstein_Hemilineage"] = coords[,"Hartenstein_Hemilineage"]
fw.neurons[,"skid"] = coords[,"skid"]
fw.neurons[,"hemibrain_match"] = coords[,"hemibrain_match"]
fw.neurons[,"dataset"] = "flywire"

# Save locally just in case
save(fw.neurons, file = paste0("/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/flywire/flywire_neurons.rda"))

# Save flywire skeletons
hemibrainr:::googledrive_upload_neuronlistfh(fw.neurons,
                                team_drive = "hemibrain",
                                file_name = "flywire_neurons_FlyWire.rds",
                                folder = "flywire_neurons",
                                subfolder = "FlyWire")
1
