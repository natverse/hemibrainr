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
# setenv PATH ${PATH}:/public/flybrain/blender-2.90.0-linux64/blender
### blender

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

## Read google sheet
tabs = googlesheets4::sheet_names("1spGSuhUX6Hhn-8HH0U_ArIWUuPpMBFNjIjeSSh_MFVY")
1
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

# Get xyz positions for FAFBv14 neurons.
# positions = gs$FAFB.xyz
# positions = positions[!is.na(positions)]
# positions = sapply(positions,strsplit,",")
# positions = positions[sapply(positions,function(p) length(p)==3)]
# positions = do.call(rbind, positions)
# positions = apply(positions,2,as.numeric)
# positions = xform_brain(positions, reference = "FlyWire", sample = "FAFB14", .parallel = TRUE, verbose = TRUE)

# Get positions for flywire neurons
positions = gs[,c("fw.x","fw.y",'fw.z')]
positions = positions[!is.na(gs$fw.x),]
positions = apply(positions,2,as.numeric)

# Batch positions for grabbing IDs from flywire
batches = split(1:nrow(positions), round(seq(from = 1, to = numCores, length.out = nrow(positions))))

# Find flywire IDs for each specified position
foreach.ids <- foreach::foreach (batch = 1:numCores) %dopar% {
  pos = positions[batches[[batch]],]
  i = fafbseg::flywire_xyz2id(pos, rawcoords = TRUE)
  if(sum(i==0)>0){
    i[i==0] =  fafbseg::flywire_xyz2id(pos[i==0,], rawcoords = FALSE)
  }
  i[!is.na(i)]
}

# Batch IDs for grabbing meshes from flywire
ids = unlist(foreach.ids)
batches = split(ids, round(seq(from = 1, to = numCores, length.out = length(ids))))

# Read neurons
foreach.skeletons <- foreach::foreach (batch = 1:numCores) %dopar% {
  fw.ids = batches[[batch]]
  j = tryCatch(fafbseg::skeletor(fw.ids, clean = FALSE), error = function(e){
    message("Batch failed:", batch)
    NULL
  })
  j[!is.na(j)]
}
fw.neurons = do.call(c, foreach.skeletons)

# Save flywire skeletons
googledrive_upload_neuronlistfh(fw.neurons,
                                team_drive = "hemibrain",
                                file_name = "flywire_neurons_FlyWire.rds",
                                folder = "flywire_neurons",
                                subfolder = "FlyWire")
1
