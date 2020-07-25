###### Script to split all neurons in the hemibrain into putative axon / dendrite #######
######## Originally run from JData5, the Jefferis lab drive on Max at the MRC LMB #######

### rsync oauth
# gargle-oauth abates$ rsync -avzh /Users/abates/.R/gargle/gargle-oauth/ max:/lmb/home/abates/.R/gargle/gargle-oauth
### rsync oauth

# Parameters
polypre = TRUE
split = "synapses"
mode = "centrifugal"
identifier = paste(ifelse(polypre,"polypre","pre"),mode,split,sep="_")
message("Saving splits with paramaters: ", identifier)

# Set things ups
message("Start")
source("/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/R/startup/packages.R")
message("packages loaded")
library(googledrive)

# Google drive location:
hemibrain = team_drive_get("hemibrain")
1
drive_hemibrain = drive_find(type = "folder", team_drive = hemibrain)
hemibrain_neurons= subset(drive_hemibrain,name=="hemibrain_neurons")
hemibrain_nblast= subset(drive_hemibrain,name=="hemibrain_nblast")
flycircuit= subset(drive_hemibrain,name=="flycircuit")
nblast.folder = subset(drive_hemibrain,name==paste0("nblast_",identifier))

# Save NBLAST files
nblast.files  = list.files(paste0("/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/NBLAST/",identifier),pattern = "rda", full.name = TRUE)
for(nblast.file in nblast.files){
  drive_put(media = nblast.file,
               path = nblast.folder,
               overwrite = TRUE,
               verbose = TRUE)
}

# Save all neuron bodyids
drive_put(media = "/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/hemibrain_all_neuron_bodyids.csv",
             path = hemibrain_neurons,
             overwrite = TRUE,
             verbose = TRUE)

# Save raw neurons
drive_put(media = "/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/hemibrain_all_neurons.rda",
             path = hemibrain_neurons,
             overwrite = TRUE,
             verbose = TRUE)

# Save pre-processed neurons
drive_put(media = "/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/hemibrain_all_neurons_flow_checked.rda",
             path = hemibrain_neurons,
             overwrite = TRUE,
             verbose = TRUE)

# Make subfolder, if it does not already exist:
if(!identifier %in% drive_hemibrain$name){
  drive_mkdir(name = identifier,
              path = hemibrain_neurons,
              overwrite = FALSE)
}
folder = subset(drive_hemibrain,name==identifier)
identifier.folder = drive_ls(path = folder[1,], type = "folder", team_drive = hemibrain)

# Split files to save
split.neurons = paste0("/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/hemibrain_all_neurons_flow_",identifier,".rda")
split.neuron.fh = paste0("/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/hemibrain_all_neurons_flow_",identifier,"/")
split.points = paste0("/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/hemibrain_all_neurons_splitpoints_",identifier,".csv")
split.synapses = paste0("/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/hemibrain_all_neurons_synapses_",identifier,".csv")
split.connections = paste0("/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/hemibrain_all_neurons_connections_",identifier,".csv")
split.elist = paste0("/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/hemibrain_all_neurons_edgelist_",identifier,".csv")
split.metrics = paste0("/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/hemibrain_all_neurons_metrics_",identifier,".csv")
split.neuron.fh.neuronlist = list.files(split.neuron.fh,pattern = "rds", full.name = TRUE)
split.neuron.fh.data = list.files(split.neuron.fh,pattern = "data", full.name = TRUE)
split.neuron.fh.data.files = list.files(split.neuron.fh.data, full.names = TRUE)

# Save split data
drive_put(media = split.neurons,
             path = folder,
             overwrite = TRUE,
             verbose = TRUE)
drive_put(media = split.neuron.fh.neuronlist,
             path = folder,
             overwrite = TRUE,
             verbose = TRUE)
drive_put(media = split.points,
             path = folder,
             overwrite = TRUE,
             verbose = TRUE)
drive_put(media = split.synapses,
             path = folder,
             overwrite = TRUE,
             verbose = TRUE)
drive_put(media = split.connections,
             path = folder,
             overwrite = TRUE,
             verbose = TRUE)
drive_put(media = split.elist,
             path = folder,
             overwrite = TRUE,
             verbose = TRUE)
drive_put(media = split.metrics,
             path = folder,
             overwrite = TRUE,
             verbose = TRUE)

# Save neuronlistfh data
identifier.folder.data = subset(identifier.folder, name == "data")
if(nrow(identifier.folder.data)){
  drive_rm(identifier.folder.data, team_drive = hemibrain)
}
drive_mkdir(name = "data",
            path = folder,
            overwrite = FALSE)
identifier.folder = drive_ls(path = folder, type = "folder", team_drive = hemibrain)
identifier.folder.data = subset(identifier.folder, name == "data")
error.files = c()
for(split.neuron.fh.data.file in split.neuron.fh.data.files){
  e = tryCatch(drive_upload(media = split.neuron.fh.data.file,
                            path = identifier.folder.data[1,],
                            overwrite = FALSE,
                            verbose = TRUE),
               error = function(e) split.neuron.fh.data.file)
  if(is.null(e)){
    error.files = c(error.files,split.neuron.fh.data.file)
  }
}
message("File not downloaded: ", length(error.files))

# Re-save missing IDs
identifier.folder.data = subset(identifier.folder, name == "data")
identifier.folder.data.files = drive_ls(identifier.folder.data, team_drive = hemibrain)
missing.files = split.neuron.fh.data.files[!gsub(".*/","",split.neuron.fh.data.files) %in% identifier.folder.data.files$name]
for(m in missing.files){
  e = tryCatch(drive_upload(media = m,
                            path = identifier.folder.data[1,],
                            overwrite = FALSE,
                            verbose = TRUE),
               error = function(e) split.neuron.fh.data.file)
  if(is.null(e)){
    error.files = c(error.files,split.neuron.fh.data.file)
  }
}


## Save transformed neurons
# make master folder
transforms = drive_ls(path = folder, type = "folder", team_drive = hemibrain)
transforms = subset(transforms, name == "transforms")
if(!nrow(transforms)){
  drive_mkdir(name = "transforms",
              path = folder,
              overwrite = TRUE)
  transforms = drive_ls(path = folder, type = "folder", team_drive = hemibrain)
  transforms = subset(transforms, name == "transforms")
}
## and sub folder for each transform
trans = c("JRCFIB2018F","FAFB","JFRC2","JRC2018F","FCWB")
for(t in trans){
  transform = drive_ls(path = transforms, type = "folder", team_drive = hemibrain)
  transform = subset(transform, name == t)
  if(!nrow(transform)){
    drive_mkdir(name = t,
                path = transforms,
                overwrite = TRUE)
    transform = drive_ls(path = transforms, type = "folder", team_drive = hemibrain)
    transform = subset(transform, name == t)
    drive_mkdir(name = "data",
                path = transform,
                overwrite = TRUE)
    t.folder.data = drive_ls(path = transform, type = "folder", team_drive = hemibrain)
    t.folder.data = subset(t.folder.data, name == "data")
  }
  t.fh = paste0("/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/",t)
  t.list = list.files(t.fh,full.names = TRUE)
  t.rds = t.list[grepl(".rds",t.list)]
  for (tr in t.rds){
    drive_put(media = tr,
              path = transform[1,],
              overwrite = TRUE,
              verbose = TRUE)
  }
  t.neuron.fh.data.files = list.files(paste0(t.fh,"/data"), full.names = TRUE)
  error.files = c()
  for(t.neuron.fh.data.file in t.neuron.fh.data.files){
    e = tryCatch(drive_upload(media = t.neuron.fh.data.file,
                              path = t.folder.data[1,],
                              overwrite = FALSE,
                              verbose = TRUE),
                 error = function(e) t.neuron.fh.data.file)
    if(is.null(e)){
      error.files = c(error.files,t.neuron.fh.data.file)
    }
  }
}


