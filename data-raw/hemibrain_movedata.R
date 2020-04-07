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
drive_hemibrain = drive_find(type = "folder", team_drive = hemibrain)
hemibrain_neurons= subset(drive_hemibrain,name=="hemibrain_neurons")

# Save all neuron bodyids
drive_upload(media = "/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/hemibrain_all_neuron_bodyids.csv",
             path = hemibrain_neurons,
             overwrite = TRUE,
             verbose = TRUE)

# Save raw neurons
drive_upload(media = "/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/hemibrain_all_neurons.rda",
             path = hemibrain_neurons,
             overwrite = TRUE,
             verbose = TRUE)

# Save pre-processed neurons
drive_upload(media = "/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/hemibrain_all_neurons_flow_checked.rda",
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
drive_upload(media = split.neurons,
             path = folder,
             overwrite = TRUE,
             verbose = TRUE)
drive_upload(media = split.neuron.fh.neuronlist,
             path = folder,
             overwrite = TRUE,
             verbose = TRUE)
drive_upload(media = split.points,
             path = folder,
             overwrite = TRUE,
             verbose = TRUE)
drive_upload(media = split.synapses,
             path = folder,
             overwrite = TRUE,
             verbose = TRUE)
drive_upload(media = split.connections,
             path = folder,
             overwrite = TRUE,
             verbose = TRUE)
drive_upload(media = split.elist,
             path = folder,
             overwrite = TRUE,
             verbose = TRUE)
drive_upload(media = split.metrics,
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
               error = function(e) NULL)
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
               error = function(e) NULL)
  if(is.null(e)){
    error.files = c(error.files,split.neuron.fh.data.file)
  }
}
