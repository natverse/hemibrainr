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
nblast.folder = subset(drive_hemibrain,name=="hemibrain_nblast")

# Save NBLAST files
nblast.files  = list.files(paste0("/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/NBLAST/",identifier),pattern = "rda", full.name = TRUE)
nblast.files  = c(nblast.files,list.files("/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/NBLAST/",pattern = "rda", full.name = TRUE))
for(nblast.file in nblast.files){
  hemibrainr:::google_drive_place(media = nblast.file,
               path = nblast.folder,
               verbose = TRUE)
}

# Save all neuron bodyids
hemibrainr:::google_drive_place(media = "/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/hemibrain_all_neuron_bodyids.csv",
             path = hemibrain_neurons,
             verbose = TRUE)

# Save raw neurons
hemibrainr:::google_drive_place(media = "/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/hemibrain_all_neurons.rda",
             path = hemibrain_neurons,
             verbose = TRUE)

# Save pre-processed neurons
hemibrainr:::google_drive_place(media = "/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/hemibrain_all_neurons_flow_checked.rda",
             path = hemibrain_neurons,
             verbose = TRUE)

# Split files to save
split.neurons = paste0("/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/hemibrain_all_neurons_flow_",identifier,".rda")
split.points = paste0("/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/hemibrain_all_neurons_splitpoints_",identifier,".csv")
split.synapses = paste0("/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/hemibrain_all_neurons_synapses_",identifier,".csv")
split.connections = paste0("/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/hemibrain_all_neurons_connections_",identifier,".csv")
split.elist = paste0("/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/hemibrain_all_neurons_edgelist_",identifier,".csv")
split.metrics = paste0("/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/hemibrain_all_neurons_metrics_",identifier,".csv")
hemibrainr:::google_drive_place(media = split.neurons,
             path = hemibrain_neurons,
             verbose = TRUE)
hemibrainr:::google_drive_place(media = split.points,
             path = hemibrain_neurons,
             verbose = TRUE)
hemibrainr:::google_drive_place(media = split.synapses,
             path = hemibrain_neurons,
             verbose = TRUE)
hemibrainr:::google_drive_place(media = split.connections,
             path = hemibrain_neurons,
             verbose = TRUE)
hemibrainr:::google_drive_place(media = split.elist,
             path = hemibrain_neurons,
             verbose = TRUE)
hemibrainr:::google_drive_place(media = split.metrics,
             path = hemibrain_neurons,
             verbose = TRUE)

## Save transformed neurons
trans = c("hemibrain_all_neurons_flow_polypre_centrifugal_synapses",
          "hemibrain_all_neurons_flow_polypre_centrifugal_synapses/dotprops",
          "JRCFIB2018F","FAFB14","JFRC2","JRC2018F","FCWB")
for(t in trans){
  if(grepl("hemibrain_all_neurons_flow_polypre_centrifugal_synapses",t)){
    gfolder = "JRCFIB2018Fraw"
  }else{
    gfolder = t
  }
  transform = drive_ls(path = hemibrain_neurons, type = "folder", team_drive = hemibrain)
  transform = subset(transform, name == gfolder)
  if(!nrow(transform)){
    drive_mkdir(name = gfolder,
                path = hemibrain_neurons,
                overwrite = TRUE)
    transform = drive_ls(path = hemibrain_neurons, type = "folder", team_drive = hemibrain)
    transform = subset(transform, name == gfolder)
    drive_mkdir(name = "data",
                path = transform,
                overwrite = TRUE)
    t.folder.data = drive_ls(path = transform, type = "folder", team_drive = hemibrain)
    t.folder.data = subset(t.folder.data, name == "data")
  }
  t.folder.data = drive_ls(path = transform, type = "folder", team_drive = hemibrain)
  t.folder.data = subset(t.folder.data, name == "data")
  if(!nrow(t.folder.data)){
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
    hemibrainr:::google_drive_place(media = tr,
              path = transform[1,],
              verbose = TRUE)
  }
  t.neuron.fh.data.files = list.files(paste0(t.fh,"/data"), full.names = TRUE)
  error.files = c()
  for(t.neuron.fh.data.file in t.neuron.fh.data.files){
    e = tryCatch(hemibrainr:::google_drive_place(media = t.neuron.fh.data.file,
                              path = t.folder.data[1,],
                              verbose = TRUE),
                 error = function(e) t.neuron.fh.data.file)
    if(is.null(e)){
      error.files = c(error.files,t.neuron.fh.data.file)
    }
  }
}


