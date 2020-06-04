############## Script to split NBALST all neurons and their compartments  ###############
######## Originally run from JData5, the Jefferis lab drive on Max at the MRC LMB #######

# Set things ups
message("Start")
source("/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/R/startup/packages.R")
message("packages loaded")

# Get registrations
library(lhns)
library(nat.jrcbrains)
library(googledrive)
# download_saalfeldlab_registrations()

# Get functions
source("/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/R/startup/functions.R")
message("functions loaded")

# A folder to work in
lhfolder = paste0("/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/light_level/lhns/")
flycircuitfolder = paste0("/net/flystore3/jdata/jdata5/JPeople/Alex/FIBSEM/data/neurons/fibsem/light_level/flycircuit/")
dir.create(lhfolder, recursive = TRUE)
dir.create(flycircuitfolder, recursive = TRUE)

# Find hdf5r
Sys.setenv(LD_LIBRARY_PATH=paste0("/public/flybrain/hdf5/lib:", Sys.getenv('LD_LIBRARY_PATH')))

# Parallelise
numCores <- detectCores()
doMC::registerDoMC(numCores/8)
message("Using ", numCores/8, " cores")

# Get some stuff I have already put on google drive
# Google drive location:
hemibrain = team_drive_get("hemibrain")
drive_hemibrain = drive_find(type = "folder", team_drive = hemibrain)
light_level= subset(drive_hemibrain,name=="light_level")
lhns = subset(drive_hemibrain,name=="lhns")
flycircuit= subset(drive_hemibrain,name=="flycircuit")

# Download transformed most.lhns
mlhns = drive_ls(lhns)
for(r in 1:nrow(mlhns)){
  nam = as.character(mlhns[r,"name"])
  drive_download(file = mlhns[r,],
                 path = paste0(lhfolder,nam),
                 type = NULL,
                 overwrite = TRUE,
                 verbose = TRUE)
}

# Download flycircuit neurons
fcs = drive_ls(flycircuit)
for(r in 1:nrow(fcs)){
  nam = as.character(fcs[r,"name"])
  drive_download(file = fcs[r,],
                 path = paste0(flycircuitfolder,nam),
                 type = NULL,
                 overwrite = TRUE,
                 verbose = TRUE)
}

# Convert most.lhns
jlns = lhns::jfw.lhns
jlns[,] = lhns::jfw.lhns[,intersect(colnames(lhns::most.lhns),colnames(lhns::jfw.lhns))]
jlns.dps = lhns::jfw.lhns.dps
jlns.dps[,] = lhns::jfw.lhns[,intersect(colnames(lhns::most.lhns.dps),colnames(lhns::jfw.lhns.dps))]
most.lhns = c(lhns::most.lhns, jlns)
most.lhns.dps = c(lhns::most.lhns, jlns.dps)
most.lhns.hemi.dps = nat.templatebrains::xform_brain(most.lhns.dps, reference= "JRCFIB2018F", sample="FCWB", .parallel = TRUE, OmitFailures = TRUE)
most.lhns.hemi = nat.templatebrains::xform_brain(most.lhns, reference= "JRCFIB2018F", sample="FCWB", .parallel = TRUE, OmitFailures = TRUE)
saveRDS(most.lhns.hemi.dps, file = paste0(lhfolder,"most_lhns_dps_JRCFIB2018F.rds"))
saveRDS(most.lhns.hemi, file = paste0(lhfolder,"most_lhns_JRCFIB2018F.rds"))

# Convert most.lhins
most.lhins.hemi.dps = nat.templatebrains::xform_brain(lhns::most.lhins.dps, reference= "JRCFIB2018F", sample="FCWB", .parallel = TRUE, OmitFailures = TRUE)
most.lhins.hemi = nat.templatebrains::xform_brain(lhns::most.lhins, reference= "JRCFIB2018F", sample="FCWB", .parallel = TRUE, OmitFailures = TRUE)
saveRDS(most.lhins.hemi.dps, file = paste0(lhfolder,"most_lhins_dps_JRCFIB2018F.rds"))
saveRDS(most.lhins.hemi, file = paste0(lhfolder,"most_lhins_JRCFIB2018F.rds"))

# Convert flycircuit neurons
fcs = flycircuit_neurons(savedir = flycircuitfolder,
                         folder="",
                         local = TRUE,
                         cable = "all",
                         data = "neuronlist")
fcs.dps = flycircuit_neurons(savedir = flycircuitfolder,
                         folder="",
                         local = TRUE,
                         cable = "all",
                         data = "dps")

# Convert
fcs.hemi.dps = nat.templatebrains::xform_brain(fcs.dps, reference= "JRCFIB2018F", sample="FCWB", .parallel = TRUE, OmitFailures = TRUE)
fcs.hemi = nat.templatebrains::xform_brain(fcs, reference= "JRCFIB2018F", sample="FCWB", .parallel = TRUE, OmitFailures = TRUE)
saveRDS(fcs.hemi.dps, file = paste0(flycircuitfolder,"FlyCircuit_all_neurons_dps_JRCFIB2018F.rds"))
saveRDS(fcs.hemi, file = paste0(flycircuitfolder,"FlyCircuit_all_neurons_JRCFIB2018F.rds"))

# Save to google drive
## Save lh files
lh.files  = list.files(lhfolder,pattern = "rds", full.name = TRUE)
for(lh.file in lh.files){
  drive_upload(media = lh.file,
               path = lhns,
               overwrite = TRUE,
               verbose = TRUE)
}
fc.files  = list.files(flycircuitfolder,pattern = "rds", full.name = TRUE)
for(fc.file in fc.files){
  drive_upload(media = fc.file,
               path = flycircuit,
               overwrite = TRUE,
               verbose = TRUE)
}




