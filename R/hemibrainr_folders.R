
# hidden
hemibrainr_folder_structure <- function(path = file.path(getwd(),"hemibrainr_data/")){
  # Make local folders
  brains = c("FAFB14", "FlyWire", "JRCFIB2018F","JRCFIB2018Fraw","JFRC2","JRC2018F","FCWB","JFRC2013")
  hemibrain.neurons = file.path(path,"hemibrain_neurons",brains,"data/")
  flywire.neurons = file.path(path,"flywire_neurons",brains,"data/")
  for(f in c(hemibrain.neurons,flywire.neurons)){
    dir.create(f, recursive = TRUE)
  }
  dir.create(paste0(path,"hemibrain_nblast", recursive = TRUE))
}

# hidden
# hemibrainr_gfolder_structure <- function(path){
#   brains = c("FAFB14", "FlyWire", "JRCFIB2018F","JRCFIB2018F","FAFB","JFRC2","JRC2018F","FCWB", "JFRC2013")
#   hemibrain.neurons = paste0(path,"hemibrain_neurons","/",brains,"/data/")
#   flywire.neurons = paste0(path,"flywire_neurons","/",brains,"/data/")
# }
