
# hidden
hemibrainr_folder_structure <- function(start = options()$hemibrain_data,
                                        team.drive = "hemibrain"){

  # Make local folders
  brains = c("FAFB14", "FlyWire", "JRCFIB2018F","JRCFIB2018Fraw","JFRC2","JRC2018F","FCWB","JFRC2013")
  hemibrain.neurons = paste0(start,"hemibrain_neurons","/",brains,"/data/")
  flywire.neurons = paste0(start,"flywire_neurons","/",brains,"/data/")
  for(f in c(hemibrain.neurons,flywire.neurons)){
    dir.create(f, recursive = TRUE)
  }
  dir.create(paste0(start,"hemibrain_nblast", recursive = TRUE))
}

# hidden
hemibrainr_gfolder_structure <- function(){
  brains = c("FAFB14", "FlyWire", "JRCFIB2018F","JRCFIB2018F","FAFB","JFRC2","JRC2018F","FCWB", "JFRC2013")
  hemibrain.neurons = paste0(start,"hemibrain_neurons","/",brains,"/data/")
  flywire.neurons = paste0(start,"flywire_neurons","/",brains,"/data/")




}
