# set options
.onAttach <- function(libname, pkgname){
  options(hemibrain_data = paste0(getwd(),"/data-raw/hemibrainr_data/"))
  if(file.exists("/Volumes/GoogleDrive/Shared\ drives/hemibrain/")){
    options(Gdrive_hemibrain_data = "/Volumes/GoogleDrive/Shared\ drives/hemibrain/")
    # Set google sheets of interest
    options(hemibrainr_matching_gsheet = "1OSlDtnR3B1LiB5cwI5x5Ql6LkZd8JOS5bBr-HTi0pOw")
    options(flywire_lineages_gsheet = "1spGSuhUX6Hhn-8HH0U_ArIWUuPpMBFNjIjeSSh_MFVY")
    options(flywire_flagged_gsheet = "1rzG1MuZYacM-vbW7100aK8HeA-BY6dWAVXQ7TB6E2cQ")
    options(hemibrainr_gsheets = c(options()$hemibrainr_matching_gsheet,
                                   options()$flywire_lineages_gsheet,
                                   options()$flywire_flagged_gsheet))
  }else{
    options(Gdrive_hemibrain_data = "/Volumes/GoogleDrive/Shared\ drives/hemibrainr/")
  }
}

