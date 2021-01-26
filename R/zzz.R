# set options
.onAttach <- function(libname, pkgname){
  options(hemibrain_data = file.path(getwd(),"hemibrainr_data"))

  paths=c(
    nullToNA(getOption("Gdrive_hemibrain_data")),
    file.path("/Volumes/GoogleDrive/Shared\ drives",
    c("hemibrainr", "hemibrain")))
  pathsok=file.exists(paths)

  if(any(pathsok)){
    options(Gdrive_hemibrain_data = paths[which(pathsok)[1]])
  }else{
    # fallback
    options(Gdrive_hemibrain_data = file.path(getwd(),"hemibrainr_data"))

    packageStartupMessage("
  ####################### remote data ########################
  To use data grabbing functions in hemibrainr please point
  hemibrainr to a shared drive, mounted with Google
  filestream OR rclone, by setting the option
      options(Gdrive_hemibrain_data='/path/to/drive')

  To mount with rclone, just run in R:
      hemibrainr_rclone()
  And unmount with:
      hemibrainr_rclone_unmount()
  You will need to have installed and configured rclone. See
  https://natverse.github.io/hemibrainr/articles/google_filestream.html
  or ?hemibrainr_set_drive
  ####################### remote data ########################
  ")
  }

  # Set googlesheets API key
  key = options()$hemibrainr_ghseets_api_key
  if(!is.null(key)){
    message("Setting goglesheets API key for hemibrainr_ghseets")
    googlesheets4::gs4_auth_configure(api_key = key)
    googlesheets4::gs4_deauth()
  }

  # Set Google sheets of interest
  hemibrainr_matching_gsheet = ifelse(!is.null(getOption("hemibrainr_matching_gsheet")),getOption("hemibrainr_matching_gsheet"),"1OSlDtnR3B1LiB5cwI5x5Ql6LkZd8JOS5bBr-HTi0pOw")
  flywire_lineages_gsheet = ifelse(!is.null(getOption("flywire_lineages_gsheet")),getOption("flywire_lineages_gsheet"),"1QyuHFdqz705OSxXNynD9moIsLvZGjjBjylx5sGZP2Yg")
  flywire_flagged_gsheet = ifelse(!is.null(getOption("flywire_flagged_gsheet")),getOption("flywire_flagged_gsheet"),"1rzG1MuZYacM-vbW7100aK8HeA-BY6dWAVXQ7TB6E2cQ")
  options(hemibrainr_matching_gsheet = hemibrainr_matching_gsheet)
  options(flywire_lineages_gsheet = flywire_lineages_gsheet)
  options(flywire_flagged_gsheet = flywire_flagged_gsheet)
  options(hemibrainr_gsheets = c(options()$hemibrainr_matching_gsheet,
                                 options()$flywire_lineages_gsheet,
                                 options()$flywire_flagged_gsheet))
}





