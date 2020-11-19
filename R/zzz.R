# set options
.onAttach <- function(libname, pkgname){
  options(hemibrain_data = paste0(getwd(),"/data-raw/hemibrainr_data/"))

  paths=c(
    getOption("Gdrive_hemibrain_data"),
    file.path("/Volumes/GoogleDrive/Shared\ drives",
    c("hemibrain", "hemibrainr")))
  pathsok=file.exists(paths)

  if(any(pathsok)){
    options(Gdrive_hemibrain_data = which(pathsok)[1])
  }else{
    # fallback
    options(Gdrive_hemibrain_data = paste0(getwd(),"/data-raw/hemibrainr_data/"))

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
  https://flyconnectome.github.io/hemibrainr/articles/google_filestream.html
  or ?hemibrainr_set_drive
  ####################### remote data ########################
  ")
  }

  # Set Google sheets of interest
  options(hemibrainr_matching_gsheet = "1OSlDtnR3B1LiB5cwI5x5Ql6LkZd8JOS5bBr-HTi0pOw")
  options(flywire_lineages_gsheet = "1spGSuhUX6Hhn-8HH0U_ArIWUuPpMBFNjIjeSSh_MFVY")
  options(flywire_flagged_gsheet = "1rzG1MuZYacM-vbW7100aK8HeA-BY6dWAVXQ7TB6E2cQ")
  options(hemibrainr_gsheets = c(options()$hemibrainr_matching_gsheet,
                                 options()$flywire_lineages_gsheet,
                                 options()$flywire_flagged_gsheet))
}

