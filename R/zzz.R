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
  }else if (file.exists("/Volumes/GoogleDrive/Shared\ drives/hemibrainr/")){
    options(Gdrive_hemibrain_data = "/Volumes/GoogleDrive/Shared\ drives/hemibrainr/")
    # Set google sheets of interest
    options(hemibrainr_matching_gsheet = "1OSlDtnR3B1LiB5cwI5x5Ql6LkZd8JOS5bBr-HTi0pOw")
    options(flywire_lineages_gsheet = "1spGSuhUX6Hhn-8HH0U_ArIWUuPpMBFNjIjeSSh_MFVY")
    options(flywire_flagged_gsheet = "1rzG1MuZYacM-vbW7100aK8HeA-BY6dWAVXQ7TB6E2cQ")
    options(hemibrainr_gsheets = c(options()$hemibrainr_matching_gsheet,
                                   options()$flywire_lineages_gsheet,
                                   options()$flywire_flagged_gsheet))
  }else{
    packageStartupMessage("
          ########################remote data#########################
          You do not have a shared Google drive mounted that has data
          for hemibrainr. To use data grabbing functions in hemibrainr
          please mount hemibrainr to your drive, mounted with Google
          filestream OR rclone, by setting the option
                            Gdrive_hemibrain_data
          To the right location. To mount with rclone, just run in R:
                            hemibrainr_rclone()
          And unmount with:
                            hemibrainr_rclone_unmount()
          You will need to have installed and configre rclone
          ########################remote data#########################
          ")
    # Set google sheets of interest
    options(Gdrive_hemibrain_data = paste0(getwd(),"/data-raw/hemibrainr_data/"))
    options(hemibrainr_matching_gsheet = "1OSlDtnR3B1LiB5cwI5x5Ql6LkZd8JOS5bBr-HTi0pOw")
    options(flywire_lineages_gsheet = "1spGSuhUX6Hhn-8HH0U_ArIWUuPpMBFNjIjeSSh_MFVY")
    options(flywire_flagged_gsheet = "1rzG1MuZYacM-vbW7100aK8HeA-BY6dWAVXQ7TB6E2cQ")
    options(hemibrainr_gsheets = c(options()$hemibrainr_matching_gsheet,
                                   options()$flywire_lineages_gsheet,
                                   options()$flywire_flagged_gsheet))
  }
}

