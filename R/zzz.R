# set options
.onAttach <- function(libname, pkgname){
  options(hemibrain_data = paste0(getwd(),"/data-raw/hemibrain_data/"))
  options(Gdrive_hemibrain_data = "/Volumes/GoogleDrive/Shared\ drives/hemibrain/")
}

