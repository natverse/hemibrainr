# set options
.onAttach <- function(libname, pkgname){
  options(hemibrain_data = paste0(getwd(),"/data-raw/hemibrain_data/"))
  if(file.exists("/Volumes/GoogleDrive/Shared\ drives/hemibrain/")){
    options(Gdrive_hemibrain_data = "/Volumes/GoogleDrive/Shared\ drives/hemibrain/")
  }else{
    options(Gdrive_hemibrain_data = "/Volumes/GoogleDrive/Shared\ drives/hemibrainr/")
  }
}

