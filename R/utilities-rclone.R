# Use hemibrainr with rclone

#' @name hemibrainr_set_drive
#' @export
hemibrainr_rclone <- function(Gdrive = "hemibrainr", path=paste0(getwd(),"hemibrainr_rclone/")){
  os = get_os()
  if(!dir.exists(path)){
    dir.create(path)
  }
  empty=list.files(path)
  if(length(empty)){
    stop("path must lead to an empty directory.
         This directory will be 'replaced' by the mounted drive.
         You chose: ",
         path)
  }
  if(os=="osx"){
    command = sprintf("rclone cmount %s: %s", Gdrive, path)
  }else{
    command = sprintf("rclone mount %s: %s", Gdrive, path)
  }
  options(Gdrive_hemibrain_data = path)
  message("Google drive path set to: ", options()$Gdrive_hemibrain_data)
  system(command, wait = FALSE)
}

#' @name hemibrainr_set_drive
#' @export
hemibrainr_rclone_unmount <- function(path=paste0(getwd(),"hemibrainr_rclone/")){
  os = get_os()
  if(os=="osx"){
    # OS X
    system(sprintf("umount %s", path))
  }else{
    # Linux
    system(sprintf("fusermount -u %s", path))
  }
  options(Gdrive_hemibrain_data = "/Volumes/GoogleDrive/Shared\ drives/hemibrain/")
}

# hidden
get_os <- function(){
  sysinf <- Sys.info()
  if (!is.null(sysinf)){
    os <- sysinf['sysname']
    if (os == 'Darwin')
      os <- "osx"
  } else { ## mystery machine
    os <- .Platform$OS.type
    if (grepl("^darwin", R.version$os))
      os <- "osx"
    if (grepl("linux-gnu", R.version$os))
      os <- "linux"
  }
  tolower(os)
}
