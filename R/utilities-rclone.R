# Use hemibrainr with rclone

hemibrainr_rclone <- function(remote = "hemibrainr", path=paste0(getwd(),"hemibrainr_rclone/")){
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
    command = sprintf("rclone cmount %s: %s", remote, path)
  }else{
    command = sprintf("rclone mount %s: %", remote, path)
  }
  options(Gdrive_hemibrain_data = path)
  message("option Gdrive_hemibrain_data set to ", path)
  system(command, wait = FALSE)
}

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
