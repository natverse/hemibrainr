# Use hemibrainr with rclone

#' @name hemibrainr_set_drive
#' @export
hemibrainr_rclone <- function(drive = "hemibrainr", path=file.path(getwd(),"hemibrainr_data/")){
  os = get_os()
  if(os=="osx"){
    if(!dir.exists(path)){
      dir.create(path)
    }
    empty=list.files(path)
    if(length(empty)){
      stop("path must lead to an empty directory.
    This directory will be 'replaced' by the mounted drive.
    You chose: ",
           path,"
    Perhaps try a new path, like: ",
           file.path(getwd(),"hemibrainr_rclone/"))
    }
    command = sprintf("rclone cmount %s: %s", drive, path)
  }else if (os == "windows"){
    command = sprintf("rclone mount %s: %s", drive, path)
  }else{
    if(!dir.exists(path)){
      dir.create(path)
    }
    empty=list.files(path)
    if(length(empty)){
      stop("path must lead to an empty directory.
    This directory will be 'replaced' by the mounted drive.
    You chose: ",
           path,"
    Perhaps try a new path, like: ",
           file.path(getwd(),"hemibrainr_rclone/"))
    }
    command = sprintf("rclone mount %s: %s", drive, path)
  }
  options(Gdrive_hemibrain_data = path)
  message("Google drive path set to: ", options()$Gdrive_hemibrain_data)
  system(command, wait = FALSE)
  try(hemibrainr_google_login(),silent=TRUE)
  invisible()
}

#' @name hemibrainr_set_drive
#' @export
hemibrainr_rclone_unmount <- function(path=file.path(getwd(),"hemibrainr_data/")){
  os = get_os()
  if(os=="osx"){
    # OS X
    system(sprintf("umount %s", path))
  }else{
    # Linux
    system(sprintf("fusermount -u %s", path))
  }
  options(Gdrive_hemibrain_data = "/Volumes/GoogleDrive/Shared\ drives/hemibrainr/")
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

# hidden
hemibrainr_google_login <- function(path = NULL){
  if(is.null(path)){
    path = file.path(options()$Gdrive_hemibrain_data,"annotations")
  }
  hemibrainr_service_account_key = file.path(path,"hemibrainr-571dd013f664.json")
  if(file.exists(hemibrainr_service_account_key)){
    googledrive::drive_auth(path = hemibrainr_service_account_key)
    googlesheets4::gs4_auth(
      scopes = 'https://www.googleapis.com/auth/spreadsheets',
      path = hemibrainr_service_account_key,
      use_oob = TRUE
    )
  }else{
    hemibrainr_ghseets_api_key = Sys.getenv("hemibrainr_ghseets_api_key")
    hemibrainr_clientid = Sys.getenv("hemibrainr_clientid")
    hemibrainr_secret = Sys.getenv("hemibrainr_secret")
    if(is.null(hemibrainr_ghseets_api_key)||hemibrainr_ghseets_api_key==""){
      hemibrainr_ghseets_api_key = options()$hemibrainr_ghseets_api_key
      hemibrainr_clientid = options()$hemibrainr_clientid
      hemibrainr_secret = options()$hemibrainr_secret
    }
    if(all(sapply(c(hemibrainr_ghseets_api_key,
                    hemibrainr_clientid,
                    hemibrainr_secret),is.null))){
      if(requireNamespace('httr', quietly = TRUE)) {
        message("We can access Google resources using the Google Cloud  Platform app 'hemibrainr'")
        # bring your own app via client id (aka key) and secret
        google_app <- httr::oauth_app(
          "hemibrainr",
          key = hemibrainr_clientid,
          secret = hemibrainr_secret
        )
        googlesheets4::gs4_auth_configure(app = google_app, api_key = hemibrainr_ghseets_api_key)
        googledrive::drive_auth_configure(app = google_app, api_key = hemibrainr_ghseets_api_key)
        # googlesheets4::gs4_oauth_app()
        # googlesheets4::gs4_api_key()
      }
    }
  }
}
