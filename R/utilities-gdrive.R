# Gooogle drive based utilities

# google upload neuronlistfh
## x is a neuronlist of a path to folder with saved neuronlistfh
googledrive_upload_neuronlistfh <- function(x,
                                            clean = FALSE,
                                            team_drive = hemibrainr_team_drive(),
                                            file_name = "neurons.rds",
                                            folder = "flywire_neurons",
                                            subfolder = NULL,
                                            numCores = 1){
  # don't exhaust rate limit
  numCores = ifelse(numCores>10,10,numCores)

  # Get drive
  td = googledrive::team_drive_get(team_drive)
  drive_td = googledrive::drive_find(type = "folder", team_drive = td)
  gfolder= subset(drive_td,drive_td$name==folder)[1,]
  if(!is.null(subfolder)){
    sub = googledrive::drive_ls(path = gfolder, type = "folder", team_drive = td)
    gfolder = subset(sub, sub$name ==subfolder )[1,]
  }
  sub = googledrive::drive_ls(path = gfolder, team_drive = td)

  # Get data folder
  t.folder.data = googledrive::drive_ls(path = gfolder, type = "folder", team_drive = td)
  t.folder.data = subset(t.folder.data, t.folder.data$name == "data")[1,]
  if(is.na(t.folder.data$name)){
    googledrive::drive_mkdir(name = "data",
                             path = gfolder,
                             overwrite = TRUE)
    t.folder.data = googledrive::drive_ls(path = gfolder, type = "folder", team_drive = td)
    t.folder.data = subset(t.folder.data, t.folder.data$name == "data")[1,]
  }

  # Save locally
  if(nat::is.neuronlist(x)){
    temp = tempfile()
    temp.data = paste0(temp,"/data")
    dir.create(temp.data)
    on.exit(unlink(temp, recursive=TRUE))
    temp.rds = paste0(temp,"/",file_name)
    nl = nat::as.neuronlistfh(x, dbdir= temp.data, WriteObjects = "yes")
    nat::write.neuronlistfh(nl, file= temp.rds, overwrite=TRUE)
    local.path = NULL
  }else{
    temp.data = paste0(x,"data")
    temp.rds = list.files(x, pattern = ".rds$", full.name =TRUE)
    message("data folder: ", temp.data)
    message(".rds files: ", paste(temp.rds,collapse=", "))
    local.path = x
  }

  # upload
  t.list.master = list.files(temp.data,full.names = TRUE)
  error.files = upload = c()
  sub.data = googledrive::drive_ls(path = t.folder.data, team_drive = td)
  if(numCores>1){
    batch = 1
    batches = split(t.list.master, round(seq(from = 1, to = numCores, length.out = length(t.list.master))))
    foreach.upload <- foreach::foreach (batch = 1:length(batches)) %dopar% {
      t.list = batches[[batch]]
      for(t.neuron.fh.data.file in t.list){
        t = basename(t.neuron.fh.data.file)
        if(t%in%sub.data$name){
          save.data = googledrive::as_id(subset(sub.data, sub.data$name==t)[1,]$id)
          # If exists, do not replace
        }else{
          save.data =  t.folder.data
          upload = tryCatch(google_drive_place(media = t.neuron.fh.data.file,
                                               path = save.data,
                                               verbose = FALSE,
                                               check = FALSE),
                            error = function(e){
                              message(e)
                              NA
                            } )
          if(is.na(upload)){
            error.files = c(error.files,t.neuron.fh.data.file)
          }
        }
      }
    }
  }else{
    pb <- progress::progress_bar$new(
      format = "  uploading [:bar] :current/:total eta: :eta",
      total = length(t.list.master), clear = FALSE, show_after = 1)
    for(t.neuron.fh.data.file in t.list.master){
      pb$tick()
      t = basename(t.neuron.fh.data.file)
      if(t%in%sub.data$name){
        save.data = googledrive::as_id(subset(sub.data, sub.data$name==t)[1,]$id)
        # If exists, do not replace
      }else{
        save.data =  t.folder.data
        upload = tryCatch(google_drive_place(media = t.neuron.fh.data.file,
                                             path = save.data,
                                             verbose = FALSE,
                                             check = FALSE),
                          error = function(e){
                            message(e)
                            NA
                          } )
        if(is.na(upload)){
          error.files = c(error.files,t.neuron.fh.data.file)
        }
      }
    }
  }
  if(length(error.files)){
    warning("Failed to upload: ", length(error.files)," files")
  }

  # upload .rds
  for(rds in temp.rds){
    file_name = basename(rds)
    gfile = subset(sub, sub$name == file_name )
    if(nrow(gfile)){
      save.position = googledrive::as_id(gfile[1,]$id)
    }else{
      save.position = gfolder
    }
    google_drive_place(media = rds,
                       path = save.position,
                       verbose = TRUE)
  }

  # Clean
  if(clean){
    message("Deleting unused file hash files for neurons ...")
    hemibrainr:::googledrive_clean_neuronlistfh(local.path = local.path,
                                   team_drive = hemibrainr_team_drive(),
                                   folder = folder,
                                   subfolder = subfolder)
  }
}

# Google drive
google_drive_place <- function(media,
                               path,
                               verbose = TRUE,
                               check = TRUE,
                               ...){
  # If path is folder, check contents to duplicates
  if(check){
    f = tryCatch(googledrive::is_folder(path), error = function(e) FALSE)
    if(f){
      ls = googledrive::drive_ls(path, ...)
      p = subset(ls, ls$name == basename(media))$id
      if(length(p)){
        path = googledrive::as_id(p[1])
      }
    }else{
      path = googledrive::as_id(path)
    }
  }

  if("drive_id"%in%class(path)){
    googledrive::drive_update(media = media,
                              file = path,
                              verbose = verbose,
                              ...)
  }else{
    googledrive::drive_put(media = media,
                           path = path,
                           verbose = verbose,
                           ...)
  }
}

# Remove unused filehash files
googledrive_clean_neuronlistfh <- function(team_drive = hemibrainr_team_drive(),
                                           local.path = NULL,
                                           folder = NULL,
                                           subfolder = NULL){
  # don't exhaust rate limit
  numCores = ifelse(numCores>10,10,numCores)
  if(!is.null(local.path)&is.null(subfolder)){
    stop("Either specify a subfolder on the ", team_drive, ' drive to clean,
         which matches the given local path, or set loca.path to NULL
         and attemot to clean the whole drive/folder on drive')
  }
  if(!is.null(local.path)&is.null(folder)&is.null(subfolder)){
    stop("A folder (e.g. flywire_neurons) and a
         subfolder (for brainspace, e.g. FCWB) must be given
         that match your local.path: ", local.path)
  }

  # Save .rds files locally, and find excess filehash
  find_excess <- function(sub,
                          local.path = NULL,
                          team_drive = NULL){
    temp = tempfile()
    temp.data = paste0(temp,"/RDS/")
    dir.create(temp.data, recursive = TRUE)
    on.exit(unlink(temp, recursive=TRUE))
    remove = googledrive::as_dribble()
    if(sum(grepl("rds$",sub$name))>1){
      message("Cleaning ", subset(drive_td,drive_td$id==folder)$name)
      fsub = subset(sub, grepl("rds$",sub$name))
      fdata = subset(sub, grepl("data$",sub$name))[1,]
      all.keys = c()
      if(is.null(local.path)){
        for(file in fsub$id){
          fnam = subset(fsub, fsub$id==file)$name
          path = paste0(temp.data,fnam)
          googledrive::drive_download(file = googledrive::as_id(file),
                                      path = path,
                                      overwrite = TRUE,
                                      verbose = FALSE)
          a = readRDS(path)
          b = attributes(a)
          keys = b$keyfilemap
          all.keys = c(all.keys,keys)
        }
      }else{
        files = list.files(local.path, pattern  = ".rds$",  full.names = TRUE)
        for(file in files){
          a = readRDS(file)
          b = attributes(a)
          keys = b$keyfilemap
          all.keys = c(all.keys,keys)
        }
      }
      data = googledrive::drive_ls(path = googledrive::as_id(fdata$id), team_drive = team_drive)
      data = data[!grepl("\\.",data$name),]
      delete = data[!data$name%in%all.keys,]
      remove = rbind(remove, delete)
    }
    remove = remove[!duplicated(remove$id),]
    remove
  }

  # Get drive
  td = googledrive::team_drive_get(team_drive)
  drive_td = googledrive::drive_find(type = "folder", team_drive = td)
  if(!is.null(folder)){
    gfolder= subset(drive_td,drive_td$name==folder)[1,]
    if(!is.null(subfolder)){
      sub = googledrive::drive_ls(path = gfolder, type = "folder", team_drive = td)
      gfolder = subset(sub, sub$name ==subfolder )[1,]
    }
    drive_td = googledrive::drive_ls(path = gfolder, team_drive = td)
  }

  # Find files to delete
  message("Searching for files we might want to remove ...")
  if(is.null(local.path)&is.null(subfolder)){
    # Iterate to find .rds files
    remove = googledrive::as_dribble()
    for(folder in drive_td$id){
      sub = googledrive::drive_ls(path = googledrive::as_id(folder), type = "folder", team_drive = td)
      if("data" %in% sub$name){
        sub = googledrive::drive_ls(path = googledrive::as_id(folder), team_drive = td)
        remove = rbind(remove, find_excess(sub,
                                           local.path = NULL,
                                           team_drive = td))
      }
    }
  }else{
    remove = find_excess(drive_td,
                         local.path = local.path,
                         team_drive = td)
  }

  # rm
  if(nrow(remove)){
    pb <- progress::progress_bar$new(
        format = "  deleting [:bar] :current/:total eta: :eta",
        total = length(rid), clear = FALSE, show_after = 1)
    for(r in remove$id){
        pb$tick()
        e = tryCatch(googledrive::drive_rm(googledrive::as_id(r), verbose = FALSE),
                     error = function(e){
                       message(e)
                       NULL
                     } )
    }
  }
}

# upload nblast matrix
googledrive_upload_nblast<- function(x,
                                     team_drive = hemibrainr_team_drive(),
                                     folder = "hemibrain_nblast",
                                     subfolder = NULL,
                                     file = NULL,
                                     threshold=-0.5,
                                     digits=3,
                                     format=c("rda", "rds"),
                                     ...){
  format=match.arg(format)

  # Get drive
  td = googledrive::team_drive_get(team_drive)
  drive_td = googledrive::drive_find(type = "folder", team_drive = td)
  gfolder= subset(drive_td,drive_td$name==folder)[1,]
  if(!is.null(subfolder)){
    sub = googledrive::drive_ls(path = gfolder, type = "folder", team_drive = td)
    gfolder = subset(sub, sub$name ==subfolder )[1,]
  }

  # Compress
  objname=deparse(substitute(x))
  colnames(x) = gsub("_m$","",colnames(x)) # factor in mirrored neurons
  rownames(x) = gsub("_m$","",rownames(x))
  x = tryCatch(apply(x, 2, function(i) tapply(x, rownames(x), sum, na.rm = TRUE)), error = function(e) x)
  x = tryCatch(t(apply(t(x), 2, function(i) tapply(i, colnames(x), sum, na.rm = TRUE))), error = function(e) x)
  x[x<threshold]=threshold
  x=round(x, digits=digits)

  # Save in temporary location
  newobjname <- paste0(objname, ".compressed")
  temp <- tempfile()
  on.exit(unlink(temp, recursive=TRUE))
  fname <- paste0(temp, "/",paste0(file, newobjname), ".", format)
  if(format=="rds") {
    saveRDS(x, file=fname, ...)
  } else {
    message(newobjname)
    assign(newobjname, x)
    save(list = newobjname, file = fname, ...)
  }

  # Does the file already exist?
  file_name = basename(fname)
  sub = googledrive::drive_ls(path = gfolder, team_drive = td)
  gfile = subset(sub, sub$name ==file_name)
  if(nrow(gfile)){
    save.position = googledrive::as_id(gfile[1,]$id)
  }else{
    save.position = gfolder
  }

  # Upload
  google_drive_place(media = fname,
                     path = save.position,
                     verbose = TRUE)
}
