# Gooogle drive based utilities

# hidden
good_savedir <- function(local=FALSE, Verbose = TRUE){
  if(isFALSE(local)){
    if(is.null(options()$remote_connectome_data)){
      hemibrainr_set_drive()
    }
    savedir = options()$remote_connectome_data
    if(dir.exists(savedir)){
      if(Verbose) message("Reading from Google Team Drive: ", hemibrainr_team_drive())
    }else{
      stop("Could not find hemibrainr data folder. Do you have Google Filestream or rclone active? or try local = TRUE: ", savedir)
    }
  }else if(isTRUE(local)){
    savedir = options()$hemibrain_data
    if(Verbose) message("Looking for data in local folder: ", savedir)
  }else{
    savedir = local
    if(Verbose) message("Looking for data in local folder: ", savedir)
  }
  savedir
}

# google upload neuronlistfh
## x is a neuronlist or a path to folder with saved neuronlistfh
googledrive_upload_neuronlistfh <- function(x,
                                            clean = FALSE,
                                            team_drive = hemibrainr_team_drive(),
                                            file_name = "neurons.rds",
                                            folder = "flywire_neurons",
                                            subfolder = NULL,
                                            WriteObjects = c("missing","yes"),
                                            dbClass = c("RDS", "RDS2", "DB1", "ZIP"),
                                            numCores = 1){
  # don't exhaust rate limit
  WriteObjects = match.arg(WriteObjects)
  dbClass = match.arg(dbClass)
  if(grepl("\\.zip$",file_name)){
    dbClass = "ZIP"
  }
  numCores = ifelse(numCores>10,10,numCores)
  if(dbClass%in%c("DB1","ZIP")){
    clean = FALSE
  }

  # Get drive
  td = googledrive::team_drive_get(team_drive)
  drive_td = googledrive::drive_find(type = "folder", team_drive = td)
  gfolder= subset(drive_td,drive_td$name==folder)[1,]
  if(!is.null(subfolder)){
    sub = googledrive::drive_ls(path = gfolder, type = "folder", team_drive = td)
    gsubfolder = subset(sub, sub$name ==subfolder )[1,]
    if(sum(is.na(gsubfolder$id))==length(gsubfolder$id)){
      googledrive::drive_mkdir(name = subfolder,
                               path = gfolder,
                               overwrite = TRUE)
      sub = googledrive::drive_ls(path = gfolder, type = "folder", team_drive = td)
      gsubfolder = subset(sub, sub$name ==subfolder )[1,]
    }
    gfolder = gsubfolder
  }
  sub = googledrive::drive_ls(path = gfolder, team_drive = td)

  # Save locally
  temp.data.csv = "no_data"
  if(nat::is.neuronlist(x)){
    temp = tempfile()
    if(!dbClass%in%c("DB1","ZIP")){
      temp.data.csv = file.path(temp,"data.csv")
      temp.data = file.path(temp,"data")
      dir.create(temp.data)
      on.exit(unlink(temp, recursive=TRUE))
    }else{
      temp.data = NULL
    }
    temp.nl = paste0(temp,"/",file_name)
    if(dbClass=="ZIP"){
      nat::write.neurons(nl, dir= temp.nl, format='qs', include.data.frame = TRUE, Force = TRUE)
    }else{
      nl = nat::as.neuronlistfh(x, dbdir = temp.data, dbClass = dbClass, WriteObjects = "yes")
      nat::write.neuronlistfh(nl, file = temp.nl, overwrite=TRUE)
    }
    local.path = NULL
  }else{
    if(!dbClass%in%c("DB1","ZIP")){
      temp.data = file.path(x,"data")
      temp.data.csv = file.path(x,"/data.csv")
      message("data folder: ", temp.data)
    }else{
      temp.data = NULL
    }
    temp.nl = list.files(x, pattern = ifelse(dbClass=="ZIP",".zip$",".rds$"), full.names =TRUE)
    message(".rds files: ", paste(temp.nl,collapse=", "))
    local.path = x
  }

  # Get data folder
  if(!dbClass%in%c("DB1","ZIP")){
    t.folder.data = googledrive::drive_ls(path = gfolder, type = "folder", team_drive = td)
    t.folder.data = subset(t.folder.data, t.folder.data$name == "data")[1,]
    if(is.na(t.folder.data$name)){
      googledrive::drive_mkdir(name = "data",
                               path = gfolder,
                               overwrite = TRUE)
      t.folder.data = googledrive::drive_ls(path = gfolder, type = "folder", team_drive = td)
      t.folder.data = subset(t.folder.data, t.folder.data$name == "data")[1,]
    }
    t.folder.csv= googledrive::drive_ls(path = gfolder, type = "csv", team_drive = td)
    t.folder.csv = subset(t.folder.csv, t.folder.csv$name == "data.csv")[1,]
    if(is.na(t.folder.csv)[1]){
      data.csv = data.frame()
    }else{
      googledrive::drive_download(file = t.folder.csv, path = temp.data.csv, overwrite = TRUE, verbose = FALSE)
      data.csv = tryCatch(readr::read_csv(temp.data.csv), error = function(e){message(e); data.frame()})
    }
  }

  # Upload
  if(!dbClass%in%c("DB1","ZIP")){
    t.list.master = list.files(temp.data,full.names = TRUE)
    error.files = upload = c()
  }else{
    t.list.master = NULL
  }
  if(WriteObjects=="missing" & !dbClass%in%c("DB1","ZIP")){
    if(nrow(data.csv)){
      t.list.master = t.list.master[! basename(t.list.master) %in% data.csv$data]
    }else{
      sub.data = googledrive::drive_ls(path = t.folder.data, team_drive = td, n_max = 50000)
      t.list.master = t.list.master[! basename(t.list.master) %in% sub.data$name]
    }
  }
  retry = 1
  while(length(t.list.master)&retry<=5){
    if(numCores>1){
      batch = 1
      batches = split(t.list.master, round(seq(from = 1, to = numCores, length.out = length(t.list.master))))
      foreach.upload <- foreach::foreach (batch = 1:length(batches)) %dopar% {
        t.list = batches[[batch]]
        for(t.neuron.fh.data.file in t.list){
          t = basename(t.neuron.fh.data.file)
          save.data =  t.folder.data
          upload = tryCatch({gsheet_manipulation(googledrive::drive_upload,
                                                 media = t.neuron.fh.data.file,
                                                 path = save.data,
                                                 verbose = FALSE)},
                            error = function(e){
                              message(e)
                              error.files <<- c(error.files,t.neuron.fh.data.file)
                              try(Sys.sleep(5))
                              NA
                            } )
        }
      }
    }else{
      pb <- progress::progress_bar$new(
        format = "  uploading [:bar] :current/:total eta: :eta",
        total = length(t.list.master), clear = FALSE, show_after = 1)
      for(t.neuron.fh.data.file in t.list.master){
        pb$tick()
        t = basename(t.neuron.fh.data.file)
        save.data =  t.folder.data
        upload = tryCatch({gsheet_manipulation(googledrive::drive_upload,
                                               media = t.neuron.fh.data.file,
                                               path = save.data,
                                               verbose = FALSE)},
                          error = function(e){
                            message(e)
                            error.files <<- c(error.files,t.neuron.fh.data.file)
                            try(Sys.sleep(5))
                            NA
                          } )
      }
    }
    success = basename(setdiff(t.list.master,error.files))
    if(nrow(data.csv)){
      data.csv = rbind(data.csv, data.frame(data = success))
    }else{
      data.csv =  data.frame(data = success)
    }
    if(length(error.files)){
      if(retry==5){
        message("Failed to upload: ", length(error.files), ' retries exhausted')
      }else{
        message("Failed to upload: ", length(error.files)," files, retry: ", retry)
        t.list.master = error.files
      }
      retry = retry + 1
      try(Sys.sleep(60))
    }
    readr::write_csv(x = data.csv, file = temp.data.csv)
  }

  # upload neuronlist
  for(nl.save in temp.nl){
    file_name = basename(nl.save)
    gfile = subset(sub, sub$name == file_name )
    if(nrow(gfile)){
      save.position = googledrive::as_id(gfile[1,]$id)
    }else{
      save.position = gfolder
    }
    google_drive_place(media = nl.save,
                       path = save.position,
                       verbose = TRUE)
    # upload db
    if(dbClass=="DB1"){
      nlfh = nat::read.neuronlistfh(nl.save)
      datafile = attributes(nlfh)$db@datafile
      file_name = basename(datafile)
      gfile = subset(sub, sub$name == file_name )
      if(nrow(gfile)){
        save.position = googledrive::as_id(gfile[1,]$id)
      }else{
        save.position = gfolder
      }
      google_drive_place(media = datafile,
                         path = save.position,
                         verbose = TRUE)
    }else{
      if(file.exists(temp.data.csv)){
        google_drive_place(media = temp.data.csv,
                           path = gfolder,
                           verbose = TRUE)
      }
    }
  }

  # Clean
  if(clean & !dbClass%in%c("DB1","ZIP")){
    message("Deleting unused file hash files for neurons ...")
    googledrive_clean_neuronlistfh(local.path = local.path,
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

# upload function
remote_upload <- function(file,
                          path,
                          remote = c("gdrive","dropbox"),
                          check = TRUE,
                          verbose = FALSE){

  # Upload
  try({
    if(remote=="gdrive"){
      if(check){
      }else{
        googledrive::drive_put(media = file,
                               path = path,
                               verbose = verbose)
      }
    }else if(remote=="dropbox"){
      if(!requireNamespace('rdrop2', quietly = T))
        stop("Please install suggested package rdrop2 to use Dropbox")
      rdrop2::drop_upload(file = file,
                          path = path,
                          mode = "overwrite",
                          mute = TRUE,
                          verbose = verbose)
    }
  }, silent = TRUE)

  # return
  invisible()
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

  # Get csv of data entries
  temp = tempdir(check = TRUE)
  temp.data.csv = paste0(temp,"/data.csv")
  on.exit(unlink(temp, recursive=TRUE))

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
    t.folder.csv= googledrive::drive_ls(path = gfolder, type = "csv", team_drive = td)
    t.folder.csv = subset(t.folder.csv, t.folder.csv$name == "data.csv")[1,]
    if(is.na(t.folder.csv)){
      data.csv = data.frame()
    }else{
      if(file.exists(temp.data.csv)){
        googledrive::drive_download(file = t.folder.csv, path = temp.data.csv,  overwrite = TRUE, verbose = FALSE)
        data.csv = readr::read_csv(temp.data.csv)
      }else{
        data.csv = data.frame()
      }
    }
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
      if(nrow(data.csv)){
        remove = try(setdiff(data.csv$data, all.keys))
      }
      if(!nrow(data.csv)||inherits("try-error",remove)){
        data = googledrive::drive_ls(path = googledrive::as_id(fdata$id), team_drive = team_drive)
        data = data[!grepl("\\.",data$name),]
        delete = data[!data$name%in%all.keys,]
        remove = rbind(remove, delete)
        remove = remove[!duplicated(remove$id),]
        message("cleaning out ", nrow(remove), " files")
      }else{
        message("cleaning out ", length(remove), " files")
      }
    }
    remove
  }

  # Find files to delete
  message("searching for files we might want to remove ...")
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
  if(nrow(data.csv)){
    if(length(remove)){
      fail = c()
      pb <- progress::progress_bar$new(
        format = "  deleting [:bar] :current/:total eta: :eta",
        total = length(remove$id), clear = FALSE, show_after = 1)
      for(r in remove){
        pb$tick()
        e = tryCatch(googledrive::drive_rm(r,  verbose = FALSE),
                     error = function(e){
                       fail <<- c(fail, r)
                       message(e)
                     } )
      }
      data.csv = subset(data.csv, !data.csv$data%in%setdiff(remove,fail))
      readr::write_csv(x = data.csv, file = temp.data.csv)
      google_drive_place(media = temp.data.csv,
                         path = gfolder,
                         verbose = TRUE)
    }
  }
  if(nrow(remove)){
    pb <- progress::progress_bar$new(
        format = "  deleting [:bar] :current/:total eta: :eta",
        total = length(remove$id), clear = FALSE, show_after = 1)
    for(r in remove$id){
        pb$tick()
        e = tryCatch(googledrive::drive_rm(googledrive::as_id(r), verbose = FALSE),
                     error = function(e){
                       message(e)
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

# hidden
gsheet_manipulation <- function(FUN,
                                ...,
                                wait = 10,
                                max.tries = 25,
                                return = FALSE,
                                Verbose = TRUE){
  sleep = wait
  success = FALSE
  tries = 0
  while(!success){
    if(Verbose){
      g = tryCatch(FUN(...),
                   error = function(e){
                     return({message(e);NULL})})
    }else{
      g = tryCatch(suppressMessages(FUN(...)),
                   error = function(e){
                     return({warning(e);NULL})})
    }
    if(!is.null(g)){
      success = TRUE
    }else{
      tries = tries+1
      if(tries>max.tries){
        stop("Could not make google manipulation")
      }
      if(Verbose){
        message("Google read/write failure(s), re-trying in ", sleep," seconds ...")
      }
      try(Sys.sleep(sleep))
      sleep = sleep + wait
      if(sleep > 600){
        slep <- 600
      }
    }
  }
  try(Sys.sleep(1))
  if(return){
    if(is.data.frame(g)){
      g = try(unlist_df(g), silent = FALSE)
    }
    return(g)
  }
}

# hidden
paste_coords <- function(xyz, sep = ", ", brackets = TRUE){
  paste0(ifelse(brackets,"(",NULL),paste(xyz,sep=sep,collapse=sep),ifelse(brackets,")",NULL))
}

# write columns to gsheet
gsheet_update_cols <- function(write.cols,
                               gs,
                               selected_sheet,
                               sheet = NULL,
                               Verbose = TRUE,
                               ...
                       ){
  shared = setdiff(write.cols,colnames(write.cols))
  if(!length(shared)){
    stop("write.cols missing from given data frame, gs")
  }
  if(!nrow(gs)){
    stop("No rows in given data frame, gs")
  }
  for(wc in write.cols){
    if(is.null(gs[[wc]])){
      gs[[wc]] = NA
    }
  }
  if(Verbose){
    pb = progress::progress_bar$new(
      format = "  updated :what [:bar] :percent eta: :eta",
      clear = FALSE, total = length(write.cols))
  }
  for(wc in write.cols){
    if(Verbose){ pb$tick(tokens = list(what = wc))}
    index = match(wc,colnames(gs))
    if(is.na(index)){
      stop("Column ",wc," could not be matched to sheet with columns: ", paste(colnames(gs),collapse=","))
    }
    if(index>26){
      letter = paste0("A",LETTERS[match(wc,colnames(gs))-26])
    }else if(index>52){
      letter = paste0("B",LETTERS[match(wc,colnames(gs))-52])
    }else{
      letter = LETTERS[match(wc,colnames(gs))]
    }
    range = paste0(letter,1,":",letter,nrow(gs)+1)
    update = tryCatch(as.data.frame(nullToNA(gs[,wc]), stringsAsFactors = FALSE),
                      error = function(e) data.frame(rep(NA,nrow(gs)), stringsAsFactors = FALSE))
    colnames(update) = wc
    gsheet_manipulation(FUN = googlesheets4::range_write,
                                     ss = selected_sheet,
                                     range = range,
                                     data = update,
                                     sheet = sheet,
                                     col_names = TRUE,
                                     Verbose=FALSE,
                                     ...)
  }
}

# Reorder sheets and remove duplicates
gsheet_reorder <- function(gs, tab, selected_sheet, remove.duplicates = TRUE, field = "root_id", Verbose = TRUE){
  if(remove.duplicates){
    gs = gs[!duplicated(gs),]
    if(!field%in%colnames(gs)){
      warning(field, " must be a valid column")
    }else{
      gs[[field]] = correct_id(gs[[field]])
      keep = !duplicated(gs[[field]]) | is.na(gs[[field]]) | gs[[field]]==0 | gs[[field]]== "(NA,NA,NA)"
      gs = gs[keep,]
    }
  }
  for(column in c(field, "status","whimsy","connectivity.type","cell_type","type","cbf","cell_body_fiber","ito_lee_hemilineage","User","added_by","tracer","last_updated_by","weight","count","score","NBLAST","nblast")){
    if(column %in% colnames(gs)){
      gs = gs[order(gs[[column]], decreasing = TRUE),]
    }
  }
  if(!nrow(gs)||!ncol(gs)){
    warning("Nothing to re-write for ",tab," of ",selected_sheet)
    return(NULL)
  }
  gsheet_manipulation(FUN = googlesheets4::write_sheet,
                      data = gs[0,],
                      ss = selected_sheet,
                      sheet = tab,
                      Verbose = Verbose)
  batches = split(1:nrow(gs), ceiling(seq_along(1:nrow(gs))/500))
  if(Verbose){
    pb = progress::progress_bar$new(
      format = "  appending rows :what [:bar] :percent eta: :eta",
      clear = FALSE, total = length(batches))
  }
  for(i in batches){
    if(Verbose) { pb$tick(tokens = list(what = paste0(min(i),":",max(i)))) }
    gsheet_manipulation(FUN = googlesheets4::sheet_append,
                        data = gs[min(i):max(i),],
                        ss = selected_sheet,
                        sheet = tab,
                        Verbose = FALSE)
  }
}



###### Seatable functionality #####






