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
    temp.rds = list.files(x, pattern = ".rds$", full.names =TRUE)
    message("data folder: ", temp.data)
    message(".rds files: ", paste(temp.rds,collapse=", "))
    local.path = x
  }

  # upload
  t.list.master = list.files(temp.data,full.names = TRUE)
  error.files = upload = c()
  sub.data = googledrive::drive_ls(path = t.folder.data, team_drive = td)
  t.list.master = t.list.master[basename(t.list.master) %in% sub.data$name]
  if(numCores>1){
    batch = 1
    batches = split(t.list.master, round(seq(from = 1, to = numCores, length.out = length(t.list.master))))
    foreach.upload <- foreach::foreach (batch = 1:length(batches)) %dopar% {
      t.list = batches[[batch]]
      for(t.neuron.fh.data.file in t.list){
        t = basename(t.neuron.fh.data.file)
        save.data =  t.folder.data
        upload = tryCatch(gsheet_manipulation(googledrive::drive_upload,
                                              media = t.neuron.fh.data.file,
                                              path = save.data,
                                              verbose = FALSE),
                            error = function(e){
                            cat(as.character(e))
                            NA
                            } )
        if(is.na(upload)){
          error.files = c(error.files,t.neuron.fh.data.file)
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
      save.data =  t.folder.data
      upload = tryCatch(gsheet_manipulation(googledrive::drive_upload,
                                            media = t.neuron.fh.data.file,
                                            path = save.data,
                                            verbose = FALSE),
                          error = function(e){
                            cat(as.character(e))
                            NA
                          } )
      if(is.na(upload)){
        error.files = c(error.files,t.neuron.fh.data.file)
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
      message("Cleaning out ", nrow(remove), " files")
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
        total = length(remove$id), clear = FALSE, show_after = 1)
    for(r in remove$id){
        pb$tick()
        e = tryCatch(googledrive::drive_rm(googledrive::as_id(r), verbose = FALSE),
                     error = function(e){
                       cat(as.character(e))
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


# Obtains flywire IDs from flywire posirions on Google sheets
# The sheets must either have columns fw.x, fw.y, fw.z
# or one column, flywire.xyz with number separated by commas, in the form: x,y,z
flywire_ids_update <- function(selected_sheets = NULL,
                               chosen.columns = c("fw.x","fw.y",'fw.z', 'flywire.xyz',
                                                  "flywire.id", "skid",
                                                  "FAFB.xyz", "side",
                                                  "ItoLee_Hemilineage", "Hartenstein_Hemilineage",
                                                  "hemibrain_match"),
                               numCores = 1){
  # Read selected sheets and extract positions for flywire neurons
  # One xyz position is enough to identify a neuron
  # We do this because flywire.ids change all of the time
  if(is.null(selected_sheets)){
    selected_sheets = getOption("hemibrainr_gsheets", stop("Please set option('hemibrainr_gsheets')"))
  }
  gs = data.frame()
  for(selected_sheet in selected_sheets){
    ## Read Google sheets and extract glywire neuron positions
    tabs = gsheet_manipulation(FUN = googlesheets4::sheet_names,
                               ss = selected_sheet,
                               return = TRUE)
    for(tab in tabs){
      gs.t = gsheet_manipulation(FUN = googlesheets4::read_sheet,
                                 wait = 20,
                                 ss = selected_sheet,
                                 sheet = tab,
                                 guess_max = 3000,
                                 return = TRUE)
      used.cols = colnames(gs.t)
      if(ncol(gs.t)&&sum(grepl("fw.x|flywire.xyz",colnames(gs.t)))>0){
        # Separate x,y,z positions
        gs2 = subset(gs.t, !is.na(gs.t$flywire.xyz))
        gs1 = subset(gs.t[setdiff(rownames(gs.t),rownames(gs2)),], !is.na(gs.t$fw.x))
        gs0 = gs.t[setdiff(rownames(gs.t),c(rownames(gs1),rownames(gs2))),]
        if(nrow(gs1)){
          gs1$flywire.xyz = tryCatch(apply(gs1[,c("fw.x","fw.y",'fw.z')],1,paste,sep=";",collapse=";"),
                                     error = function(e) NA)
        }
        if(nrow(gs2)){
          positions.gs = sapply(gs2$flywire.xyz,strsplit,",|/|;")
          ruleofthree = sapply(positions.gs,function(p) length(p)==3)
          positions.gs = positions.gs[ruleofthree]
          positions.gs = do.call(rbind, positions.gs)
          positions.gs = apply(positions.gs,2,as.numeric)
          gs2[ruleofthree,c("fw.x","fw.y",'fw.z')] = positions.gs
        }
        gs.t = rbind(gs2,gs1)
        if(nrow(gs.t)){
          # Get flywire IDs from these positions
          bbx = matrix(c(5100, 1440, 16, 59200, 29600, 7062),ncol=3,byrow = TRUE)
          bbx = nat::boundingbox(scale(bbx, scale = 1/c(4, 4, 40), center = FALSE))
          if(numCores>1){
            batch = 1
            batches = split(1:nrow(gs.t), round(seq(from = 1, to = numCores, length.out = nrow(gs.t))))
            foreach.ids <- foreach::foreach (batch = 1:length(batches)) %dopar% {
              pos = gs.t[batches[[batch]],]
              pos = pos[apply(pos, 1, function(row) sum(is.na(row[c("fw.x","fw.y",'fw.z')]))==0),]
              p = nat::pointsinside(pos[,c("fw.x","fw.y",'fw.z')],bbx)
              pos = pos[p,]
              if(nrow(pos)){
                i <- tryCatch(fafbseg::flywire_xyz2id(pos[,c("fw.x","fw.y",'fw.z')], rawcoords = TRUE),
                              error = function(e){cat(as.character(e));rep("0",nrow(pos))})
                names(i) = pos$flywire.xyz
                i[is.na(i)|is.nan(i)] = 0
                i
              }
            }
          }else{
            pos = gs.t[apply(gs.t, 1, function(row) sum(is.na(row[c("fw.x","fw.y",'fw.z')]))==0),]
            p = nat::pointsinside(pos[,c("fw.x","fw.y",'fw.z')],bbx)
            pos = pos[p,]
            if(nrow(pos)){
              foreach.ids = fafbseg::flywire_xyz2id(pos[,c("fw.x","fw.y",'fw.z')], rawcoords = TRUE)
              names(foreach.ids) = pos[,"flywire.xyz"]
            }else{
              foreach.ids = NULL
            }
          }
          fids = unlist(foreach.ids)
          fids[is.na(fids)|is.nan(fids)] = "0"
          if(length(fids)){
            gs.t[match(names(fids),gs.t$flywire.xyz),"flywire.id"] = fids
            gs.t$flywire.xyz = apply(gs.t[,c("fw.x","fw.y",'fw.z')],1,paste,sep=";",collapse=";")
            # Update
            update = rbind(gs.t[,used.cols],gs0[,used.cols])
            rownames(update) = NULL
            continue =  TRUE
            while (continue) {
              gsheet_manipulation(FUN = googlesheets4::write_sheet,
                                  wait = 20,
                                  update[0,],
                                  ss = selected_sheet,
                                  sheet = tab)
              batches = split(1:nrow(update), ceiling(seq_along(1:nrow(update))/500))
              errors = c()
              for(i in batches){
                tryCatch(gsheet_manipulation(FUN = googlesheets4::sheet_append,
                                    wait = 20,
                                    data = update[min(i):max(i),],
                                    ss = selected_sheet,
                                    sheet = tab),
                         errors = function(e) errors <<- c(errors, i))
              }
              continue = ifelse(length(errors),TRUE,FALSE)
            }
            # Now continue processing
            gs.t = gs.t[,colnames(gs.t)%in%chosen.columns]
            for(col in chosen.columns){
              if(is.null(gs.t[[col]])){
                gs.t[[col]] = NA
              }
            }
            gs = plyr::rbind.fill(gs.t,gs)
          }
        }
      }
    }
    if(is.null(gs$flywire.xyz)){
      gs$flywire.xyz = apply(gs[,c("fw.x","fw.y",'fw.z')],1,paste,sep=";",collapse=";")
    }
  }
  # Make this unique, but keep row with most information
  master = gs
  master = master[!is.na(master$fw.x),]
  master = master[!is.na(master$fw.y),]
  master = master[!is.na(master$fw.z),]
  master$filled = apply(master, 1, function(r) sum(!is.na(r)))
  master = master[order(master$filled,decreasing = TRUE),]
  master = master[!duplicated(master$flywire.xyz),]
  rownames(master) = master$flywire.xyz
  master = subset(master, !is.na(master$flywire.id))
  master$filled = NULL
  master
}


# hidden
gsheet_manipulation <- function(FUN,
                                wait = 10,
                                ...,
                                return = FALSE){
  sleep = wait
  success = FALSE
  while(!success){
    g = tryCatch(FUN(...),
                 error = function(e){
                   cat(as.character(e))
                   return(NULL)
                 })
    if(!is.null(g)){
      success = TRUE
    }else{
      message("Google read/write failure(s), re-trying in ", sleep," seconds ...")
      Sys.sleep(sleep)
      sleep = sleep + wait
      if(sleep > 600){
        slep <- 600
      }
    }
  }
  if(return){
    if(is.data.frame(g)){
      g = unlist_df(g)
    }
    return(g)
  }
}


# Skeletonise neurons in parallel from a folder of obj files
skeletor_batch <- function(obj, swc, numCores = 1, max.file.size = 1000000000, ...){
  obj.files = list.files(obj, pattern = "obj$", full.names = TRUE)
  ids = obj.files[sapply(obj.files, file.size) < max.file.size]
  big = setdiff(obj.files,ids)
  if(length(big)){
    warning("Dropping ", length(big), " .obj files larger than ", max.file.size, " bytes")
  }
  batches = split(ids, round(seq(from = 1, to = numCores, length.out = length(ids))))
    foreach.skeletons <- foreach::foreach (batch = seq_along(batches)) %dopar% {
      neuron.ids = batches[[batch]]
      j = tryCatch({
        skels = fafbseg::skeletor(neuron.ids, ...)
        nat::write.neurons(skels, dir=swc, format='swc', Force = FALSE)
        skels},
                   error = function(e){
                     cat(as.character(e))
                     message(paste(fw.ids,collapse=","))
                     NULL
                   })
    }
    isnl = sapply(foreach.skeletons, nat::is.neuronlist)
    do.call(c, foreach.skeletons[isnl])
}

# hidden
download_neuron_obj_batch <- function(ids, numCores = 1, ratio = 1, save.obj = "obj"){
    batches = split(ids, round(seq(from = 1, to = numCores, length.out = length(ids))))
    foreach.skeletons <- foreach::foreach (batch = seq_along(batches)) %dopar% {
      neuron.ids = batches[[batch]]
      j = tryCatch(fafbseg::download_neuron_obj(segments = neuron.ids,
                                                ratio = ratio,
                                                save.obj = save.obj),
                   error = function(e){
                     cat(as.character(e))
                     NULL
                   })
    }
}

# hidden
update.neuronlistfh <- function(x, rds, ...){
  data = file.path(dirname(rds),"data")
  old.neurons = nat::read.neuronlistfh(rds)
  combined.neurons = union(old.neurons, x)
  given.neurons = nat::as.neuronlistfh(combined.neurons, dbdir= data, WriteObjects = "missing", ...)
  nat::write.neuronlistfh(given.neurons, file=rds, overwrite=TRUE, ...)
}


