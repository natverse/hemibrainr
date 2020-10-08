#############################################################################
################################ Utilities ##################################
#############################################################################

# google upload neuronlistfh
googledrive_upload_neuronlistfh <- function(x,
                                            team_drive = "hemibrain",
                                            file_name = "neurons.rds",
                                            folder = "flywire_neurons",
                                            subfolder = NULL,
                                            numCores = 1){
  # Get drive
  td = googledrive::team_drive_get(team_drive)
  drive_td = googledrive::drive_find(type = "folder", team_drive = td)
  gfolder= subset(drive_td,drive_td$name==folder)[1,]
  if(!is.null(subfolder)){
    sub = googledrive::drive_ls(path = gfolder, type = "folder", team_drive = td)
    gfolder = subset(sub, sub$name ==subfolder )[1,]
  }
  sub = googledrive::drive_ls(path = gfolder, team_drive = td)
  gfile = subset(sub, sub$name ==file_name )
  if(nrow(gfile)){
    save.position = googledrive::as_id(gfile[1,]$id)
  }else{
    save.position = gfolder
  }

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
  temp = tempdir(check=TRUE)
  temp.data = paste0(temp,"/data")
  dir.create(temp.data)
  on.exit(unlink(temp.data, recursive=TRUE))
  temp.rds = paste0(temp,"/",file_name)
  nl = nat::as.neuronlistfh(x, dbdir= temp.data, WriteObjects = "yes")
  nat::write.neuronlistfh(nl, file= temp.rds, overwrite=TRUE)

  # upload
  t.list.master = list.files(temp.data,full.names = TRUE)
  error.files = c()
  sub.data = googledrive::drive_ls(path = t.folder.data, team_drive = td)
  if(numCores>1){
    batch = 1
    batches = split(t.list.master, round(seq(from = 1, to = numCores, length.out = length(t.list.master))))
    foreach.skeletons <- foreach::foreach (batch = 1:numCores) %dopar% {
      t.list = batches[[batch]]
      for(t.neuron.fh.data.file in t.list){
        t = basename(t.neuron.fh.data.file)
        if(t%in%sub.data$name){
          save.data = googledrive::as_id(subset(sub.data, sub.data$name==t)[1,]$id)
        }else{
          save.data =  t.folder.data
        }
        e = tryCatch(google_drive_place(media = t.neuron.fh.data.file,
                                        path = save.data,
                                        verbose = FALSE),
                     error = function(e){
                       message(e)
                       NULL
                     } )
        if(is.null(e)){
          warning(e)
          error.files = c(error.files,t.neuron.fh.data.file)
        }
      }
    }
  }else{
    pb <- progress::progress_bar$new(
      format = "  downloading [:bar] :current/:total eta: :eta",
      total = length(t.list), clear = FALSE, show_after = 1)
    for(t.neuron.fh.data.file in t.list.master){
      pb$tick()
      t = basename(t.neuron.fh.data.file)
      if(t%in%sub.data$name){
        save.data = googledrive::as_id(subset(sub.data, sub.data$name==t)[1,]$id)
      }else{
        save.data =  t.folder.data
      }
      e = tryCatch(google_drive_place(media = t.neuron.fh.data.file,
                                      path = save.data,
                                      verbose = FALSE),
                   error = function(e){
                     message(e)
                     NULL
                   } )
      if(is.null(e)){
        warning(e)
        error.files = c(error.files,t.neuron.fh.data.file)
      }
    }
  }
  if(length(error.files)){
    warning("Failed to upload: ", length(error.files)," files")
  }
  # upload
  google_drive_place(media = temp.rds,
            path = save.position,
            verbose = TRUE)
}

# Google drive
google_drive_place <- function(media,
                               path,
                               verbose = TRUE,
                               ...){
  # If path is folder, check contents to duplicates
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

# upload nblast matrix
googledrive_upload_nblast<- function(x,
                                     team_drive = "hemibrain",
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
  temp <- tempdir(check=TRUE)
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


# collapse matrix by names
collapse_matrix_by_names <- function(M, FUN = mean, ...){
  M = apply(M, 2, function(x) tapply(x, rownames(M), FUN, ...))
  M = t(apply(t(M), 2, function(x) tapply(x, colnames(M), FUN, ...)))
  M
}


# Add missing columns to a df
add_blanks <- function(df, missing, add = ""){
  for(m in missing){
    df[[m]] <- add
  }
  df
}

# hidden
nullToNA <- function(x) {
  if(is.list(x)){
    x[sapply(x, is.null)] <- NA
  }else{
    x = sapply(x, function(y) ifelse(is.null(y)|!length(y), NA, y))
    if(!length(x)){
      x = NA
    }
  }
  x
}

# hidden
break_into_subtrees <- function(x, prune = FALSE){
  if(!nat::is.neuron(x)){
    stop("x must be a neuron object")
  }else{
    y = x
  }
  if(prune){
    nulls = subset(rownames(x$d), x$d$Label %in% c(4,7))
    y = nat::prune_vertices(x, verticestoprune = as.numeric(nulls), invert = FALSE)
  }
  nlist = nat::neuronlist()
  if(y$nTrees>1){
    for(i in 1:y$nTrees){
      segs = y$SubTrees[[i]]
      seg.points = unique(unlist(segs))
      d = y$d[seg.points,]
      if(nrow(d)>1){
        rownames(d) = 1:nrow(d)
        n = nat::as.neuron(d)
        n$orig.PointNo = d[,"PointNo"]
        n$orig.indices = match(n$orig.PointNo, x$d$PointNo)
        nlist  = c(nlist, nat::as.neuronlist(n))
      }
    }
    nlist
  }else{
    nat::as.neuronlist(y)
  }
}

# hidden
is.hxsurf <- function(x){
  "hxsurf"%in%class(x)
}

#' @importFrom nat progress_natprogress
#' @export
nat::progress_natprogress

# Get the Nth highest number
maxN <- function(x, N=2){
  len <- length(x)
  if(N>len){
    warning('N greater than length(x).  Setting N=length(x)')
    N <- length(x)
  }
  sort(x,partial=len-N+1)[len-N+1]
}


# hidden
maxout <- function(x, max){
  x[x>max] = max
  x[is.infinite(x)|is.null(x)|is.na(x)] = 0
  x
}

# hidden
lengthnorm <- function(x){
  sum(x, na.rm = TRUE)/length(x)
}

# hidden
hemibrain_neuron_class <- function (x){
  if(nat::is.neuronlist(x)){
    x = nat::nlapply(x,hemibrain_neuron_class)
  } else {
    our_classes=c("neuprintneuron","catmaidneuron","neuron","list")
    # we assume that any extra classes besides those 4 should be at the front
    # of the list
    extra_classes=setdiff(class(x), our_classes)
    class(x) = c(extra_classes, our_classes)
  }
  x
}

# hidden
which.consecutive <- function(Vec,
                              only.jumps = FALSE,
                              run = c("all","min","max","minmax")){
  run = match.arg(run)
  if(is.logical(Vec)){
    Vec = which(Vec)
  }else if(!is.integer(Vec)){
    Vec = as.integer(Vec)
  }
  if(only.jumps){
    if(length(which.consecutive(Vec)) == length(Vec)){
      return(NULL)
    }
  }
  Breaks <- c(0, which(diff(Vec) != 1), length(Vec))
  if(run!="all"){
    cons <- lapply(seq(length(Breaks) - 1),
                  function(i) c(Vec[(Breaks[i] + 1):Breaks[i+1]]))
    if(run=="minmax"){
      unlist(lapply(cons, function(c) c(min(c),max(c))))
    }else if (run=="min"){
      unlist(lapply(cons, function(c) min(c)))
    }else if (run=="max"){
      unlist(lapply(cons, function(c) max(c)))
    }
  }else{
    unlist(lapply(seq(length(Breaks) - 1),
                  function(i) Vec[(Breaks[i] + 1):Breaks[i+1]]))
  }
}

# hidden
change_points <- function(x, v, only.jumps = FALSE, run = "min"){
  eps = nat::endpoints(x)
  segs = x$SegList
  df=data.frame(node=unlist(segs), seg=rep(seq_along(segs), sapply(segs, length)))
  bb=by(df$node%in%v, df$seg, function(x) any(x))
  segs.d=segs[bb]
  s.d = unique(unlist(lapply(segs.d, function(seg) seg[which.consecutive(seg %in% v, only.jumps = only.jumps, run = run)])))
  s.d = setdiff(s.d, eps)
}

# hidden
is.issue <- function(x){
  x = x[1]
  if(length(x)){
    if(!is.na(x)){
      if(!is.nan(x)){
        if(x!=""){
          FALSE
        }else{TRUE}
      }else{TRUE}
    }else{TRUE}
  }else{TRUE}
}

# hidden
carryover_tags <- function(x, y){
  y$tags = x$tags
  y
}

# hidden
carryover_labels <- function(x, y){
  y$d$Label = x$d$Label[match(y$d$PointNo,x$d$PointNo)]
  y$connectors$Label = x$connectors$Label[match(y$connectors$PointNo,x$d$PointNo)]
  y
}

# hidden
purify <- function(x){
  as.character(unname(unlist(nullToNA(c(x)))))
}

# hidden
remove_duplicates <- function(manual){
  delete = c()
  for(bi in unique(manual$bodyid)){
    m = manual[manual$bodyid==bi,]
    dupe = which(duplicated(m$bodyid))
    droot = which(m[dupe,]$point=="root")
    del = rownames(m)[dupe>=droot]
    delete = c(delete,del)
  }
  manual[setdiff(rownames(manual),delete),]
}

# hidden
suppress <- function(x, ...){
  suppressWarnings(suppressMessages(x, ...), ...)
}

# hidden
replace_with_none <- function(x, FUN = is.na){
  x[FUN(x)] = "none"
  x
}

# hidden
unlist_df <- function(df){
  data = as.data.frame(df, stringsAsFactors = FALSE)
  if(nrow(df)&ncol(df)){
    data = as.data.frame(df, stringsAsFactors = FALSE)
    data = apply(data,2,function(c) unlist(nullToNA(c)))
    data = as.data.frame(unlist(data), stringsAsFactors = FALSE)
    dimnames(data) = dimnames(df)
    data
  }
  data[] <- lapply(data, as.character)
  data
}

# hidden
very_simple_connectivity <- function(conn.df){
  conn = apply(conn.df, 1, function(row) c(
    name = row["name"],
    bodyid = ifelse(is.na(row["input"]),row["output"],row["input"]),
    type = row["type"],
    total.weight = sum(as.numeric(row[grepl("weight",names(row))]), na.rm = TRUE) ))
  conn = as.data.frame(t(conn), stringsAsFactors = FALSE)
  colnames(conn) = gsub( "\\..*", "",colnames(conn) )
  conn
}

# We can make another version of the matrix that contains essentially all
# the useful information, but compresses much smaller (~16x at time of writing)
save_compressed_nblast_mat <- function(x,
                                       file = NULL,
                                       threshold=-0.5,
                                       digits=3,
                                       format=c("rda", "rds"),
                                       ...) {
  objname=deparse(substitute(x))
  format=match.arg(format)
  colnames(x) = gsub("_m$","",colnames(x)) # factor in mirrored hemibrain neurons
  rownames(x) = gsub("_m$","",rownames(x))
  x = apply(x, 2, function(i) tapply(x, rownames(x), sum, na.rm = TRUE))
  x = t(apply(t(x), 2, function(i) tapply(i, colnames(x), sum, na.rm = TRUE)))
  x[x<threshold]=threshold
  x=round(x, digits=digits)
  newobjname <- paste0(objname, ".compressed")
  fname <- paste0(paste0(file, newobjname), ".", format)
  message("Resaving a compressed version of ", objname, " to ", fname)
  if(format=="rds") {
    saveRDS(x, file=fname, ...)
  } else {
    # the object being saved must have the name that you would like it
    # to have when it is loaded again
    assign(newobjname, x)
    save(list = newobjname, file = fname, ...)
  }
}

# Transform into template brainspace
xform.mesh3d <- function(mesh3d, reg = reg){
  points = t(mesh3d$vb)[,1:3]
  nat::xyzmatrix(mesh3d)  = nat::xformpoints(reg = reg,
                                             points = points, transformtype = c("warp"), direction = NULL, FallBackToAffine = FALSE)
  mesh3d
}
xform_brain.mesh3d <- function(mesh3d,
                               sample = nat.templatebrains::regtemplate(mesh3d),
                               reference){
  points = t(mesh3d$vb)[,1:3]
  nat::xyzmatrix(mesh3d)  = nat.templatebrains::xform_brain(x = points, sample=sample, reference = reference)
  mesh3d
}
xform.shapelist <- function(shapelist, reg = reg){
  shapelist.transformed = lapply(shapelist,xform.mesh3d, reg = reg)
  class(shapelist.transformed) = c("shape3d","shapelist3d")
  shapelist.transformed
}
xform_brain.shapelist <- function(shapelist,
                                  sample = nat.templatebrains::regtemplate(shapelist),
                                  reference){
  shapelist.transformed = lapply(shapelist,xform_brain.mesh3d, sample = sample, reference = reference)
  class(shapelist.transformed) = c("shape3d","shapelist3d")
  shapelist.transformed
}
scale_points <-function(x, scaling = (8/1000)){
  nat::xyzmatrix(x) = nat::xyzmatrix(x)*scaling
  x
}
scale_points.shapelist <- function(shapelist,
                                   scaling = (8/1000)){
  shapelist.transformed = lapply(shapelist,scale_points, scaling = scaling)
  class(shapelist.transformed) = c("shape3d","shapelist3d")
  shapelist.transformed
}

# Cut neurons by bounding box
subbbx <- function(n.dps, min.x, max.x, min.y, max.y, min.z, max.z, ret='inside'){
  points = nat::xyzmatrix(n.dps)
  inside = (
    min.x < points[, "X"] & max.x > points[, "X"]
    & min.y < points[, "Y"] & max.y > points[, "Y"]
    & min.z < points[, "Z"] & max.z > points[, "Z"]
  )
  # or subset(n.dps, inside, invert=ret != 'inside')
  if (ret == 'inside') {
    subset(n.dps, inside)
  }
  else {
    subset(n.dps,!inside)
  }
}

# Subset FAFB neurons to FIB bounding box
fafb_hemibrain_cut <- function(x, ...){
  x = nat::nlapply(x, subbbx,
                                min.x=312048/1e3, max.x=601733/1e3,
                                min.y=71265/1e3, max.y=319018/1e3,
                                min.z=4315/1e3, max.z=270859/1e3,
                                ret='inside', ...)
  # Remove a chunk of left antennal lobe that's missing in hemibrain
  x = nat::nlapply(x, subbbx,
                                min.x=519300/1e3, max.x=601733/1e3,
                                min.y=259570/1e3, max.y=319018/1e3,
                                min.z=0/1e3, max.z=300000/1e3,
                                ret='outside', ...)
  # drop any FAFB neurons that don't have at least 5 vertices left
  x = x[nat::nvertices(x)>=5]
  x
}

# use foreach to process in parallel
#' @importFrom nat.templatebrains regtemplate xform_brain
xform_brain_parallel <- function(x,
                                 numCores = 2,
                                 sample = regtemplate(x),
                                 reference,
                                 ...){
  batch = 1
  batches = split(x, round(seq(from = 1, to = numCores, length.out = length(x))))
  foreach.nl <- foreach::foreach (batch = 1:numCores) %dopar% {
    y = batches[[batch]]
    j = java_xform_brain(y,
                    reference =reference, sample = sample,
                    .parallel = FALSE,...)
  }
  neurons = do.call(c, foreach.nl)
  neurons
}

# Strip meshes
strip_meshes<-function(x){
  strip_mesh.neuron <- function(x){
    x$mesgh3d = NULL
    class(x) = setdiff(class(x),"neuronmesh")
    x
  }
  nat::nlapply(x, strip_mesh.neuron)
}

# Remove unused filehash files
googledrive_clean_neuronlistfh <- function(team_drive = "hemibrain"){
  # Save .rds files locally
  temp = tempdir(check=TRUE)
  temp.data = paste0(temp,"/RDS/")
  dir.create(temp.data)
  on.exit(unlink(temp.data, recursive=TRUE))

  # Iterate to find .rds files
  td = googledrive::team_drive_get(team_drive)
  drive_td = googledrive::drive_find(type = "folder", team_drive = td)
  remove = c()
  for(folder in drive_td$id){
    sub = googledrive::drive_ls(path = googledrive::as_id(folder), type = "folder", team_drive = td)
    if("data" %in% sub$name){
      sub = googledrive::drive_ls(path = googledrive::as_id(folder), team_drive = td)
      if(sum(grepl("rds$",sub$name))>1){
        message("Cleaning ", subset(drive_td,drive_td$id==folder)$name)
        fsub = subset(sub, grepl("rds$",sub$name))
        fdata = subset(sub, grepl("data$",sub$name))[1,]
        all.keys = c()
        for(file in fsub$id){
          fnam = subset(fsub, fsub$id==file)$name
          path = paste0(temp.data,fnam)
          googledrive::drive_download(file = googledrive::as_id(file), path = path, overwrite = TRUE, verbose = FALSE)
          a = readRDS(path)
          b = attributes(a)
          keys = b$keyfilemap
          all.keys = c(all.keys,keys)
        }
        data = googledrive::drive_ls(path = googledrive::as_id(fdata$id), team_drive = td)
        delete = setdiff(data$name,all.keys)
        delete = googledrive::as_id(subset(data, data$name%in%delete)$id)
        remove = unique(c(remove, delete))
      }
    }
  }
  # rm
  if(length(remove)){
    message("Removing files ...")
    googledrive::drive_rm(remove, verbose = TRUE)
  }
}

# Obtains flywire IDs from flywire posirions on google sheets
# The sheets must either have columns fw.x, fw.y, fw.z
# or one column, flywire.xyz with number separated by commas, in the form: x,y,z
flywire_flagged <- function(selected_sheets = c("1OSlDtnR3B1LiB5cwI5x5Ql6LkZd8JOS5bBr-HTi0pOw",
                                                "1rzG1MuZYacM-vbW7100aK8HeA-BY6dWAVXQ7TB6E2cQ",
                                                "1spGSuhUX6Hhn-8HH0U_ArIWUuPpMBFNjIjeSSh_MFVY"),
chosen.columns = c("fw.x","fw.y",'fw.z', 'flywire.xyz',
                   "flywire.id", "skid",
                   "FAFB.xyz", "side",
                   "ItoLee_Hemilineage", "Hartenstein_Hemilineage",
                   "hemibrain_match"),
numCores = 1){
  # Read selected sheets and extract positions for flywire neurons
  # One xyz position is enough to identify a neuron
  # We do this because flywire.ids change all of the time
  gs = data.frame()
  for(selected_sheet in selected_sheets){
    ## Read google sheets and extract glywire neuron positions
    tabs = googlesheets4::sheet_names(selected_sheet)
    for(tab in tabs){
      gs.t = gsheet_manipulation(FUN = googlesheets4::read_sheet,
                                              ss = selected_sheet,
                                              sheet = tab,
                                              guess_max = 3000,
                                              return = TRUE)
      if(ncol(gs.t)&&sum(grepl("fw.x|flywire.xyz",colnames(gs.t)))>0){
        # Separate x,y,z positions
        gs1 = subset(gs.t, (is.na(gs.t$flywire.xyz)||is.null(gs.t$flywire.xyz)) && !is.na(gs.t$fw.x))
        gs2= subset(gs.t, !is.na(gs.t$flywire.xyz))
        gs1$flywire.xyz = apply(gs1[,c("fw.x","fw.y",'fw.z')],1,paste,sep=",",collapse=",")
        if(nrow(gs2)){
          positions.gs = sapply(gs2$flywire.xyz,strsplit,",|/|;")
          ruleofthree = sapply(positions.gs,function(p) length(p)==3)
          positions.gs = positions.gs[ruleofthree]
          positions.gs = do.call(rbind, positions.gs)
          positions.gs = apply(positions.gs,2,as.numeric)
          gs2[ruleofthree,c("fw.x","fw.y",'fw.z')] = positions.gs
        }
        gs.t = rbind(gs2,gs1)
        # Get flywire IDs from these positions
        batch = 1
        batches = split(1:nrow(gs.t), round(seq(from = 1, to = numCores, length.out = nrow(gs.t))))
        foreach.ids <- foreach::foreach (batch = 1:numCores) %dopar% {
          pos = gs.t[batches[[batch]],]
          j <- tryCatch({i = fafbseg::flywire_xyz2id(pos[,c("fw.x","fw.y",'fw.z')], rawcoords = TRUE)
          names(i) = pos$flywire.xyz
          if(sum(i==0)>0){
            i[i==0] =  fafbseg::flywire_xyz2id(pos[i==0,], rawcoords = FALSE)
          }
          i}, error = function(e) NULL)
        }
        fids = unlist(foreach.ids)
        fids[fids=="0"] = NA
        gs.t[match(names(fids),gs.t$flywire.xyz),"flywire.id"] = fids
        # Update
        rownames(gs.t) = NULL
        googlesheets4::write_sheet(gs.t[0,],
                                   ss = selected_file,
                                   sheet = tab)
        batches = split(1:nrow(gs.t), ceiling(seq_along(1:nrow(gs.t))/500))
        for(i in batches){
          gsheet_manipulation(FUN = googlesheets4::sheet_append,
                              data = gs.t[min(i):max(i),],
                              ss = selected_sheet,
                              sheet = tab)
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
    if(is.null(gs$flywire.xyz)){
      gs$flywire.xyz = apply(gs[,c("fw.x","fw.y",'fw.z')],1,paste,sep=",",collapse=",")
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


# Update neurons xyz locations
java_xform_brain <- function(x,
                             sample = regtemplate(x),
                             reference,
                             method = "rJava",
                             progress.rjava=TRUE,
                             ...){
  # Transform treenodes
  points = nat::xyzmatrix(x)
  t = xform_brain(points, reference = reference, sample = sample, method = method, progress.rjava=progress.rjava, ...)
  nat::xyzmatrix(x) = t
  # Transform synapses
  syns = lapply(names(x), function(n){
    conn = x[[n]]$connectors
    conn$id = n
    conn
  })
  conns = do.call(rbind, syns)
  if(nrow(conns)){
    conns.t = xform_brain(nat::xyzmatrix(conns), reference = reference, sample = sample, method = method, progress.rjava=progress.rjava, ...)
    nat::xyzmatrix(conns) = conns.t
    x = add_field_seq(x,names(x),field="id")
    x = nat::nlapply(x, function(n){
      n$connectors = subset(conns, conns$id == n$id)
      n$id = NULL
      n$connectors$id = NULL
      n
    })
  }
  x
}



