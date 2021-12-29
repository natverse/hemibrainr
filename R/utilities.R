#############################################################################
################################ Utilities ##################################
#############################################################################

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
  df=data.frame(node=unlist(segs), seg=rep(seq_along(segs), sapply(segs, length)), stringsAsFactors = FALSE)
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
is.nrowlength <- function(x){
  if(length(x)){
    if(nrow(x)){
      ans <- TRUE
    }else{
      ans <- FALSE
    }
  }else{
    ans <- FALSE
  }
  ans
}

# hidden
carryover_tags <- function(x, y){
  y$tags = x$tags
  y
}

# hidden
carryover_labels <- function(x, y = NULL){
  if(is.null(y)){
    if(nrow(x$connectors)){
      tid = ifelse("PointNo"%in%colnames(x$connectors),"PointNo","treenode_id")
      x$connectors$Label = x$d$Label[match(x$connectors[[tid]],x$d$PointNo)]
    }
    x
  }else{
    tidx = ifelse("PointNo"%in%colnames(x$connectors),"PointNo","treenode_id")
    tidy = ifelse("PointNo"%in%colnames(x$connectors),"PointNo","treenode_id")
    y$d$Label = x$d$Label[match(y$d$PointNo,x$d$PointNo)]
    y$connectors$Label = x$d$Label[match(y$connectors[[tidy]],x$d$PointNo)]
    y
  }
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
    data = apply(data,2,function(c) unlist(nullToNA(c)))
    if(nrow(df)==1){
      data = t(data)
    }
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


# use foreach to process in parallel
#' @importFrom nat.templatebrains regtemplate xform_brain
xform_brain_parallel <- function(x,
                                 numCores = 2,
                                 sample = regtemplate(x),
                                 reference,
                                 ...){
  batch = 1
  batches = split(x, round(seq(from = 1, to = numCores, length.out = length(x))))
  foreach.nl <- foreach::foreach (batch = 1:length(batches)) %dopar% {
    y = batches[[batch]]
    j = java_xform_brain(y,
                         reference =reference,
                         sample = sample,
                         .parallel = FALSE,
                         ...)
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

# Update neurons xyz locations
java_xform_brain <- function(x,
                             sample = regtemplate(x),
                             reference,
                             method = "rJava",
                             progress.rjava=TRUE,
                             ...){
  # Transform treenodes
  points = nat::xyzmatrix(x)
  t = nat.templatebrains::xform_brain(points,
                                      reference = reference,
                                      sample = sample,
                                      method = method,
                                      progress.rjava=progress.rjava,
                                      ...)
  nat::xyzmatrix(x) = t
  # Transform synapses
  syns = lapply(names(x), function(n){
    conn = x[[n]]$connectors
    if(!is.null(conn)|length(conn)){
      conn$id = n
    }
    conn
  })
  conns = do.call(plyr::rbind.fill, syns)
  xyz.good = tryCatch(nat::xyzmatrix(conns), error = function(e) NULL)
  if(!is.null(xyz.good)|length(xyz.good)){
    conns.t = nat.templatebrains::xform_brain(xyz.good, reference = reference, sample = sample, method = method, progress.rjava=progress.rjava, ...)
    nat::xyzmatrix(conns) = conns.t
    x = add_field_seq(x,names(x),field="id")
    x = nat::nlapply(x, function(n){
      if(!is.null( n$id)){
        n$connectors = subset(conns, conns$id == n$id)
        if(!is.null(n$connectors)){
          n$id = NULL
          n$connectors$id = NULL
        }
      }
      n
    })
  }
  x
}

# hidden
## Save with given name
saveit <- function(..., file) {
  x <- list(...)
  save(list=names(x), file=file, envir=list2env(x))
}

# remove filehash files
remove_unused_filehash <- function(path,
                                   dbClass = c("RDS", "RDS2", "DB1", "ZIP"),
                                   meta = NULL,
                                   id = NULL){
  dbClass = match.arg(dbClass)
  if(dbClass=="ZIP"){
    message("remove_unused_filehash not implemented for .zip files supporting neuronlistz objects")
    return(invisible())
  }
  for(p in path){
    if(dbClass=="DB1"){
      files = list.files(path, pattern = ".rds$", full.names = TRUE)
      for(rds in files){
        nlfh = nat::read.neuronlistfh(rds)
        if(!is.null(id)&!is.null(meta)){
          if(!id%in%colnames(meta)){
            stop("id must be in colnames(meta)")
          }
          inmeta = names(nlfh)%in%meta[[id]]
          if(sum(!inmeta)){
            nlfh = nlfh[inmeta]
            update.neuronlistfh(nat::as.neuronlist(nlfh),
                                file=rds,
                                dbClass = "DB1")
          }
        }
      }
    }else{
      data = file.path(p,"data")
      if(!dir.exists(data)){
        warning("data folder for neuronlistfh object does not exit at: ", data)
      }else{
        files = list.files(path, pattern = ".rds$", full.names = TRUE)
        if(!length(files)){
          warning("No .rds files found at: ", path)
        }else{
          all.keys = c()
          for(f in files){
            a = readRDS(f)
            b = attributes(a)
            keys = b$keyfilemap
            all.keys = c(all.keys,keys)
          }
          all.fh = list.files(data)
          all.fh = all.fh[!grepl("zip$|csv$|swc$|rds$",all.fh)]
          delete = setdiff(all.fh,all.keys)
          message("Deleting ", length(delete), " files")
          delete = file.path(data,delete)
          file.remove(delete)
        }
      }
    }
  }
}

# Skeletonise neurons in parallel from a folder of obj files
skeletor_batch <- function(obj, swc, numCores = 1, max.file.size = 1000000000, ...){
  if(dir.exists(obj[1])){
    obj.files = list.files(obj, pattern = "obj$", full.names = TRUE)
  }else{
    obj.files = obj
  }
  ids = obj.files[sapply(obj.files, file.size) < max.file.size]
  if(!length(ids)){
    return(NULL)
  }
  big = setdiff(obj.files,ids)
  if(length(big)){
    warning("Dropping ", length(big), " .obj files larger than ", max.file.size, " bytes")
  }
  batches = split(ids, round(seq(from = 1, to = numCores, length.out = length(ids))))
  batch = 0
  foreach.skeletons <- foreach::foreach (batch = seq_along(batches)) %dopar% {
    neuron.ids = batches[[batch]]
    j = tryCatch({
      skels = suppressMessages(fafbseg::skeletor(neuron.ids, save.obj = NULL, mesh3d = FALSE, ...))
      skels[,"id"] = names(skels) = basename(gsub("\\.obj","",names(skels)))
      nat::write.neurons(skels, dir=swc, format='swc', Force = FALSE)
      skels},
      error = function(e){
        warning(paste(as.character(e), neuron.ids,collapse=",") )
        NULL
      })
  }
  isnl = sapply(foreach.skeletons, nat::is.neuronlist)
  do.call(c, foreach.skeletons[isnl])
}

# hidden
download_neuron_obj_batch <- function(ids, numCores = 1, ratio = 1, save.obj = "obj"){
  if(!length(ids)){
    return(NULL)
  }
  batches = split(ids, round(seq(from = 1, to = numCores, length.out = length(ids))))
  batch = 0
  foreach.skeletons <- foreach::foreach (batch = seq_along(batches)) %dopar% {
    neuron.ids = batches[[batch]]
    j = tryCatch(suppressMessages(fafbseg::download_neuron_obj(segments = neuron.ids,
                                              ratio = ratio,
                                              save.obj = save.obj)),
                 error = function(e){
                   warning(as.character(e))
                   NULL
                 })
    NULL
  }
}

# hidden
update.neuronlistfh <- function(x = NULL,
                                file,
                                dbClass = c("RDS", "RDS2", "DB1","HDF5", "ZIP"),
                                remote = NULL,
                                localdir = NULL,
                                pref.meta = NULL,
                                meta = NULL,
                                compress = TRUE,
                                id = 'flywire_xyz',
                                swc = NULL,
                                ...){
  dbClass = match.arg(dbClass)
  if(grepl("\\.zip$",file)){
    dbClass = "ZIP"
  }
  if(dbClass=="DB1"){
    data = gsub("\\.rds","_datafile", file)
    if(file.exists(paste0(data,"___LOCK"))){
      file.remove(paste0(data,"___LOCK"))
    }
    WriteObjects = 'yes'
  }else if (dbClass == "HDF5"){
    data =  file
  }else{
    data = paste0(dirname(file),"/","data/")
    WriteObjects = "missing"
  }
  if(file.exists(file)){
    if(dbClass == "HDF5"){
      old.neurons = tryCatch(nat.hdf5::read.neurons.hdf5(file), error = function(e) {
        warning(e)
        message("original neuron file cannot link to data, overwriting ", file)
        NULL
      })
    }else if(dbClass == "ZIP"){
      old.neurons = tryCatch(nat::read.neurons(file), error = function(e) {
        warning(e)
        message("original neuron file cannot link to data, overwriting ", file)
        NULL
      })
    }else{
      old.neurons = tryCatch(nat::read.neuronlistfh(file, localdir = localdir), error = function(e){
        warning(e)
        message("original neuron file cannot link to data, overwriting ", file)
        NULL
      })
    }
    if(!is.null(old.neurons)&&length(old.neurons)&&nat::is.neuronlist(old.neurons)){
      if(!is.null(attr(old.neurons,"df"))){
        just.old = setdiff(names(old.neurons),names(x))
        if(length(just.old)){
          old.neurons = tryCatch(old.neurons[just.old], error = function(e){
            message(as.character(e))
            NULL
          })
          if(is.null(x)){
            x = old.neurons
          }else{
            x = nat::union(x, old.neurons)
            message(length(x), " given neurons combined with ", length(old.neurons), " old neurons from extant: ", file)
          }
        }
      }
    }
  }
  if(nat::is.neuronlist(x)){
    if(dbClass=="DB1"){try(file.remove(data))}
    if(!is.null(meta)){
      if(!all(names(x)%in%rownames(meta))){
        warning("each neuron name in x should be a rowname in meta")
      }
      if(!is.null(rownames(meta))){
        in.meta = names(x) %in% rownames(meta)
        in.meta[is.na(in.meta)] = FALSE
        x = update_metdata(x, meta = meta, id = id)
      }
    }
    if(!is.null(pref.meta)){
      shared.cols = intersect(pref.meta,colnames(x[,]))
      if(length(shared.cols)){
        x[,] = x[,shared.cols]
      }
    }
    if(compress){
      x = squeeze_neuronlist(x, OmitFailures = TRUE)
    }
    if(dbClass == "HDF5"){
      given.neurons = nat.hdf5::write.neurons.hdf5(x, file = data, ...)
    }else if(dbClass == "ZIP"){
      temp.zip =  file.path(dirname(file),"temp_neuronlist_archive.zip")
      if(file.exists(file)){
        file.copy(from = file, to = temp.zip, overwrite = TRUE)
      }
      given.neurons = try(nat::write.neurons(x, dir = file, format='qs', include.data.frame = TRUE, Force = TRUE, ...), silent = FALSE)
      if(inherits("try-class",given.neurons)){
        try(file.copy(to = file, from = temp.zip, overwrite = TRUE), silent = TRUE)
        try(suppress(file.remove(temp.zip)), silent = TRUE)
        warning("could not create neuronlistz object")
        return(given.neurons)
      }
      try(suppress(file.remove(temp.zip)), silent = TRUE)
    }else{
      tryCatch({
        given.neurons = nat::as.neuronlistfh(x, dbdir = data, dbClass = dbClass, WriteObjects = WriteObjects, remote = remote, ...)
        }, error = function(e){
          file.remove(file)
          error("Neuronlitfh could not be generated, deleting original file")
        })
      if(dbClass=="DB1"){
        saveRDS(given.neurons, file = file)
      }else{
        nat::write.neuronlistfh(given.neurons, file=file, overwrite=TRUE, ...)
      }
    }
  }else{
    warning("could not create neuronlistfh/z object")
  }
  if(!is.null(swc) & !is.null(x)){
    ##  Write to SWC
    if(!dir.exists(swc)){
      dir.create(swc)
    }
    if(!is.null(meta)){
      if(all(names(x)%in%rownames(meta))){
        if("root_id"%in%colnames(meta)){
          fw.ids = unique(as.character(meta$root_id))
          swc.files = list.files(swc, full.names = TRUE)
          swc.delete = swc.files[!basename(swc.files)%in%paste0(fw.ids,".swc")]
          message("Removing ", length(swc.delete), " outdated .swc files")
          file.remove(swc.delete)
        }
      }
    }
    nat::write.neurons(x, dir=swc, format='swc', Force = FALSE, metadata=TRUE)
  }
}

# hidden
correct_id <-function(v){
  gsub(" ","",v)
}

# hidden
is64ToChar <- function(res){
  # convert 64 bit ints to char (safer but bigger)
  is64=sapply(res, bit64::is.integer64)
  if(any(is64)) {
    for(i in which(is64)) {
      res[[i]]=as.character(res[[i]])
    }
  }
  res
}

# hidden
try_with_time_limit <- function(expr, cpu = Inf, elapsed = Inf, error = NULL){
  on.exit({
    setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)
  })
  y <- try({setTimeLimit(cpu, elapsed, transient = TRUE); expr}, silent = TRUE)
  clear <- gc()
  if(inherits(y, "try-error")){
    error
  }else{
    y
  }
}

# hidden
fafb14_to_flywire_ids_timed <- function(x, only.biggest = FALSE){
  try_with_time_limit(fafbseg::fafb14_to_flywire_ids(search=x, only.biggest=only.biggest), elapsed = 1800, error = NA)
}

# hidden
squeeze_neuronlist <- function(x, digits=6, ...) {
  nat::nlapply(x, squeeze_neuron, digits=digits, ...)
}

# hidden
squeeze_neuron <- function(x, digits=6, ...) {
  stopifnot(nat::is.neuron(x)|nat::is.dotprops(x))
  check_package_available('bitsqueezr')
  x$d=squeeze_dataframe(x$d, exclude=c("X", "Y", "Z"), digits=digits, ...)
  if(!is.null(x$connectors)) {
    x$connectors=squeeze_dataframe(x$connectors, digits=digits, ...)
  }
  x
}

# hidden
squeeze_dataframe <- function(x, exclude=NULL, digits = 6, ...) {
  numcols <- names(x)[sapply(x, function(c) is.numeric(c) && !inherits(c, 'integer64'))]
  numcols=setdiff(numcols, exclude)
  for(i in numcols) {
    col=x[[i]]
    # does it look like an int, if so, make it one
    intcol=try(checkmate::asInteger(col), silent = TRUE)
    if((sum(is.na(col))==length(col))){
      x[[i]]=col
    }else if(is.integer(intcol)){
      x[[i]]=intcol
    }
    else x[[i]]=bitsqueezr::squeeze_bits(col, digits = digits, ...)
  }
  x
}

# hidden
check_package_available <- function(pkg) {
  if(!requireNamespace(pkg, quietly = TRUE)) {
    stop("Please install suggested package: ", pkg)
  }
}

# hidden
update_metdata <- function(neurons, meta, id){
  check_package_available('dplyr')
  check_package_available('snakecase')
  df = neurons[,]
  df = matchColClasses(meta, df)
  dfn = suppress(dplyr::left_join(df, meta))
  matched = match(dfn[[id]], meta[[id]])
  matched.good = !is.na(matched)
  shared.cols = intersect(colnames(dfn[,]),colnames(meta))
  dfn[matched.good,shared.cols] = meta[matched[matched.good],shared.cols]
  colnames(dfn) = snakecase::to_snake_case(dfn)
  attr(neurons,'df') = dfn
  neurons
}

# hidden
matchColClasses <- function(df1, df2) {
  sharedColNames <- colnames(df1)[colnames(df1) %in% colnames(df2)]
  if(length(sharedColNames)==1){
    sharedColTypes <- class(df1[,sharedColNames])
  }else{
    sharedColTypes <- sapply(df1[,sharedColNames], class)
  }
  for (n in 1:length(sharedColNames)) {
    class(df2[, sharedColNames[n]]) <- sharedColTypes[n]
  }
  return(df2)
}

# Change to snakecase
find_and_replace <- function(x = c("flywire.xyz",
                                          "flywire.id",
                                          "flywire.svid",
                                          "cell.type",
                                          "FAFB.xyz",
                                          "ItoLee_Hemilineage",
                                          "Hartenstein_Hemilineage",
                                          "hemibrain.match",
                                          "hemibrain.match.quality",
                                          "FAFB.hemisphere.match",
                                          "FAFB.hemisphere.match.quality",
                                          "ItoLee_Lineage",
                                          "Hartenstein_Lineage",
                                          "dataset",
                                          "total.outputs",
                                          "total.inputs",
                                          "axon.outputs",
                                          "dend.outputs",
                                          "axon.inputs",
                                          "dend.inputs",
                                          "total.outputs.density",
                                          "total.inputs.density",
                                          "axon.outputs.density",
                                          "dend.outputs.density",
                                          "axon.inputs.density",
                                          "dend.inputs.density",
                                          "total.length",
                                          "axon.length",
                                          "dend.length",
                                          "pd.length",
                                          "cable.length",
                                          "top.nt",
                                          "top.p",
                                          "total.pre",
                                          "total.post",
                                          "axon.pre",
                                          "axon.post",
                                          "dend.post",
                                          "dend.pre",
                                          "hemibrain.match.qualiy",
                                          "putative.classic.transmitter",
                                          "putative.other.transmitter",
                                          "FAFB.match",
                                          "FAFB.match.quality",
                                          "ct.layer"
                                          #, "soma.edit",
                                          # "edited.cable",
                                          #"orig.soma",
                                          #"orig.cut",
                                          #"soma.checked"
                                          ),
                             dir = getwd()){
  if(is.data.frame(x)){
    cols = colnames(x)
  }else{
    x = cols
  }
  cols = cols[nchar(cols)>1]
  for(col in cols){
    col_new = snakecase::to_snake_case(col)
    if (col == "flywire.id"){
      col_new = "root_id"
    }
    if (col_new == "itolee_hemilineage"){
      col_new = "ito_lee_hemilineage"
    }
    if (col_new == "itolee_lineage"){
      col_new = "ito_lee_lineage"
    }
    if(col!=col_new){
      message(sprintf("Replacing %s with %s within %s", col, col_new, dir))
      xfun::gsub_dir(dir = dir, pattern = col, replacement = col_new, rw_error = FALSE, mimetype =  '^text/')
    }
  }
}

# gsheet cols to snake cases
to_snake_case_gsheets  <- function(gsheets){
  for(ss in gsheets){
    message(ss)
    tabs = gsheet_manipulation(FUN = googlesheets4::sheet_names,
                               ss = ss,
                               return = TRUE,
                               Verbose = FALSE)
    for(ws in tabs){
      message(ss, " -> ", ws)
      df = gsheet_manipulation(
        FUN = googlesheets4::range_read,
        ss = ss,
        sheet = ws,
        range = NULL,
        col_names = TRUE,
        col_types = NULL,
        na = "",
        trim_ws = TRUE,
        skip = 0,
        n_max = 1,
        .name_repair = "unique"
      )
      cols_new = snakecase::to_snake_case(colnames(df))
      for (col in 1:length(cols_new)){
        if (cols_new[col] == "flywire.id"){
          cols_new[col] = "root_id"
        }
        if (cols_new[col] == "root_id"){
          cols_new[col] = "root_id"
        }
        if (cols_new[col] == "itolee_hemilineage"){
          cols_new[col] = "ito_lee_hemilineage"
        }
        if (cols_new[col] == "itolee_lineage"){
          cols_new[col] = "ito_lee_lineage"
        }
      }
      colnames(df) = cols_new
      gsheet_manipulation(
        FUN = googlesheets4::range_write,
        ss = ss,
        data = df,
        sheet = ws,
        range = NULL,
        col_names = TRUE,
        reformat = FALSE
      )
    }
  }
}


