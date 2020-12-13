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

# We can make another version of the matrix that contains essentially all
# the useful information, but compresses much smaller (~16x at time of writing)
save_compressed_nblast_mat <- function(x,
                                       overwrite = c("combine","yes","no"),
                                       file = NULL,
                                       threshold=-0.5,
                                       digits=3,
                                       format=c("rda", "rds"),
                                       remove = NULL,
                                       ...) {
  format=match.arg(format)
  overwrite=match.arg(overwrite)
  objname=deparse(substitute(x))
  newobjname <- paste0(objname, ".compressed")
  fname <- paste0(file.path(file, newobjname), ".", format)
  combine = FALSE
  if(file.exists(fname)){
    if(overwrite=="no"){
      stop(fname, " already exists")
    }else if(overwrite=="combine"){
      combine = TRUE
      if(format=="rds"){
        old = tryCatch(readRDS(fname), function(e) NULL)
      }else{
        old = tryCatch(load_assign(fname), function(e) NULL)
      }
    }
  }
  if(is.null(old)){
    warning("No extant NBLAST at ", fname)
  }
  colnames(x) = gsub("_m$","",colnames(x)) # factor in mirrored hemibrain neurons
  rownames(x) = gsub("_m$","",rownames(x))
  x = apply(x, 2, function(i) tapply(i, rownames(x), sum, na.rm = TRUE))
  x = t(apply(t(x), 2, function(i) tapply(i, colnames(x), sum, na.rm = TRUE)))
  x[x<threshold]=threshold
  x=round(x, digits=digits)
  y = x
  if(combine){
    old = old[!rownames(old)%in%rownames(x),]
    if(!is.null(remove)){
      old = old[!rownames(old)%in%remove,!colnames(old)%in%remove]
    }
    if(!is.null(old)){
      if(nrow(old)){
        warning("combining with extant file: ", fname)
        y = plyr::rbind.fill.matrix(old, x)
        rownames(y) = c(rownames(old), rownames(x))
      }else{
        warning("original NBLAST depleted")
        y = x
      }
    }
  }
  message("Resaving a compressed version of ", objname, " to ", fname)
  if(format=="rds") {
    saveRDS(y, file=fname, ...)
  } else {
    # the object being saved must have the name that you would like it
    # to have when it is loaded again
    assign(newobjname, y)
    save(list = newobjname, file = fname, ...)
  }
}

# hidden
load_assign <- function(f){
  env <- new.env()
  nm <- load(f, env)[1]
  env[[nm]]
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
remove_unused_filehash <- function(path){
  for(p in path){
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
        delete = setdiff(all.fh,all.keys)
        message("Deleting ", length(delete), " files")
        delete = file.path(data,delete)
        file.remove(delete)
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
  big = setdiff(obj.files,ids)
  if(length(big)){
    warning("Dropping ", length(big), " .obj files larger than ", max.file.size, " bytes")
  }
  batches = split(ids, round(seq(from = 1, to = numCores, length.out = length(ids))))
  batch = 0
  foreach.skeletons <- foreach::foreach (batch = seq_along(batches)) %dopar% {
    neuron.ids = batches[[batch]]
    j = tryCatch({
      skels = fafbseg::skeletor(neuron.ids, save.obj = NULL, mesh3d = FALSE, ...)
      skels[,"id"] = names(skels) = basename(gsub("\\.obj","",names(skels)))
      nat::write.neurons(skels, dir=swc, format='swc', Force = FALSE)
      skels},
      error = function(e){
        cat(as.character(e))
        message(paste(neuron.ids,collapse=","))
        NULL
      })
  }
  isnl = sapply(foreach.skeletons, nat::is.neuronlist)
  do.call(c, foreach.skeletons[isnl])
}

# hidden
download_neuron_obj_batch <- function(ids, numCores = 1, ratio = 1, save.obj = "obj"){
  batches = split(ids, round(seq(from = 1, to = numCores, length.out = length(ids))))
  batch = 0
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
  if(file.exists(rds)){
    old.neurons = nat::read.neuronlistfh(rds)
    if(!is.null(attr(old.neurons,"df"))){
      old.neurons = old.neurons[!sapply(old.neurons, function(x) isFALSE(x))]
      old.neurons = old.neurons[setdiff(names(old.neurons),names(x))]
      x = nat::union(old.neurons, x)
    }
  }
  if(nat::is.neuronlist(x)){
    given.neurons = nat::as.neuronlistfh(x, dbdir= data, WriteObjects = "missing", ...)
    nat::write.neuronlistfh(given.neurons, file=rds, overwrite=TRUE, ...)
  }else{
    warning("Could not create neuronlistfh object")
  }
}


# hidden
correct_id <-function(v){
  gsub(" ","",v)
}


