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
  x = apply(x, 2, function(i) tapply(i, rownames(x), sum, na.rm = TRUE))
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


# Obtains flywire IDs from flywire posirions on google sheets
# The sheets must either have columns fw.x, fw.y, fw.z
# or one column, flywire.xyz with number separated by commas, in the form: x,y,z
flywire_ids_update <- function(selected_sheets = options()$hemibrainr_gsheets,
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
    tabs = gsheet_manipulation(FUN = googlesheets4::sheet_names,
                               ss = selected_sheet,
                               return = TRUE)
    for(tab in tabs){
      gs.t = gsheet_manipulation(FUN = googlesheets4::read_sheet,
                                              ss = selected_sheet,
                                              sheet = tab,
                                              guess_max = 3000,
                                              return = TRUE)
      if(ncol(gs.t)&&sum(grepl("fw.x|flywire.xyz",colnames(gs.t)))>0){
        # Separate x,y,z positions
        gs1 = subset(gs.t, is.na(gs.t$flywire.xyz))
        gs2 = subset(gs.t, !is.na(gs.t$flywire.xyz))
        gs1$flywire.xyz = apply(gs1[,c("fw.x","fw.y",'fw.z')],1,paste,sep=";",collapse=";")
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
        # batch = 1
        # batches = split(1:nrow(gs.t), round(seq(from = 1, to = numCores, length.out = nrow(gs.t))))
        # foreach.ids <- foreach::foreach (batch = 1:length(batches)) %dopar% {
        #   pos = gs.t[batches[[batch]],]
        #   j <- tryCatch({i = fafbseg::flywire_xyz2id(pos[,c("fw.x","fw.y",'fw.z')], rawcoords = TRUE)
        #   i[is.na(i)|is.nan(i)] = 0
        #   names(i) = pos$flywire.xyz
        #   # if(sum(i==0)>0){
        #   #   i[i==0] =  fafbseg::flywire_xyz2id(pos[i==0,c("fw.x","fw.y",'fw.z')], rawcoords = FALSE)
        #   # }
        #   i}, error = function(e) NULL)
        # }
        foreach.ids = fafbseg::flywire_xyz2id(gs.t[,c("fw.x","fw.y",'fw.z')], rawcoords = TRUE)
        names(foreach.ids) = gs.t[,"flywire.xyz"]
        fids = unlist(foreach.ids)
        fids[is.na(fids)|is.nan(fids)] = "0"
        gs.t[match(names(fids),gs.t$flywire.xyz),"flywire.id"] = fids
        # Update
        rownames(gs.t) = NULL
        gs.t$flywire.xyz = apply(gs.t[,c("fw.x","fw.y",'fw.z')],1,paste,sep=";",collapse=";")
        googlesheets4::write_sheet(gs.t[0,],
                                   ss = selected_sheet,
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
  xyz.good = tryCatch(nat::xyzmatrix(conns), error = function(e) NULL)
  if(!is.null(xyz.good)){
    conns.t = nat.templatebrains::xform_brain(nat::xyzmatrix(conns), reference = reference, sample = sample, method = method, progress.rjava=progress.rjava, ...)
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

# hidden
## Save with given name
saveit <- function(..., file) {
  x <- list(...)
  save(list=names(x), file=file, envir=list2env(x))
}

