# hidden
# nb = nblast_big(kcs20,kcs20[1:10],numCores = 4)
`%dopar%` <- foreach::`%dopar%`
`%:%` <- foreach::`%:%`
nblast_big <-function(query.neuronlistfh,
                      target.neuronlistfh,
                      query.addition.neuronlistfh = NULL,
                      query = names(query.neuronlistfh),
                      numCores=1,
                      smat = NULL,
                      sd = 3,
                      version = c(2, 1),
                      normalised = TRUE,
                      UseAlpha = TRUE,
                      no.points = 2,
                      compress = TRUE,
                      threshold = 0, # or -0.5?
                      digits = 6,
                      update.old = NULL,
                      outfile = "",
                      overlap = FALSE,
                      split = FALSE,
                      check_function = is_big,
                      just.leaves = TRUE,
                      normalise = FALSE
                      ){

  # Register cores
  check_package_available("nat.nblast")
  check_package_available("doParallel")
  check_package_available("doSNOW")
  batch.size = numCores #floor(sqrt(numCores))
  # cl = parallel::makeForkCluster(batch.size)
  cl = parallel::makeCluster(batch.size, outfile = outfile)
  doParallel::registerDoParallel(cl)
  doSNOW::registerDoSNOW(cl)
  if(numCores<2){
    `%go%` <- foreach::`%do%`
  }else{
    `%go%` <- foreach::`%dopar%`
  }

  # What are our query and target neurons
  query = intersect(names(query.neuronlistfh),query)
  if(!length(query)){
    stop("query neurons not in query.neuronlistfh")
  }
  target = names(target.neuronlistfh)

  # Split
  if(split){
    target.axons = paste0(target, "_axon")
    target.dendrites = paste0(target, "_dendrites")
    query.axons = paste0(query, "_axon")
    query.dendrites = paste0(query, "_dendrites")
    target.names = c(target.axons, target.dendrites)
    query.names = c(query.axons, query.dendrites)
  }else{
    query.names = query
    target.names = target
  }

  # Initialize matrix
  nblast.mat = bigstatsr::FBM(length(target.names),length(query.names), init = NA)
  if(is.null(update.old)){
    update.old = NULL
  }else  if(!file.exists(update.old)){
    update.old = NULL
  }

  # Make matrix to fill
  ## Get old matrix
  if(!is.null(update.old)){
    if(grepl("rds$",update.old)){
      old = t(readRDS(update.old))
    }else{
      old = t(load_assign(update.old))  # historically, I had the matrix the other way around :(
    }
    old = old[apply(old, 1, function(r) sum(is.na(r))==0),apply(old, 2, function(c) sum(is.na(c))==0)]
    old = old[!colnames(old)%in%query,!rownames(old)%in%target]
    if(length(old)){
      if(nrow(old)&ncol(old)){
        query.names = c(sort(colnames(old)), setdiff(query.names,colnames(old)))
        target.names = c(sort(rownames(old)), setdiff(target.names,rownames(old)))
        nblast.mat = bigstatsr::FBM(length(target.names),length(query.names), init = NA)
        nblast.mat[match(rownames(old), target.names),match(colnames(old),query.names)] = old
      }else{
        nblast.mat = bigstatsr::FBM(length(target.names),length(query.names), init = NA)
      }
    }else{
      nblast.mat = bigstatsr::FBM(length(target.names),length(query.names), init = NA)
    }
  }else{
    nblast.mat = bigstatsr::FBM(length(target.names),length(query.names), init = NA)
  }

  # Get batches to iterate over
  ## this would be a better way of doing it, but at the moment thwarted by DB1 lock files
  # batches.query = split(sample(query), round(seq(from = 1, to = batch.size, length.out = length(query))))
  # batches.target = split(sample(target), round(seq(from = 1, to = batch.size, length.out = length(target))))

  # Foreach loop
  ## this would be a better way of doing it, but at the moment thwarted by DB1 lock files
  # by.query <- foreach::foreach (chosen.query = batches.query, .combine = 'c') %:%
  #   foreach::foreach (chosen.target = batches.target, .combine = 'c') %dopar% {
  ### This is a slightly more inefficient way
  batches.query = split(sample(query), round(seq(from = 1, to = batch.size, length.out = length(query))))
  batches.target = split(sample(target), round(seq(from = 1, to = batch.size, length.out = length(target))))

  # Progress bar
  iterations <- length(batches.query)
  pb <- txtProgressBar(max = iterations, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)

  # Parallel process
  by.query <- foreach::foreach(chosen.query = batches.query, .combine = 'c', .errorhandling='pass', .options.snow = opts) %do%
     foreach::foreach(chosen.target = batches.target, .combine = 'c', .errorhandling='pass') %go% {
      ## this would be a better way of doing it, but at the moment thwarted by DB1 lock files
      query.neuronlist = query.neuronlistfh[names(query.neuronlistfh)%in%unlist(chosen.query)]
      target.neuronlist = target.neuronlistfh[names(target.neuronlistfh)%in%unlist(chosen.target)]
      chosen.query = names(query.neuronlist)
      chosen.target = names(target.neuronlist)
      if(!length(query.neuronlist)|!length(target.neuronlist)){
        message("chosen not in data")
        next
      }
      query.addition.neuronlist = NULL
      if(!is.null(query.neuronlist)&&length(query.neuronlist)){
        ### This is a slightly more inefficient way
        if(!is.null(check_function)){
          query.neuronlist = query.neuronlist[unlist(sapply(query.neuronlist, check_function , no.points))]
          target.neuronlist = target.neuronlist[unlist(sapply(target.neuronlist, check_function, no.points))]
        }
        if(!is.null(query.addition.neuronlistfh)){
          query.addition.neuronlist = query.addition.neuronlistfh[names(query.addition.neuronlistfh)%in%chosen.query]
        }else{
          query.addition.neuronlist = NULL
        }
        chosen.query = union(names(query.neuronlist),names(query.addition.neuronlist))
        if(!is.null(update.old)){ # Do not re-calculate calculated scores for the most part
          if(all(chosen.query%in%colnames(old))){
            target.neuronlist = target.neuronlist[setdiff(names(target.neuronlist),rownames(old))]
          }
        }
        chosen.target = names(target.neuronlist)
      }

      # Run NBLASTs
      if(length(chosen.query)&&length(chosen.target)){
        ### NBLAST native
        if(!overlap){
          nblast.res.1 = nat.nblast::nblast(query = query.neuronlist,
                                            target = target.neuronlist,
                                            .parallel=FALSE,
                                            normalised = normalised,
                                            smat = smat,
                                            sd = sd,
                                            version = version,
                                            UseAlpha = UseAlpha,
                                            OmitFailures = FALSE)
          nblast.res.2 = nat.nblast::nblast(query = target.neuronlist,
                                            target = query.neuronlist,
                                            .parallel=FALSE,
                                            normalised = normalised,
                                            smat = smat,
                                            sd = sd,
                                            version = version,
                                            UseAlpha = UseAlpha,
                                            OmitFailures = FALSE)
          if(is.null(dim(nblast.res.2))){
            nblast.res.2 = matrix(nblast.res.2, nrow = 1, ncol = length(nblast.res.2), dimnames = list(chosen.query,chosen.target))
          }
          nblast.res.sub = nblast.res.native = tryCatch((nblast.res.1+t(nblast.res.2))/2, error = function(e){
            message(as.character(e))
            nblast.res.1
          })
          ### NBLAST mirrored
          if(!is.null(query.addition.neuronlist)&&length(query.addition.neuronlist)){ # For an additional neuronlist of mirrored neurons
            nblast.res.3 = nat.nblast::nblast(query = query.addition.neuronlist,
                                              target = target.neuronlist,
                                              .parallel=FALSE,
                                              normalised = normalised,
                                              smat = smat,
                                              sd = sd,
                                              version = version,
                                              UseAlpha = UseAlpha,
                                              OmitFailures = FALSE)
            nblast.res.4 = nat.nblast::nblast(query = target.neuronlist,
                                              target = query.addition.neuronlist,
                                              .parallel=FALSE,
                                              normalised = normalised,
                                              smat = smat,
                                              sd = sd,
                                              version = version,
                                              UseAlpha = UseAlpha,
                                              OmitFailures = FALSE)
            if(is.null(dim(nblast.res.2))){
              nblast.res.4 = matrix(nblast.res.4, nrow = 1, ncol = length(nblast.res.2), dimnames = list(names(query.addition.neuronlist),chosen.target))
            }
            nblast.res.m = (nblast.res.3+t(nblast.res.4))/2
            nblast.res.sub = plyr::rbind.fill.matrix(t(nblast.res.native), t(nblast.res.m))
            rownames(nblast.res.sub) = c(colnames(nblast.res.native), colnames(nblast.res.sub))
            nblast.res.sub = collapse_matrix_by_names(nblast.res.sub, FUN = max)
          }
        }else{
          if(split){
            q.axons = axonic_cable(query.neuronlist)
            names(q.axons) = paste0(names(q.axons),"_axon")
            q.dendrites = dendritic_cable(query.neuronlist)
            names(q.dendrites) = paste0(names(q.dendrites),"_dendrite")
            t.axons = axonic_cable(target.neuronlist)
            names(t.axons) = paste0(names(t.axons),"_axon")
            t.dendrites = dendritic_cable(target.neuronlist)
            names(t.dendrites) = paste0(names(t.dendrites),"_dendrite")
            query.neuronlist = c(q.axons, q.dendrites)
            target.neuronlist = c(t.axons, t.dendrites)
            chosen.query = names(query.neuronlist)
            chosen.query = names(target.neuronlist)
          }
          nblast.res.sub = overlap_score_delta(query.neuronlist, target.neuronlist, just.leaves=just.leaves, normalise=normalise)
          nblast.res.sub = t(nblast.res.sub)
        }
        # Compress
        if(compress){
          nblast.res.sub[nblast.res.sub<threshold] = threshold
          nblast.res.sub = signif(nblast.res.sub, digits=digits)
        }
        nblast.mat[match(unlist(chosen.target), target),match(unlist(chosen.query), query)] = nblast.res.sub
        #try({nblast.mat[match(chosen.target, target),match(chosen.query, query)] = nblast.res.sub}, silent = FALSE)
        #nblast.res.sub
        NULL
      }
  }
  for(i in 1:length(by.query)){
    if(!is.null(by.query[[i]])){
      message(by.query[[i]])
    }
  }
  parallel::stopCluster(cl)
  clear = gc()
  nmat = matrix(nblast.mat[,], nrow = length(target), ncol = length(query))
  dimnames(nmat) = list(target, query)
  nmat
}

# hidden
is_big_dps <- function(dps, no.points = 5){
  !is.null(dps$points)&&"matrix"%in%class(dps$points)&&nrow(dps$points)>=no.points
}
is_big <- function(dps, no.points = 5){
  xyz <- nat::xyzmatrix(dps)
  !is.null(xyz)&&nrow(xyz)>=no.points
}

# We can make another version of the matrix that contains essentially all
# the useful information, but compresses much smaller (~16x at time of writing)
save_compressed_nblast_mat <- function(x,
                                       overwrite = c("combine","yes","no"),
                                       file = NULL,
                                       threshold= 0, # or -0.5?
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
  old <- NULL
  if(file.exists(fname)){
    if(overwrite=="no"){
      stop(fname, " already exists")
    }else if(overwrite=="combine"){
      combine = TRUE
      if(format=="rds"){
        old = tryCatch(readRDS(fname), error = function(e) NULL)
      }else{
        old = tryCatch(load_assign(fname), error = function(e) NULL)
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
  message("Saving a compressed version of ", objname, " to ", fname)
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


# hidden, similar function now in nat
overlap_score_big <- function(output.neurons,
                              input.neurons,
                              delta = 62.5,
                              just.leaves = TRUE,
                              max = exp(-delta^2/(2*delta^2)),
                              normalise = TRUE,
                              update.old = NULL){

  # Register cores
  check_package_available("nat.nblast")
  check_package_available("doParallel")
  batch.size = numCores #floor(sqrt(numCores))
  # cl = parallel::makeForkCluster(batch.size)
  cl = parallel::makeCluster(batch.size, outfile = outfile)
  doParallel::registerDoParallel(cl)
  if(numCores<2){
    `%go%` <- foreach::`%do%`
  }else{
    `%go%` <- foreach::`%dopar%`
  }

  # Make matrix to fill
  ## Get old matrix
  if(!is.null(update.old)){
    if(grepl("rds$",update.old)){
      old = t(readRDS(update.old))
    }else{
      old = t(load_assign(update.old))  # historically, I had the matrix the other way around :(
    }
    old = old[apply(old, 1, function(r) sum(is.na(r))==0),apply(old, 2, function(c) sum(is.na(c))==0)]
    old = old[!colnames(old)%in%query,!rownames(old)%in%target]
    if(length(old)){
      if(nrow(old)&ncol(old)){
        query = c(sort(colnames(old)), setdiff(query,colnames(old)))
        target = c(sort(rownames(old)), setdiff(target,rownames(old)))
        overlap.mat = bigstatsr::FBM(length(target),length(query), init = NA)
        overlap.mat[match(rownames(old), target),match(colnames(old),query)] = old
      }else{
        overlap.mat = bigstatsr::FBM(length(target),length(query), init = NA)
      }
    }else{
      overlap.mat = bigstatsr::FBM(length(target),length(query), init = NA)
    }
  }else{
    overlap.mat = bigstatsr::FBM(length(target),length(query), init = NA)
  }

  # Sort neurons
  output.neurons = nat::as.neuronlist(output.neurons)
  input.neurons = nat::as.neuronlist(input.neurons)
  score.matrix = matrix(0,nrow = length(output.neurons), ncol = length(input.neurons))
  rownames(score.matrix) = names(output.neurons)
  colnames(score.matrix) = names(input.neurons)
  if(just.leaves){
    input.neurons.d = nat::nlapply(input.neurons, function(x) nat::xyzmatrix(x)[nat::endpoints(x),], .progress = "none")
  }else{
    input.neurons.d = nat::nlapply(input.neurons, nat::xyzmatrix, .progress = "none")
  }
  for (n in 1:length(output.neurons)){
    if(just.leaves){
      a = nat::xyzmatrix(output.neurons[[n]])[nat::endpoints(output.neurons[[n]]),]
    }else{
      a = nat::xyzmatrix(output.neurons[[n]])
    }
    if(normalise){
      s = sapply(input.neurons.d, function(x)lengthnorm(maxout(exp(-nabor::knn(query = a, data = x,k=nrow(x))$nn.dists^2/(2*delta^2)),max=max)))
    }else{
      s = sapply(input.neurons.d, function(x)sum(maxout(exp(-nabor::knn(query = a, data = x,k=nrow(x))$nn.dists^2/(2*delta^2)),max=max))) # Score similar to that in Schlegel et al. 2015
    }
    score.matrix[n,] = s
  }
  score.matrix
}


hemibrain_flywire_match_check <-function(matchcheck = "~/Downloads/"){

  # Colors
  nt.cols = c(
    gaba = "#E6A749",
    acetylcholine = "#4B506B",
    glutamate = "#70B657",
    octopamine = "#7A4F98",
    serotonin = "#93A3CF",
    dopamine = "#CF6F6C",
    neither = "grey70",
    unknown = "grey10",
    hemibrain = hemibrain_bright_colors[["marine"]],
    flywire = hemibrain_bright_colors[["cerise"]],
    hemisphere = hemibrain_bright_colors[["green"]])

  # Get flywire meta data
  fw.meta = flywire_meta()
  matches = flytable_matches()

  # Get matches
  hb.ids = unique(matches$hemibrain_match)

  # Get hemibrain remote data
  hb = hemibrain_neurons()
  hb = hb[names(hb) %in% hb.ids]
  hb = hemibrain_add_nt.neuronlist(hb)
  hb.meta = hb[,]
  matches$cell_type = hb.meta$type[match(matches$hemibrain_match,hb.meta$bodyid)]
  exclude = c("none","unknown","NA", "", " ")
  matches.lr = subset(matches, ! hemisphere_match %in% exclude & !is.na(hemisphere_match) )
  matches.hb = subset(matches, !is.na(hemibrain_match) & !hemibrain_match %in% exclude)

  # Get nblast scores
  nb.flywire.mirrored = hemibrain_nblast("flywire-mirror")
  nb.flywire.hemibrain = hemibrain_nblast("hemibrain-flywire")

  # NBLAST LR matches
  hits = unlist(apply(nb.flywire.mirrored[matches.lr$root_id,], 1, function(x) colnames(nb.flywire.mirrored)[which.max(x)]))
  hits.value = unlist(apply(nb.flywire.mirrored[matches.lr$root_id,], 1, function(x) max(x)))
  match.value = unlist(sapply(seq_along(matches.lr$root_id), function(x)  nullToNA(nb.flywire.mirrored[matches.lr$root_id[x],matches.lr$hemisphere_match[x]])))
  matches.lr$hemisphere_nblast = match.value

  # NBLAST hemibrain_matches
  hits.hb = unlist(apply(nb.flywire.hemibrain[matches.hb$hemibrain_match,], 1, function(x) colnames(nb.flywire.hemibrain)[which.max(x)]))
  hits.hb.value = unlist(nb.flywire.hemibrain[matches.hb$hemibrain_match,], 1, function(x) max(x))
  match.value = unlist(sapply(seq_along(matches.hb$hemibrain_match), function(x) nullToNA(nb.flywire.hemibrain[matches.hb$root_id[x],matches.hb$hemibrain_match[x]]) ))
  matches.hb$hemibrain_nblast = match.value

  # Add nt information
  matches.hb$hb_top_nt = hb.meta$top_nt[match(matches.hb$hemibrain_match, hb.meta$bodyid)]
  matches.hb$fw_top_nt = fw.meta$top_nt[match(matches.hb$root_id, fw.meta$root_id)]
  matches.lr$n_top_nt = fw.meta$top_nt[match(matches.lr$root_id, fw.meta$root_id)]
  matches.lr$m_top_nt = fw.meta$top_nt[match(matches.lr$hemisphere_match, fw.meta$root_id)]

  # Create directories
  matchcheck.lr = file.path(matchcheck,"lr_match_check")
  matchcheck.hb = file.path(matchcheck,"hemibrain_match_check")
  dir.create(matchcheck.lr, showWarnings = FALSE)
  dir.create(matchcheck.hb, showWarnings = FALSE)

  # Plot for flywire, local
  fafb14_view()
  for (i in 1:nrow(matches.lr)){
    try({
      a = matches.lr[i,"root_id"]
      a.nt = matches.lr[i,"n_top_nt"]
      b = matches.lr[i,"hemisphere_match"]
      b.nt = matches.lr[i,"m_top_nt"]
      v = matches.lr[i,"hemisphere_nblast"]
      v = ifelse(is.na(v),NA,round(v, 4))
      an = read_cloudvolume_meshes(a)
      bn = read_cloudvolume_meshes(b)
      plot3d(elmr::FAFB14, alpha = 0.1)
      plot3d(an, col = nt.cols[["flywire"]])
      plot3d(bn, col = nt.cols[["hemisphere"]])
      par3d(cex=3.0)
      a.nt = ifelse(is.na(a.nt),"unknown",a.nt)
      b.nt = ifelse(is.na(b.nt),"unknown",b.nt)
      rgl::text3d(x = 350000, y = 400000, z = 100000, texts = paste0("flywire: ", a.nt), col = nt.cols[[a.nt]])
      rgl::text3d(x = 700000, y = 400000, z = 100000, texts = paste0("hemisphere: ", b.nt), col = nt.cols[[b.nt]])
      shot = file.path(matchcheck.lr, sprintf("%s_%s_%s_flywire_LR_%s_%s.png",a.nt, b.nt, v,a,b))
      rgl::rgl.snapshot(shot)
    })
    print(i)
    rgl::clear3d()
  }

  # Plot for hemibrain, local
  fafb14_view()
  hb.fafb14 = hemibrain_neurons(brain = "FAFB14")
  for (i in 1:nrow(matches.hb)){
    try({
      ct = matches.hb[i,"cell_type"]
      a = matches.hb[i,"root_id"]
      a.nt = matches.hb[i,"fw_top_nt"]
      b = matches.hb[i,"hemibrain_match"]
      b.nt = matches.hb[i,"hb_top_nt"]
      v = matches.hb[i,"hemibrain_nblast"]
      v = ifelse(is.na(v),NA,round(v, 4))
      an = read_cloudvolume_meshes(a)
      bn = hb.fafb14[names(hb.fafb14)%in%b]
      plot3d(elmr::FAFB14, alpha = 0.1)
      plot3d(an, col = nt.cols[["flywire"]])
      plot3d(bn, col = nt.cols[["hemibrain"]])
      rgl::par3d(cex=3.0)
      a.nt = ifelse(is.na(a.nt),"unknown",a.nt)
      b.nt = ifelse(is.na(b.nt),"unknown",b.nt)
      rgl::text3d(x = 350000, y = 400000, z = 100000, texts = paste0("flywire: ", a.nt), col = nt.cols[[a.nt]])
      rgl::text3d(x = 700000, y = 400000, z = 100000, texts = paste0("hemibrain: ", b.nt), col = nt.cols[[b.nt]])
      shot = file.path(matchcheck.hb, sprintf("%s_%s_%s_%s_flywire_hemibrain_%s_%s.png", a.nt, b.nt, ct, v, a,b))
      rgl::rgl.snapshot(shot)
    })
    print(i)
    rgl::clear3d()
  }

}

# flywire view
fafb14_view <- function(){
  open3d(userMatrix = structure(c(0.994450092315674, 0.0210675243288279,
                                  -0.103079274296761, 0, 0.0466670729219913, -0.966418564319611,
                                  0.252699792385101, 0, -0.0942939221858978, -0.256107956171036,
                                  -0.962038159370422, 0, 0, 0, 0, 1), .Dim = c(4L, 4L)), zoom = 0.644609212875366,
         windowRect = c(1440L, 45L, 3790L, 1416L))
}

