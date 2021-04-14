# hidden
# nb = nblast_big(kcs20,kcs20[1:10],numCores = 4)
`%dopar%` <- foreach::`%dopar%`
`%:%` <- foreach::`%:%`
nblast_big <-function(query.neuronlistfh, target.neuronlistfh,
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
                      digits = 3,
                      update.old = NULL){

  # Register cores
  if(!requireNamespace('nat.nblast'))
    stop("Please install suggested nat.nblast package!")
  batch.size = numCores #floor(sqrt(numCores))
  cl = parallel::makeForkCluster(batch.size)
  check_package_available("doParallel")
  doParallel::registerDoParallel(cl)

  # What are our query and target neurons
  query = intersect(names(query.neuronlistfh),query)
  if(!length(query)){
    stop("query neurons not in query.neuronlistfh")
  }
  target = names(target.neuronlistfh)

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
        nblast.mat = bigstatsr::FBM(length(target),length(query))
        nblast.mat[match(rownames(old), target),match(colnames(old),query)] = old
      }else{
        nblast.mat = bigstatsr::FBM(length(target),length(query))
      }
    }else{
      nblast.mat = bigstatsr::FBM(length(target),length(query))
    }
  }else{
    nblast.mat = bigstatsr::FBM(length(target),length(query))
  }

  # Get batches to iterate over
  ## this would be a better way of doing it, but at the moment thwarted by DB1 lock files
  # batches.query = split(sample(query), round(seq(from = 1, to = batch.size, length.out = length(query))))
  # batches.target = split(sample(target), round(seq(from = 1, to = batch.size, length.out = length(target))))
  chosen.query = chosen.target = NULL

  # Foreach loop
  ## this would be a better way of doing it, but at the moment thwarted by DB1 lock files
  # by.query <- foreach::foreach (chosen.query = batches.query, .combine = 'c') %:%
  #   foreach::foreach (chosen.target = batches.target, .combine = 'c') %dopar% {
  ### This is a slightly more inefficient way
  batches.query = split(sample(query.neuronlistfh[query]), round(seq(from = 1, to = batch.size, length.out = length(query))))
  batches.target = split(sample(target.neuronlistfh), round(seq(from = 1, to = batch.size, length.out = length(target.neuronlistfh))))
  by.query <- foreach::foreach (query.neuronlist = batches.query, .combine = 'c', .errorhandling='pass') %:%
    foreach::foreach (target.neuronlist = batches.target, .combine = 'c', .errorhandling='pass') %dopar% {
      ## this would be a better way of doing it, but at the moment thwarted by DB1 lock files
      # query.neuronlist = query.neuronlistfh[chosen.query]
      # target.neuronlist = target.neuronlistfh[chosen.target]
      if(!is.null(query.neuronlist)&&length(query.neuronlist)){
        ### This is a slightly more inefficient way
        query.neuronlist = query.neuronlist[unlist(sapply(query.neuronlist, is_big_dps,no.points=no.points))]
        target.neuronlist = target.neuronlist[unlist(sapply(target.neuronlist, is_big_dps,no.points=no.points))]
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
        nblast.res.native = (nblast.res.1+t(nblast.res.2))/2
        ### NBLAST mirrored
        if(!is.null(query.addition.neuronlist)&&!length(query.addition.neuronlist)){
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
        }else{
          nblast.res.sub = nblast.res.native
        }
        # Compress
        if(compress){
          nblast.res.sub[nblast.res.sub<threshold] = threshold
          nblast.res.sub = round(nblast.res.sub, digits=digits)
        }
        try({nblast.mat[match(chosen.target, target),match(chosen.query, query)] = nblast.res.sub}, silent = FALSE)
        NULL
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
