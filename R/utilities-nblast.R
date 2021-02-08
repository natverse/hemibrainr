# hidden
# nb = nblast_big(kcs20,kcs20[1:10],numCores = 4)
`%dopar%` <- foreach::`%dopar%`
`%:%` <- foreach::`%:%`
nblast_big <-function(query.neuronlistfh, target.neuronlistfh,
                      query.neuronlistfh.addition = NULL,
                      query = names(query.neuronlistfh),
                      numCores=1,
                      smat = NULL,
                      sd = 3,
                      version = c(2, 1),
                      normalised = TRUE,
                      UseAlpha = TRUE,
                      no.points = 2){

  # Register cores
  cl = parallel::makeForkCluster(numCores)
  doParallel::registerDoParallel(cl)

  # What are our query and target neurons
  query = intersect(names(query.neuronlistfh),query)
  if(!length(query)){
    stop("query neurons not in query.neuronlistfh")
  }
  target = names(target.neuronlistfh)

  # Make matrix to fill
  nblast.mat = bigstatsr::FBM(length(target),length(query))
  batch.size = numCores #floor(sqrt(numCores))

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
  by.query <- foreach::foreach (query.neuronlist = batches.query, .combine = 'c') %:%
    foreach::foreach (target.neuronlist = batches.target, .combine = 'c') %dopar% {
      ## this would be a better way of doing it, but at the moment thwarted by DB1 lock files
      # query.neuronlist = query.neuronlistfh[chosen.query]
      # target.neuronlist = target.neuronlistfh[chosen.target]
      if(is.null(query.neuronlist)||!length(query.neuronlist)){
        return(NULL)
      }
      ### This is a slightly more inefficient way
      query.neuronlist = query.neuronlist[unlist(sapply(query.neuronlist,hemibrainr:::is_big_dps,no.points=no.points))]
      target.neuronlist = target.neuronlist[unlist(sapply(target.neuronlist,hemibrainr:::is_big_dps,no.points=no.points))]
      chosen.query = names(query.neuronlist)
      chosen.target = names(target.neuronlist)
      if(!length(chosen.query)||!length(chosen.target)){
        return(NULL)
      }

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
      if(!is.null(query.neuronlistfh.addition)){
        query.addition.neuronlist = query.neuronlistfh.addition[names(query.neuronlistfh.addition)%in%chosen.query]
        if(length(query.addition.neuronlist)!=length(query.neuronlist)){
          break
        }
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
      nblast.mat[match(chosen.target, target),match(chosen.query,query)] = nblast.res.sub
      NULL
    }
  parallel::stopCluster(cl)
  nmat = nblast.mat[,]
  dimnames(nmat) = list(target, query)
  nmat
}

# hidden
is_big_dps <- function(dps, no.points = 5){
  !is.null(dps$points)&&"matrix"%in%class(dps$points)&&nrow(dps$points)>=no.points
}
