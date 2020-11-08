#####################################################################################
################################ FIB-FAFB Matching ##################################
#####################################################################################

#' Match up neurons between the hemibrain and FAFB
#'
#' @description Match up neurons between the hemibrain, FAFB and light level data and save the
#'   result using a
#'   \href{https://docs.google.com/spreadsheets/d/1OSlDtnR3B1LiB5cwI5x5Ql6LkZd8JOS5bBr-HTi0pOw/edit#gid=0}{Google
#'   Sheet} on the hemibrain Google Team Drive operated by the flyconnectome
#'   group at the University of Cambridge. You must have access to the Team Drive
#'   in order to use this function. This function makes use of the Google
#'   Filestream application, which should be installed on your machine. Further,
#'   note that neurons are read from the FAFB CATMAID project, and you must have
#'   login details for this project recorded in your .Renviron for these
#'   functions to work.
#'
#' @param ids body IDs for hemibrain neurons present in the
#'   \href{https://docs.google.com/spreadsheets/d/1OSlDtnR3B1LiB5cwI5x5Ql6LkZd8JOS5bBr-HTi0pOw/edit#gid=0}{Google
#'    Sheet}, for which the user will attempt to make a match if one has not
#'   been made already. Else, LM neuron IDs in the tab 'lm' when using \code{lm_matching}.
#' @param hemibrain.nblast a FAFB (rows) - hemibrain (columns) normalised
#'   NBLAST matrix. By default this is read from the flyconnectome Team Drive.
#' @param mirror.nblast a flywire (rows) - flywire-mirrored (columns) normalised
#' NBLAST matrix. By default this is read from the hemibrain Team Drive.
#' @param selected_file the Google Sheet database to read and write from. For
#'   now, defaults to this
#'   \href{https://docs.google.com/spreadsheets/d/1OSlDtnR3B1LiB5cwI5x5Ql6LkZd8JOS5bBr-HTi0pOw/edit#gid=0}{Google
#'    Sheet}.
#' @param batch_size the number of FAFB top matches to read from CATMAID in one
#'   go.
#' @param db Either a neuronlist or the name of a character vector naming a
#'   neuronlist. Defaults to the value of \code{\link{hemibrain_neurons}()}.
#' @param repository whether to match up FAFB skeletons from
#' CATMAID (using \code{catmaid::read.neurons.catmaid}) or flywire (using \code{flywire_neurons}) for matching. Alternatively,light level
#' skeletons  (\code{"lm"}) from Dolan et al. and Frechter et al. 2019 (eLife), stored in the package \code{lhns} as
#' \code{most.lhns}.
#' @param query a neuronlist of light level neurons to match against. Should correspond to the given NBLAST matrix.
#' Defaults to reading a transformed \code{most.lhns} from the Hemibrain Google Team Drive.
#' @param overwrite Whether or not to overwrite matches already made.
#' The neurons you could possibly be looking at are selected through
#' the arguments: \code{ids}, \code{column}, \code{field}. If \code{ids} is \code{NULL} you will be given all neurons that have been assigned to your user on the
#' \href{https://docs.google.com/spreadsheets/d/1OSlDtnR3B1LiB5cwI5x5Ql6LkZd8JOS5bBr-HTi0pOw/edit#gid=0}{Google sheet}.
#' If \code{overwrite} is set to \code{FALSE}, you will not overwrite any matches that have already been made from among the selected neurons.
#' If \code{'mine'} you will re-examine and overwrite any matches you have made out of the selected neurons.
#' With \code{'mine_empty'} the same happens, but you will also retain neurons that have no match, and have been assigned to any other user.
#' If \code{FALSE} (be careful!) then you overwrite made matches among the selected neurons.
#' @param column defaults to \code{NULL}, no further subsetting. Else, you can select a column from the Google sheet.
#' Only  neurons with a certain value (\code{field}) in that column will be chosen for matching.
#' @param field defaults to \code{NULL}, no further subsetting. Else, it is a value in \code{column}.
#'
#' @details Currently, the
#'   \href{https://docs.google.com/spreadsheets/d/1OSlDtnR3B1LiB5cwI5x5Ql6LkZd8JOS5bBr-HTi0pOw/edit#gid=0}{Google
#'    Sheet} is set up with  limited number of users, each of whom have been
#'   assigned a number of neurons to match up. In order to add yourself as a
#'   user, simply open this Google Sheet in your browser and add your initials
#'   to neurons of your choosing on the rightmost column 'Users'. Once a match
#'   is recorded, the user selects a quality for that match. There can be no
#'   match (none), a poor match (poor) an okay match (okay) or an exact match (good). As a
#'   rule of thumb, a poor match could be a neuron from a very similar same cell
#'   type or a highly untraced neuron that may be the correct cell type. An okay
#'   match should be a neuron that looks to be from the same morphological cell
#'   type but there may be some discrepancies in its arbour. A good match is a
#'   neuron that corresponds well between FAFB and the hemibrain data.
#'
#' @examples
#' \donttest{
#' \dontrun{
#'
#' # install package to bridge neurons between FAFB14 and hemibrain space
#' if (!requireNamespace("remotes")) install.packages("remotes")
#' remotes::install_github('natverse/nat.jrcbrains')
#' nat.jrcbrains::download_saalfeldlab_registrations()
#'
#' # Load precomputed NBLAST from the flyconnectome Team Drive
#' load(file.path("/Volumes/GoogleDrive/Shared drives/flyconnectome/",
#'   "fafbpipeline/fib.fafb.crossnblast.twigs5.mean.compress.rda"))
#' # Sometimes the filenames are changed by drive
#' # load(file.path("/Volumes/GoogleDrive/Shared drives/flyconnectome",
#' # "fafbpipeline/fib.fafb.crossnblast.twigs5.mean.compress (1).rda"))
#'
#' # Match!
#' hemibrain_matching(hemibrain.nblast = fib.fafb.crossnblast.twigs5.mean.compress)
#' }}
#' @rdname hemibrain_matching
#' @export
#' @seealso \code{\link{hemibrain_adjust_saved_split}}
hemibrain_matching <- function(ids = NULL,
                         hemibrain.nblast = NULL,
                         selected_file = options()$hemibrainr_matching_gsheet,
                         batch_size = 10,
                         db=hemibrain_neurons(), # brain="FAFB"
                         repository = c("CATMAID", "flywire", "lm"),
                         query = NULL,
                         overwrite = c("FALSE","mine","mine_empty","TRUE"),
                         column = NULL,
                         field = NULL){
  repository = match.arg(repository)
  # Other packages
  if(!requireNamespace("nat.jrcbrains", quietly = TRUE)) {
    stop("Please install nat.jrcbrains using:\n", call. = FALSE,
         "remotes::install_github('natverse/nat.jrcbrains')")
  }
  if(!requireNamespace("elmr", quietly = TRUE)) {
    stop("Please install elmr using:\n", call. = FALSE,
         "remotes::install_github('natverse/elmr')")
  }
  # Motivate!
  nat::nopen3d()
  plot_inspirobot()
  unsaved = c()
  message("
          #######################Colours#######################
          black = hemibrain neuron,
          red = potential FAFB matches based on NBLAST score,
          green = a chosen FAFB neuron during scanning,
          blue = your selected FAFB match.
          #######################Colours#######################
          ")
  ## Get NBLAST
  if(is.null(hemibrain.nblast) & repository == "CATMAID"){
    hemibrain.nblast = hemibrain_nblast("hemibrain-fafb14")
  }
  if(is.null(hemibrain.nblast) & repository == "flywire"){
    hemibrain.nblast = t(hemibrain_nblast("hemibrain-flywire"))
  }
  if(is.null(hemibrain.nblast) & repository =="LM"){
    hemibrain.lhns.mean.compressed=NULL
    matname="hemibrain.lhns.mean.compressed"
    message("Loading LM-FIB NBLAST ", matname,
            " from hemibrain Google Team Drive using Google Filestream: ")
    load(sprintf("/Volumes/GoogleDrive/Shared\ drives/hemibrain/hemibrain_nblast/%s.rda", matname))
    hemibrain.nblast = get(matname)
    hemibrain.nblast = t(hemibrain.nblast)
    rm("hemibrain.lhns.mean.compressed")
  }
  # FlyWire or CATMAID?
  if(repository=="FAFB"){
    id = "skid"
  }else{
    id = "flywire.id"
  }
  # Read the Google Sheet
  gs = hemibrain_match_sheet(sheet = "hemibrain", selected_file = selected_file)
  gs$User[is.na(gs$User)] = ""
  # Get hemibrain neurons
  if(missing(db)) {
    # this means we weren't told to use a specific neuronlist, so
    # we'll use the default. force() means evaluate hemibrain_neurons() now.
    db=tryCatch(force(db), error=function(e) {
        message("Unable to use `hemibrain_neurons()`. ",
                "I will read neurons from neuPrint, but this will be slower!")
      })
  } else if(is.character(db)) {
    db=tryCatch(get(db), error=function(e) stop("Unable to find neuronlist: ", db))
  }
  if(repository=="LM" & is.null(query)){
    fafb = hemibrain_lm_lhns(brainspace = c("JRCFIB2018F"))
  }else if(repository=="LM"){
    fafb = query
  }
  # How much is done?
  match.field = if(repository=="CATMAID"){
    "FAFB.match"
  }else if(repository == "flywire"){
    "flywire.xyz"
  }else{
    paste0(repository,".match")
  }
  quality.field = if(repository%in%c("CATMAID","flywire")){
    "FAFB.match.quality"
  }else{
    paste0(repository,".match.quality")
  }
  done = subset(gs, !is.na(gs[[match.field]]))
  message("Neuron matches: ", nrow(done))
  print(table(gs[[quality.field]]))
  # Choose user
  message("Users: ", paste(sort(unique(gs$User)),collapse = " "))
  initials = must_be("Choose a user : ", answers = sort(unique(gs$User)))
  say_hello(initials)
  rgl::bg3d("white")
  # choose ids
  selected = id_selector(gs=gs, ids=ids, id=id, overwrite = overwrite,
                         quality.field = quality.field, match.field = match.field,
                         initials = initials, column = column, field = field)
  ids = unique(selected[[id]])
  # Read neuron meta data
  meta = neuprintr::neuprint_get_meta(ids)
  meta = meta[order(meta$type),]
  meta = subset(meta, meta$bodyid%in%ids)
  # choose brain
  if(repository%in%c("flywire","CATMAID")){
    brain = elmr::FAFB14.surf
  }else{
    brain = hemibrainr::hemibrain_microns.surf
  }
  # Make matches!
  for(n in meta$bodyid){
    # Get bodyid
    n = as.character(n)
    # Remove neurons with matches
    if(n%in%donotdo$bodyid | !n%in%unlist(dimnames(hemibrain.nblast))){
      next
    }
    # Read hemibrain neuron
    if(is.null(db)){
      lhn  = neuprintr::neuprint_read_neurons(n, all_segments = TRUE, heal = FALSE)
      lhn = scale_neurons.neuronlist(lhn, scaling = (8/1000))
      lhn = suppressWarnings(nat.templatebrains::xform_brain(lhn, reference = "FAFB14", sample = "JRCFIB2018F"))
    } else {
      lhn = tryCatch(db[as.character(n)], error = function(e) {
        warning("Cannot read neuron: ", n, " from local db; fetching from neuPrint!")
        nat::as.neuronlist(neuprintr::neuprint_read_neurons(n, all_segments = TRUE, heal = FALSE))
        })
      lhn = scale_neurons.neuronlist(lhn, scaling = (8/1000))
      lhn = suppressWarnings(java_xform_brain(lhn, reference = "FAFB14", sample = "JRCFIB2018F"))
    }
    # Transform hemibrain neuron to FAFB space
    message("Hemibrain body ID: ", lhn[n,"bodyid"])
    message("Hemibrain-assigned cell type : ",lhn[n,"type"])
    # Read top 10 FAFB matches
    if(repository=="CATMAID"){
      message(sprintf("Reading the top %s FAFB matches from CATMAID",batch_size))
      r = tryCatch(sort(hemibrain.nblast[,as.character(n)],decreasing = TRUE), error = function(e) NULL)
      if(is.null(r)){
        r = tryCatch(sort(hemibrain.nblast[as.character(n),],decreasing = TRUE), error = function(e) NULL)
      }
      if(is.null(r)){
        message(n, " not in NBLAST matrix, skipping ...")
        next
      }
      fafb = catmaid::read.neurons.catmaid(names(r)[1:batch_size], OmitFailures = TRUE)
      j = batch_size
    }else if(repository=="flywire"){
      message(sprintf("Reading the top %s FAFB matches from flywire",batch_size))
      r = tryCatch(sort(hemibrain.nblast[,as.character(n)],decreasing = TRUE), error = function(e) NULL)
      if(is.null(r)){
        r = tryCatch(sort(hemibrain.nblast[as.character(n),],decreasing = TRUE), error = function(e) NULL)
      }
      if(is.null(r)){
        message(n, " not in NBLAST matrix, skipping ...")
        next
      }
      fafb = tryCatch(flywire_neurons(), error = function(e) NULL)
      if(is.null(fafb)){
        message("Flywire skeletons not found on drive, using fafbseg::skeletor ...")
        fafbseg::choose_segmentation("flywire")
        fafb = fafbseg::skeletor(names(r)[1:batch_size])
        fafb = flywire_basics(fafb)
        flywire.good = FALSE
      }else{
        flywire.good = TRUE
      }
      j = batch_size
    }else{
      r = sort(hemibrain.nblast[,n],decreasing = TRUE)
      j = length(fafb)
    }
    sel = c("go","for","it")
    k = 1
    # Plot brain
    rgl::clear3d()
    plot3d(brain, alpha = 0.1, col ="grey")
    plot3d(lhn[n], lwd = 2, soma = TRUE, col = "black")
    # Cycle through potential matches
    while(length(sel)>1){
      sel = sel.orig = nat::nlscan(fafb[names(r)[1:j]], col = "red", lwd = 2, soma = TRUE)
      if(length(sel)>1){
        message("Note: You selected more than one neuron")
      }
      prog = hemibrain_choice(sprintf("You selected %s neurons. Are you happy with that? ",length(sel)))
      while(length(sel)>1){
        message("Choose single best match: ")
        sel = nat::nlscan(fafb[as.character(sel.orig)], col = "orange", lwd = 2, soma = TRUE)
        message(sprintf("You selected %s neurons", length(sel)))
        if(!length(sel)){
          noselection = hemibrain_choice("You selected no neurons. Are you happy with that? ")
          if(!noselection){
            sel = sel.orig
          }
        }
      }
      if(!prog){
        sel = c("go","for","it")
        prog = hemibrain_choice(sprintf("Do you want to read %s more FAFB neurons? ", batch_size))
        if(prog){
          k = j
          j = j + batch_size
          if(repository=="CATMAID"){
          fafb = nat::union(fafb, catmaid::read.neurons.catmaid(names(r)[(k+1):j]), OmitFailures = TRUE)
          }else if(!flywire.good){
            fafb = nat::union(fafb, fafbseg::skeletor(names(r)[(k+1):j], OmitFailures = TRUE))
            fafb = flywire_basics(fafb)
          }
        }
      }
    }
    if(length(sel) & repository == "flywire"){
      sel = as.character(sel)
      if(!is.na(fafb[sel,"flywire.xyz"])){
        hit = fafb[sel,"flywire.xyz"]
      }else{
        fixed = flywire_basics(fafb[sel])
        hit = fixed[,"flywire.xyz"]
      }
    }else{
      hit = "none"
    }
    gs[n,match.field] = hit
    # Assign match and its quality
    if(length(sel)){
      rgl::plot3d(fafb[sel],col="blue",lwd=2,soma=TRUE)
      quality = must_be("What is the quality of this match? good(e)/okay(o)/poor(p)/tract-only(t) ", answers = c("e","o","p","t"))
    }else{
      quality = "n"
    }
    quality = standardise_quality(quality)
    gs[n,quality.field] = quality
    gs = gs[!duplicated(gs$bodyid),]
    unsaved = c(unsaved, n)
    message(length(unsaved), " unsaved matches")
    print(knitr::kable(gs[unsaved,c("bodyid","cell.type",match.field,quality.field)]))
    p = must_be("Continue (enter) or save (s)? ", answers = c("","s"))
    if(p=="s"){
      plot_inspirobot()
      say_encouragement(initials)
      # Read!
      gs2 = hemibrain_match_sheet(selected_file = selected_file, sheet = "hemibrain")
      gs2[match(gs[[id]],gs2[[id]]),match.field]=gs[[match.field]]
      gs2[match(gs[[id]],gs2[[id]]),quality.field]=gs[[quality.field]]
      # Write!
      write_matches(gs=gs,
                 ids = unsaved,
                 column = match.field)
      write_matches(gs=gs,
                 ids = unsaved,
                 column = quality.field)
      unsaved = c()
      gs = gs2
      rgl::bg3d("white")
    }
  }
  if(length(unsaved)){
    plot_inspirobot()
    say_encouragement(initials)
    # Read!
    gs2 = hemibrain_match_sheet(selected_file = selected_file, sheet = "hemibrain")
    gs2[match(gs[[id]],gs2[[id]]),match.field]=gs[[match.field]]
    gs2[match(gs[[id]],gs2[[id]]),quality.field]=gs[[quality.field]]
    # Write!
    write_matches(gs=gs,
                  ids = unsaved,
                  column = match.field)
    write_matches(gs=gs,
                  ids = unsaved,
                  column = quality.field)
  }
}

# hidden
## Google sheet write function
write_matches <- function(gs,
                       ids = NULL,
                       id.field = "bodyid",
                       ws = "hemibrain",
                       selected_file = options()$hemibrainr_matching_gsheet,
                       column = colnames(gs)){
  if((!identical(colnames(gs),column) & length(column) > 1)|(sum(column%in%colnames(gs))<1)){
    stop("Column must be one column of the Google sheet, or all the columns")
  }
  rows = (1:nrow(gs))+1
  rownames(gs) = rows
  if(!is.null(ids)){
    gs = subset(gs, gs[[id.field]] %in% ids)
    message("Updating ", nrow(gs), " entries")
    rows = rownames(gs)
  }
  for(r in rows){
    if(length(column)==1){
      letter = LETTERS[match(column,colnames(gs))]
      range = paste0(letter,r)
    }else{
      range = paste0("A",r,":",LETTERS[ncol(gs)],r)
    }
    gsheet_manipulation(FUN = googlesheets4::range_write,
                                     ss = selected_file,
                                     range = range,
                                     data = as.data.frame(gs[as.character(r),column], stringsAsFactors = FALSE),
                                     sheet = ws,
                                     col_names = FALSE)
  }
}

#' @rdname hemibrain_matching
#' @export
lm_matching <- function(ids = NULL,
                        hemibrain.nblast = NULL,
                        selected_file = options()$hemibrainr_matching_gsheet,
                        batch_size = 50,
                        db=hemibrain_neurons(),
                        query = NULL,
                        overwrite = c("FALSE","mine","mine_empty","TRUE")){
  # Motivate!
  nat::nopen3d()
  plot_inspirobot()
  message("
          #######################Colours##########################
          black = LM neuron,
          red = potential hemibrain matches based on NBLAST score,
          green = a chosen hemibrain neuron during scanning,
          blue = your selected hemibrain match.
          #######################Colours##########################
          ")
  ## Get NBLAST
  if(is.null(hemibrain.nblast)){
      message("Loading LM-FIB NBLAST from hemibrain Google Team Drive using Google Filestream: ")
    hemibrain.lhns.mean.compressed = hemibrain.dolan.mean.compressed = hemibrain.dolan.mean.compressed = hemibrain.lhins.mean.compressed = NULL
      load(sprintf("/Volumes/GoogleDrive/Shared\ drives/hemibrain/hemibrain_nblast/%s.rda", "hemibrain.lhns.mean.compressed"))
      load(sprintf("/Volumes/GoogleDrive/Shared\ drives/hemibrain/hemibrain_nblast/%s.rda", "hemibrain.dolan.mean.compressed"))
      load(sprintf("/Volumes/GoogleDrive/Shared\ drives/hemibrain/hemibrain_nblast/%s.rda", "hemibrain.lhins.mean.compressed"))
      m <- merge(hemibrain.dolan.mean.compressed, hemibrain.lhns.mean.compressed, by = 0, all = TRUE)
      m <- m[match(rownames(hemibrain.lhns.mean.compressed), m[, "Row.names"]), -1]
      rownames(m) <- rownames(hemibrain.lhns.mean.compressed)
      m2 <- merge(hemibrain.lhins.mean.compressed, m, by = 0, all = TRUE)
      rnams <- m2[, "Row.names"]
      hemibrain.nblast <- m2[, -1]
      rownames(hemibrain.nblast) <- rnams
      rm("hemibrain.lhns.mean.compressed")
      rm("hemibrain.dolan.mean.compressed")
      rm("hemibrain.lhins.mean.compressed")
  }
  if(is.null(query)){
    q1 = hemibrain_lm_lhns(brainspace = c("JRCFIB2018F"))
    q2 = hemibrain_lm_lhns(cable = "lhins", brainspace = c("JRCFIB2018F"))
    q3 = hemibrain_lm_lhns(cable = "lines", brainspace = c("JRCFIB2018F"))
    query = nat::union(q1, q2)
    query = nat::union(query, q3)
  }
  # Read the Google Sheet
  gs = gsheet_manipulation(FUN = googlesheets4::read_sheet,
                           ss = selected_file,
                           sheet = "lm",
                           guess_max = 3000,
                           return = TRUE)
  gs$id = correct_id(gs$id)
  rownames(gs) = gs$id
  gs[gs==""] = NA
  gs$User[is.na(gs$User)] = ""
  # Get hemibrain neurons
  if(missing(db)) {
    # this means we weren't told to use a specific neuronlist, so
    # we'll use the default. force() means evaluate hemibrain_neurons() now.
    db=tryCatch(force(db), error=function(e) {
      stop("Unable to use `hemibrain_neurons()`. ",
              "You must load the hemibrain Google Team Drive")
    })
  }
  # How much is done?
  match.field = paste0("hemibrain",".match")
  quality.field = paste0("hemibrain",".match.quality")
  done = subset(gs, !is.na(gs[[match.field]]))
  message("Neuron matches: ", nrow(done), "/", nrow(gs))
  print(table(gs[[quality.field]]))
  # Choose user
  message("Users: ", paste(sort(unique(gs$User)),collapse = " "))
  initials = must_be("Choose a user : ", answers = sort(unique(gs$User)))
  say_hello(initials)
  rgl::bg3d("white")
  # choose ids
  selected = id_selector(gs=gs, ids=ids, id=id, overwrite = overwrite,
                         quality.field = quality.field, match.field = match.field,
                         initials = initials, column = column, field = field)
  ids = unique(selected[[id]])
  # choose brain
  brain = hemibrainr::hemibrain.surf
  unsaved = c()
  # Make matches!
  for(n in selected$id){
    # Get id
    n = as.character(n)
    end = n==selected$id[length(selected$id)]
    # Remove neurons with matches
    if(!n%in%colnames(hemibrain.nblast)){
      message(n, " not found in NBLAST matrix, skipping")
      next
    }
    # Plot brain
    rgl::clear3d()
    plot3d(brain, alpha = 0.1, col ="grey")
    # Transform hemibrain neuron to FAFB space
    lhn = query[n]
    lhn = tryCatch( scale_neurons(lhn, scaling = (1000/8)),
                    error = function(e) lhn*(1000/8))
    if(nat::is.dotprops(lhn[[1]])){
      points3d(nat::xyzmatrix(lhn), col = "black")
    }else{
      plot3d(lhn, lwd = 2, soma = TRUE, col = "black")
    }
    message("ID: ", n)
    message("cell type : ",lhn[n,"cell.type"])
    # Read top 10 FAFB matches
    message(sprintf("Reading the top %s hemibrain hits",batch_size))
    r = rownames(hemibrain.nblast)[order(hemibrain.nblast[,n],decreasing = TRUE)]
    r = r[!is.na(r)]
    # Read hemibrain neurons
    if(is.null(db)){
      hemi = neuprintr::neuprint_read_skeletons((r[1:batch_size]),all_segments = TRUE, heal = FALSE)
    } else {
      batch = r[1:batch_size]
      batch.in = intersect(batch, names(db))
      hemi = tryCatch(db[match(batch.in,names(db))], error = function(e) NULL)
      if(is.null(hemi)|length(batch.in)!=length(batch)){
        message("Cannot read neurons from local db; fetching from neuPrint!")
        batch.out = setdiff(batch, names(hemi))
        hemi = c(hemi,tryCatch(neuprintr::neuprint_read_skeletons(batch.out, all_segments = TRUE, heal = FALSE),error=function(e) NULL))
        hemi = hemi[as.character(batch)]
      }
    }
    if(is.null(hemi)){
      next
    }
    sel = c("go","for","it")
    k = 1
    j = batch_size
    # Cycle through potential matches
    while(length(sel)>1){
      cycle = hemi[r[1:j]]
      cycle = cycle[!is.na(names(cycle))]
      sel = sel.orig = nat::nlscan(cycle, col = "red", lwd = 2, soma = TRUE)
      if(length(sel)>1){
        message("Note: You selected more than one neuron")
      }
      if(length(sel) > 0){
        rgl::plot3d(hemi[sel], lwd = 2, soma = TRUE)
      }
      prog = hemibrain_choice(sprintf("You selected %s neurons. Are you happy with that? ",length(sel)))
      if(length(sel)>0){
        nat::npop3d()
      }
      if(!prog){
        sel = c("go","for","it")
          prog = hemibrain_choice(sprintf("Do you want to read %s more neurons? ", batch_size))
          if(prog){
            k = j
            j = j + batch_size
            if(is.null(db)){
              hemi2  = neuprintr::neuprint_read_skeletons((r[(k+1):j]), all_segments = TRUE, heal = FALSE)
            } else {
              hemi2 = tryCatch(db[(names(r)[(k+1):j])], error = function(e) {
                warning("Cannot read neuron: ", n, " from local db; fetching from neuPrint!")
                neuprintr::neuprint_read_skeletons((r[(k+1):j]), all_segments = TRUE, heal = FALSE)
              })
            }
            hemi = nat::union(hemi, hemi2)
          }
      }else{
        while(length(sel)>1){
          message("Choose single best match: ")
          sel = nat::nlscan(hemi[as.character(sel.orig)], col = "orange", lwd = 2, soma = TRUE)
          message(sprintf("You selected %s neurons", length(sel)))
          if(!length(sel)){
            noselection = hemibrain_choice("You selected no neurons. Are you happy with that? ")
            if(!noselection){
              sel = sel.orig
            }
          }
        }
      }
    }
    # Assign match and its quality
    gs[n,match.field] = ifelse(length(sel)==0,'none',sel)
    if(length(sel)){
      rgl::plot3d(hemi[sel],col="blue",lwd=2,soma=TRUE)
      quality = must_be("What is the quality of this match? good(e)/okay(o)/poor(p)/tract-only(t) ", answers = c("e","o","p","t"))
    }else{
      quality = "n"
    }
    quality = standardise_quality(quality)
    gs[n,quality.field] = quality
    gs = gs[!duplicated(gs$id),]
    unsaved = c(unsaved, n)
    message(length(unsaved), " unsaved matches")
    print(knitr::kable(gs[unsaved,c("id","cell.type",match.field,quality.field)]))
    p = must_be("Continue (enter) or save (s)? ", answers = c("","s"))
    if(p=="s"|end){
      plot_inspirobot()
      say_encouragement(initials)
      # Read!
      gs2 = gsheet_manipulation(FUN = googlesheets4::read_sheet,
                                ss = selected_file,
                                sheet = "lm",
                                guess_max = 3000,
                                return = TRUE)
      gs2$id = correct_id(gs2$id)
      rownames(gs2) = gs2$id
      gs = gs[rownames(gs2),]
      # Write!
      write_matches(gs=gs,
                    ids = unsaved,
                    id.field = "id",
                    column = match.field,
                    ws = "lm")
      write_matches(gs=gs,
                    ids = unsaved,
                    id.field = "id",
                    column = quality.field,
                    ws = "lm")
      unsaved = c()
      gs = gs2
      rgl::bg3d("white")
    }
  }
  if(length(unsaved)){
    plot_inspirobot()
    say_encouragement(initials)
    # Read!
    gs2 = gsheet_manipulation(FUN = googlesheets4::read_sheet,
                              ss = selected_file,
                              sheet = "lm",
                              guess_max = 3000,
                              return = TRUE)
    gs2$id = correct_id(gs2$id)
    rownames(gs2) = gs2$id
    gs = gs[rownames(gs2),]
    # Write!
    write_matches(gs=gs,
                  ids = unsaved,
                  id.field = "id",
                  column = match.field,
                  ws = "lm")
    write_matches(gs=gs,
                  ids = unsaved,
                  id.field = "id",
                  column = quality.field,
                  ws = "lm")
  }
  say_encouragement(initials)
}


#' @rdname hemibrain_matching
#' @export
fafb_matching <- function(ids = NULL,
                        repository = c("CATMAID", "flywire"),
                        hemibrain.nblast = NULL,
                        selected_file = options()$hemibrainr_matching_gsheet,
                        batch_size = 50,
                        db=hemibrain_neurons(),
                        query = NULL,
                        overwrite = c("FALSE","mine","mine_empty","TRUE"),
                        column = NULL,
                        field = NULL){
  repository = match.arg(repository)
  # Packages
  if(!requireNamespace("nat.jrcbrains", quietly = TRUE)) {
    stop("Please install nat.jrcbrains using:\n", call. = FALSE,
         "remotes::install_github('natverse/nat.jrcbrains')")
  }
  if(!requireNamespace("elmr", quietly = TRUE)) {
    stop("Please install elmr using:\n", call. = FALSE,
         "remotes::install_github('natverse/elmr')")
  }
  # Motivate!
  nat::nopen3d()
  plot_inspirobot()
  unsaved = c()
  message("
          #######################Colours##########################
          black = FAFB CATMAID neuron,
          dark grey = flywire neuron,
          light grey = mirrrored flywire neuron,
          red = potential hemibrain matches based on NBLAST score,
          green = a chosen hemibrain neuron during scanning,
          blue = your selected hemibrain match.
          #######################Colours##########################
          ")
  ## Get NBLAST
  if(is.null(hemibrain.nblast) & repository == "CATMAID"){
    hemibrain.nblast = hemibrain_nblast("hemibrain-fafb14")
  }
  if(is.null(hemibrain.nblast) & repository == "flywire"){
    hemibrain.nblast = t(hemibrain_nblast("hemibrain-flywire"))
  }
  # FlyWire or CATMAID?
  if(repository=="CATMAID"){
    id = "skid"
  }else{
    id = "flywire.id"
    if(is.null(ids)){
      ids = names(flywire_neurons())
    }
  }
  # Read the Google Sheet
  gs = hemibrain_match_sheet(selected_file = selected_file, sheet = repository)
  gs$User[is.na(gs$User)] = ""
  # Get hemibrain neurons
  if(missing(db)) {
    # this means we weren't told to use a specific neuronlist, so
    # we'll use the default. force() means evaluate hemibrain_neurons() now.
    db=tryCatch(force(db), error=function(e) {
      stop("Unable to use `hemibrain_neurons()`. ",
           "You must load the hemibrain Google Team Drive")
    })
  }
  # How much is done?
  match.field = paste0("hemibrain",".match")
  quality.field = paste0("hemibrain",".match.quality")
  done = subset(gs, !is.na(gs[[match.field]]))
  message("Neuron matches: ", nrow(done), "/", nrow(gs))
  print(table(gs[[quality.field]]))
  # Choose user
  message("Users: ", paste(sort(unique(gs$User)),collapse = " "))
  initials = must_be("Choose a user : ", answers = unique(gs$User))
  say_hello(initials)
  rgl::bg3d("white")
  # choose ids
  selected = id_selector(gs=gs, ids=ids, id=id, overwrite = overwrite,
                         quality.field = quality.field, match.field = match.field,
                         initials = initials, column = column, field = field)
  ids = unique(selected[[id]])
  # choose brain
  brain = hemibrainr::hemibrain.surf
  # Load flywire neurons, if needed
  if(repository == 'flywire'){
    fw.neurons = flywire_neurons()
    fw.neurons.m = flywire_neurons(mirror=TRUE)
  }
  # Subset gs
  message("We'll look at ", nrow(selected)," ", repository, " neurons sequentially.")
  # Make matches!
  for(n in selected[[id]]){
    # Get id
    n = as.character(n)
    end = n==selected[[id]][length(selected[[id]])]
    # Plot brain
    rgl::clear3d()
    hemibrain_view()
    plot3d(brain, alpha = 0.1, col ="grey")
    # Transform hemibrain neuron to FAFB space
    if(!is.null(query)){
      lhn = query[n]
    }else if(repository == "CATMAID"){
      lhn = catmaid::read.neurons.catmaid(n, OmitFailures = TRUE)
    }else if(repository == "flywire"){
      sk = gs[n,]$skid[1]
      if(!is.na(sk)){
        lhn = tryCatch(catmaid::read.neurons.catmaid(sk, OmitFailures = TRUE), error = function(e) NULL)
      }else{
        lhn = NULL
      }
      fw.n = tryCatch({
        fw.n = fw.neurons[as.character(n)]
        fw.n = suppressWarnings(nat.templatebrains::xform_brain(fw.n, sample = "FAFB14", reference = "JRCFIB2018F"))
        scale_neurons(fw.n, scaling = (1000/8))
      }, error = function(e) {NULL})
      if(is.null(fw.n)){
        fw.n = tryCatch({
          message("Neuron not found on Google drive, attempting to read from flywire ...")
          if(!requireNamespace("fafbseg", quietly = TRUE)) {
            stop("Please install fafbseg using:\n", call. = FALSE,
                 "remotes::install_github('natverse/fafbseg')")
          }
          fw.n = fafbseg::skeletor(n)
          fw.n = suppressWarnings(nat.templatebrains::xform_brain(fw.n, sample = "FAFB14", reference = "JRCFIB2018F"))
          scale_neurons(fw.n, scaling = (1000/8))
        }, error = function(e) {NULL})
      }
      fw.m = tryCatch({
        fw.m = fw.neurons.m[names(fw.neurons.m) %in% as.character(n)]
        fw.m = suppressWarnings(nat.templatebrains::xform_brain(fw.m, sample = "FAFB14", reference = "JRCFIB2018F"))
        scale_neurons(fw.m, scaling = (1000/8))
      }, error = function(e) {NULL})
    }
    if(!is.null(lhn)){
      lhn = suppressWarnings(nat.templatebrains::xform_brain(lhn, sample = "FAFB14", reference = "JRCFIB2018F"))
      lhn = scale_neurons(lhn, scaling = (1000/8))
      plot3d(lhn, lwd = 2, soma = TRUE, col = "black")
    }
    if(repository=="flywire"){
      if(!is.null(fw.n)&length(is.null(fw.n))) {plot3d(fw.n, lwd = 3, soma = TRUE, col = "grey30")}
      if(!is.null(fw.m)&length(fw.m)) {plot3d(fw.m, lwd = 3, soma = TRUE, col = "grey50")}
    }
    message("ID: ", n)
    message("ItoLee_Hemilineage : ",lhn[n,"ItoLee_Hemilineage"])
    # Read top 10 FIB matches
    r = tryCatch(sort(hemibrain.nblast[,n],decreasing = TRUE), error = function(e) NULL)
    if(is.null(r)){
      message(n, " not in NBLAST matrix, skipping ...")
      next
    }
    message(sprintf("Reading the top %s hemibrain hits",batch_size))
    # Read hemibrain neurons
    if(is.null(db)){
      hemi  = neuprintr::neuprint_read_skeletons((names(r)[1:batch_size]), all_segments = TRUE, heal = FALSE)
    } else {
      batch = names(r)[1:batch_size]
      batch.in = intersect(batch, names(db))
      hemi = tryCatch(db[match(batch.in,names(db))], error = function(e) NULL)
      if(is.null(hemi)|length(batch.in)!=length(batch)){
        message("Cannot read neurons from local db; fetching from neuPrint!")
        batch.out = setdiff(batch, names(hemi))
        hemi2 =(tryCatch( nat::as.neuronlist(neuprintr::neuprint_read_skeletons(batch.out, all_segments = TRUE, heal = FALSE)), error=function(e)NULL) )
        if(!is.null(hemi2)){
          hemi = nat::union(hemi, hemi2)
          hemi = hemi[as.character(batch)]
        }
      }
    }
    sel = c("go","for","it")
    k = 1
    j = batch_size
    # Cycle through potential matches
    while(length(sel)>1){
      sel = sel.orig = tryCatch(nat::nlscan(hemi[names(r)[1:j]], col = "red", lwd = 2, soma = TRUE), error = function(e) NULL)
      if(is.null(sel)){
        next
      }
      if(length(sel)>1){
        message("Note: You selected more than one neuron")
      }
      if(length(sel) > 0){
        rgl::plot3d(hemi[sel], lwd = 2, soma = TRUE)
      }
      prog = hemibrain_choice(sprintf("You selected %s neurons. Are you happy with that? ",length(sel)))
      if(length(sel)>0){
        nat::npop3d()
      }
      if(!prog){
        sel = c("go","for","it")
        prog = hemibrain_choice(sprintf("Do you want to read %s more neurons? ", batch_size))
        if(prog){
          k = j
          j = j + batch_size
          if(is.null(db)){
            hemi2  = neuprintr::neuprint_read_skeletons((names(r)[(k+1):j]), all_segments = TRUE, heal = FALSE)
          } else {
            hemi2 = tryCatch(db[(names(r)[(k+1):j])], error = function(e) {
              warning("Cannot read neuron: ", n, " from local db; fetching from neuPrint!")
              neuprintr::neuprint_read_skeletons((names(r)[(k+1):j]), all_segments = TRUE, heal = FALSE)
            })
          }
          hemi = nat::union(hemi, hemi2)
        }
      }else{
        while(length(sel)>1){
          message("Choose single best match: ")
          sel = nat::nlscan(hemi[as.character(sel.orig)], col = "orange", lwd = 2, soma = TRUE)
          message(sprintf("You selected %s neurons", length(sel)))
          if(!length(sel)){
            noselection = hemibrain_choice("You selected no neurons. Are you happy with that? ")
            if(!noselection){
              sel = sel.orig
            }
          }
        }
      }
    }
    # Assign match and its quality
    gs[n,match.field] = ifelse(length(sel)==0,'none',sel)
    if(length(sel)){
      rgl::plot3d(hemi[sel],col="blue",lwd=2,soma=TRUE)
      quality = must_be("What is the quality of this match? good(e)/okay(o)/poor(p)/tract-only(t) ", answers = c("e","o","p","t"))
    }else{
      quality = "n"
    }
    quality = standardise_quality(quality)
    gs[n,quality.field] = quality
    unsaved = c(unsaved, n)
    message(length(unsaved), " unsaved matches")
    print(knitr::kable(gs[unsaved,c("skid","ItoLee_Hemilineage",match.field,quality.field)]))
    p = must_be("Continue (enter) or save (s)? ", answers = c("","s"))
    if(p=="s"|end){
      plot_inspirobot()
      say_encouragement(initials)
      # Read!
      gs2 = hemibrain_match_sheet(selected_file = selected_file, sheet = "flywire")
      gs2[match(gs[[id]],gs2[[id]]),match.field]=gs[[match.field]]
      gs2[match(gs[[id]],gs2[[id]]),quality.field]=gs[[quality.field]]
      # Write!
      write_matches(gs=gs2,
                    ids = unsaved,
                    id.field = id,
                    column = match.field,
                    ws = "FAFB")
      write_matches(gs=gs2,
                    ids = unsaved,
                    id.field = id,
                    column = quality.field,
                    ws = "FAFB")
      unsaved = c()
      gs = gs2
      rgl::bg3d("white")
    }
  }
  if(length(unsaved)){
    plot_inspirobot()
    say_encouragement(initials)
    # Read!
    gs2 = hemibrain_match_sheet(selected_file = selected_file, sheet = "flywire")
    gs2[match(gs[[id]],gs2[[id]]),match.field]=gs[[match.field]]
    gs2[match(gs[[id]],gs2[[id]]),quality.field]=gs[[quality.field]]
    # Write!
    write_matches(gs=gs2,
                  ids = unsaved,
                  id.field = id,
                  column = match.field,
                  ws = "FAFB")
    write_matches(gs=gs2,
                  ids = unsaved,
                  id.field = id,
                  column = quality.field,
                  ws = "FAFB")
  }
  say_encouragement(initials)
}


#' Retrieve matched up neurons between the hemibrain and FAFB
#'
#' @description Many neurons in the hemibrain data have been matched one to one with FAFB neurons (hemibrain->FAFB) and
#' FAFB neurons have been matched with their hemibrain counterparts (FAFB->hemibrain). These matches have been done by the
#' Flyconnectome Group at the University of Cambridge, and are recorded on a
#' \href{https://docs.google.com/spreadsheets/d/1OSlDtnR3B1LiB5cwI5x5Ql6LkZd8JOS5bBr-HTi0pOw/edit#gid=0}{Google
#'   Sheet} on our hemibrain Google Team Drive, to which you will need access
#' through an authenticated account to view and use this function. Matches have three levels of 'quality', largely
#' dependent on the degree of manual tracing for FAFB neurons - good (could be the same cell), medium (same cell type) and poor (could be the same or similar cell type).
#'
#' @param priority whether to use FAFB->hemibrain matches (FAFB) or hemibrain->FAFB matches (hemibrain) in order to ascribe
#' cell type names to FAFB neurons. In both cases, cell type names are attached to hemibrain bodyids, and propagated to their FAFB matches.
#'
#' @return  a \code{data.frame} which includes neuron's ID (either its CATMAID skeleton ID or neuprint body ID), the data set from which it comes,
#' its putative cell type and connectivity type, and its match in the other dataset.
#'
#' @details Currently, the
#'   \href{https://docs.google.com/spreadsheets/d/1OSlDtnR3B1LiB5cwI5x5Ql6LkZd8JOS5bBr-HTi0pOw/edit#gid=0}{Google
#'    Sheet} is set up with  limited number of users, each of whom have been
#'   assigned a number of neurons to match up. In order to add yourself as a
#'   user, simply open this Google Sheet in your browser and add your initials
#'   to neurons of your choosing on the rightmost column 'Users'.
#'
#' @examples
#' \donttest{
#' \dontrun{
#'
#' # Get matches
#' matched = hemibrain.matches()
#'
#' }}
#' @rdname hemibrain_matches
#' @export
#' @seealso \code{\link{hemibrain_matching}}
#' @importFrom stats ave
hemibrain_matches <- function(priority = c("FAFB","hemibrain")){
  priority = match.arg(priority)

  # Get matches
  selected_file = options()$hemibrainr_matching_gsheet
  hemibrain.matches = gsheet_manipulation(FUN = googlesheets4::read_sheet,
                                                       ss = selected_file,
                                                       sheet = "hemibrain",
                                                       return = TRUE)
  hemibrain.matches$bodyid = correct_id(hemibrain.matches$bodyid)
  hemibrain.matches = hemibrain.matches[!duplicated(hemibrain.matches$bodyid),]
  hemibrain.matches = hemibrain.matches[hemibrain.matches$bodyid!="",]
  hemibrain.matches$ItoLee_Lineage = gsub("_.*","",hemibrain.matches$ItoLee_Hemilineage)
  hemibrain.matches$dataset = "hemibrain"
  hemibrain.matches = subset(hemibrain.matches, !is.na(hemibrain.matches$bodyid))
  hemibrain.matches = hemibrain.matches[!duplicated(hemibrain.matches$bodyid),]
  rownames(hemibrain.matches) = hemibrain.matches$bodyid
  #hemibrain.matches$cell.type = hemibrain.matches$connectivity.type

  # Get FAFB matches
  fafb.matches = gsheet_manipulation(FUN = googlesheets4::read_sheet,
                                                  ss = selected_file,
                                                  sheet = "FAFB",
                                                  return = TRUE)
  fafb.matches$skid = correct_id(fafb.matches$skid)
  fafb.matches = subset(fafb.matches, !grepl("l",fafb.matches$side))
  fafb.matches = fafb.matches[!duplicated(fafb.matches$skid),]
  fafb.matches = fafb.matches[fafb.matches$skid!="",]
  fafb.matches = subset(fafb.matches, !is.na(fafb.matches$skid))
  fafb.matches = fafb.matches[!duplicated(fafb.matches$skid),]
  fafb.matches$dataset = "FAFB"
  rownames(fafb.matches) = fafb.matches$skid

  # Unmatched
  unmatched = is.na(fafb.matches$hemibrain.match)|fafb.matches$hemibrain.match==""
  fafb.matches$hemibrain.match.quality[unmatched] = "none"
  fafb.matches$hemibrain.match[unmatched] = "none"

  # Add in neuprint types
  meta = neuprintr::neuprint_get_meta(hemibrain.matches$bodyid)
  types = meta$type
  names(types) = meta$bodyid
  hemibrain.matches[,"cell.type"] = types[as.character(hemibrain.matches$bodyid)]

  # Address FAFB cell types
  fafb.matches$cell.type = NA
  for(q in c("good","medium","poor")){
    if(is.na(q)){
      next
    }
    if(priority=="FAFB"){
      # FAFB -> Hemibrain
      inq = as.character(subset(fafb.matches, fafb.matches$hemibrain.match.quality==q)$skid)
      inq = inq[!is.na(inq)]
      isna = is.na(fafb.matches[inq,]$cell.type)
      cts = hemibrain.matches$cell.type[match(fafb.matches[inq,]$hemibrain.match,hemibrain.matches$bodyid)]
      fafb.matches[inq,]$cell.type[isna] = cts[isna]
      # Hemibrain -> FAFB
      inq = as.character(subset(hemibrain.matches, hemibrain.matches$FAFB.match.quality==q)$FAFB.match)
      inq = inq[!is.na(inq)]
      isna = is.na(fafb.matches[inq,]$cell.type)
      fafb.matches[inq,]$cell.type[isna] = hemibrain.matches$cell.type[match(fafb.matches[inq,]$skid[isna],hemibrain.matches$FAFB.match)]
    }else{
      # Hemibrain -> FAFB
      inq = as.character(subset(hemibrain.matches, hemibrain.matches$FAFB.match.quality==q)$FAFB.match)
      inq = inq[!is.na(inq)]
      isna = is.na(fafb.matches[inq,]$cell.type)
      fafb.matches[inq,]$cell.type[isna] = hemibrain.matches$cell.type[match(fafb.matches[inq,]$skid[isna],hemibrain.matches$FAFB.match)]
      # FAFB -> Hemibrain
      inq = as.character(subset(fafb.matches, fafb.matches$hemibrain.match.quality==q)$skid)
      inq = inq[!is.na(inq)]
      isna = is.na(fafb.matches[inq,]$cell.type)
      cts = hemibrain.matches$cell.type[match(fafb.matches[inq,]$hemibrain.match,hemibrain.matches$bodyid)]
      fafb.matches[inq,]$cell.type[isna] = cts[isna]
    }
  }

  # Add in neuprint types
  ntotype = fafb.matches$hemibrain.match[is.na(fafb.matches$cell.type)]
  ntotype = ntotype[ntotype!="none" & !is.na(ntotype)]
  meta = neuprintr::neuprint_get_meta(ntotype)
  types = meta$type
  names(types) = meta$bodyid
  types = types[!is.na(types)]
  fafb.matches[match(names(types),fafb.matches$hemibrain.match),"cell.type"] = types

  # Work out lineages
  for(id in as.character(fafb.matches$skid)){
    if(is.na(id)){
      next
    }
    ct = fafb.matches[id,"cell.type"]
    if(is.na(ct)){
      fafb.matches[id,"cell.type"] = hemibrain.matches$cell.type[match(id,hemibrain.matches$FAFB.match)]
      if(is.na(fafb.matches[id,"cell.type"])){
        fafb.matches[id,"cell.type"] = "uncertain"
      }else{
        fafb.matches[id,"ItoLee_Hemilineage"] = hemibrain.matches$ItoLee_Hemilineage[match(id,hemibrain.matches$FAFB.match)]
      }
    }else{
      fafb.matches[id,"ItoLee_Hemilineage"] = hemibrain.matches$ItoLee_Hemilineage[match(ct,hemibrain.matches$cell.type)]
    }
  }

  # Add cell body fiber infor for FAFB cells
  fafb.matches$cellBodyFiber = hemibrain.matches$cellBodyFiber[match(fafb.matches$hemibrain.match,hemibrain.matches$bodyid)]

  # Rename cells
  fafb.matches$cell = paste0(fafb.matches$cell.type,"#",ave(fafb.matches$cell.type,fafb.matches$cell.type,FUN= seq.int))
  hemibrain.matches$cell = paste0(hemibrain.matches$cell.type,"#",ave(hemibrain.matches$cell.type,hemibrain.matches$cell.type,FUN= seq.int))

  # Make matching data frame
  hemibrain.matches$dataset = "hemibrain"
  fafb.matches$dataset = "FAFB"
  matched.h = hemibrain.matches[,c("bodyid", "cell.type", "cell", "cellBodyFiber", "ItoLee_Hemilineage",
                                   "FAFB.match", "FAFB.match.quality", "LM.match", "LM.match.quality", "dataset")]
  matched.f = fafb.matches[,c("skid",  "cell.type",  "cell", "cellBodyFiber", "ItoLee_Hemilineage",
                              "hemibrain.match", "hemibrain.match.quality", "LM.match", "LM.match.quality","dataset")]
  colnames(matched.h) = colnames(matched.f) = c("id","cell.type", "cell","cellBodyFiber","ItoLee_Hemilineage","match","quality", "LM.match", "LM.match.quality","dataset")
  matched = rbind(matched.h,matched.f)
  matched$quality[is.na(matched$match)] = "none"
  matched$match[is.na(matched$match)] = "none"
  matched$ItoLee_Lineage = gsub("_.*","",matched$ItoLee_Hemilineage)

  # Sort out types
  matched$connectivity.type = matched$cell.type
  matched$cell.type = gsub("_.*","",matched$cell.type)
  matched$cell.type[is.na(matched$cell.type)] = "uncertain"
  matched$connectivity.type[is.na(matched$connectivity.type)] = "uncertain"

  # Return
  matched
}

#' @rdname hemibrain_matches
#' @export
lm_matches <- function(priority = c("hemibrain","lm")){
  priority = match.arg(priority)

  # Get matches
  selected_file = options()$hemibrainr_matching_gsheet
  hemibrain.matches = gsheet_manipulation(FUN = googlesheets4::read_sheet,
                                          ss = selected_file,
                                          sheet = "hemibrain",
                                          return = TRUE)
  hemibrain.matches$bodyid = correct_id(hemibrain.matches$bodyid)
  hemibrain.matches = hemibrain.matches[!duplicated(hemibrain.matches$bodyid),]
  hemibrain.matches = hemibrain.matches[hemibrain.matches$bodyid!="",]
  hemibrain.matches$ItoLee_Lineage = gsub("_.*","",hemibrain.matches$ItoLee_Hemilineage)
  hemibrain.matches$dataset = "hemibrain"
  hemibrain.matches = subset(hemibrain.matches, !is.na(hemibrain.matches$bodyid))
  hemibrain.matches = hemibrain.matches[!duplicated(hemibrain.matches$bodyid),]
  rownames(hemibrain.matches) = hemibrain.matches$bodyid
  hemibrain.matches$cell.type = hemibrain.matches$connectivity.type

  # Get lm matches
  lm.matches = gsheet_manipulation(FUN = googlesheets4::read_sheet,
                                     ss = selected_file,
                                     sheet = "lm",
                                     return = TRUE)
  lm.matches$id = correct_id(lm.matches$id)
  lm.matches = lm.matches[!duplicated(lm.matches$id),]
  lm.matches = lm.matches[lm.matches$id!="",]
  lm.matches = subset(lm.matches, !is.na(lm.matches$id))
  lm.matches = lm.matches[!duplicated(lm.matches$id),]
  lm.matches$dataset = "lm"
  rownames(lm.matches) = lm.matches$id

  # Unmatched
  unmatched = is.na(lm.matches$hemibrain.match)|lm.matches$hemibrain.match==""
  lm.matches$hemibrain.match.quality[unmatched] = "none"
  lm.matches$hemibrain.match[unmatched] = "none"

  # Add in neuprint types
  meta = neuprintr::neuprint_get_meta(hemibrain.matches$bodyid)
  types = meta$type
  names(types) = meta$bodyid
  hemibrain.matches[,"cell.type"] = types[as.character(hemibrain.matches$bodyid)]

  # Address lm cell types
  lm.matches$cell.type = NA
  for(q in c("good","medium","poor")){
    if(is.na(q)){
      next
    }
    if(priority=="lm"){
      # lm -> Hemibrain
      inq = as.character(subset(lm.matches, lm.matches$hemibrain.match.quality==q)$id)
      inq = inq[!is.na(inq)]
      isna = is.na(lm.matches[inq,]$cell.type)
      cts = hemibrain.matches$cell.type[match(lm.matches[inq,]$hemibrain.match,hemibrain.matches$bodyid)]
      lm.matches[inq,]$cell.type[isna] = cts[isna]
      # Hemibrain -> lm
      inq = as.character(subset(hemibrain.matches, hemibrain.matches$lm.match.quality==q)$lm.match)
      inq = inq[!is.na(inq)]
      isna = is.na(lm.matches[inq,]$cell.type)
      lm.matches[inq,]$cell.type[isna] = hemibrain.matches$cell.type[match(lm.matches[inq,]$id[isna],hemibrain.matches$lm.match)]
    }else{
      # Hemibrain -> lm
      inq = as.character(subset(hemibrain.matches, hemibrain.matches$lm.match.quality==q)$lm.match)
      inq = inq[!is.na(inq)]
      isna = is.na(lm.matches[inq,]$cell.type)
      lm.matches[inq,]$cell.type[isna] = hemibrain.matches$cell.type[match(lm.matches[inq,]$id[isna],hemibrain.matches$lm.match)]
      # lm -> Hemibrain
      inq = as.character(subset(lm.matches, lm.matches$hemibrain.match.quality==q)$id)
      inq = inq[!is.na(inq)]
      isna = is.na(lm.matches[inq,]$cell.type)
      cts = hemibrain.matches$cell.type[match(lm.matches[inq,]$hemibrain.match,hemibrain.matches$bodyid)]
      lm.matches[inq,]$cell.type[isna] = cts[isna]
    }
  }

  # Add in neuprint types
  ntotype = lm.matches$hemibrain.match[is.na(lm.matches$cell.type)]
  ntotype = ntotype[ntotype!="none" & !is.na(ntotype)]
  meta = neuprintr::neuprint_get_meta(ntotype)
  types = meta$type
  names(types) = meta$bodyid
  types = types[!is.na(types)]
  lm.matches[match(names(types),lm.matches$hemibrain.match),"cell.type"] = types

  # Work out lineages
  for(id in as.character(lm.matches$id)){
    if(is.na(id)){
      next
    }
    ct = lm.matches[id,"cell.type"]
    if(is.na(ct)){
      lm.matches[id,"cell.type"] = hemibrain.matches$cell.type[match(id,hemibrain.matches$lm.match)]
      lm.matches[id,"ItoLee_Hemilineage"] = hemibrain.matches$ItoLee_Hemilineage[match(id,hemibrain.matches$lm.match)]
    }else{
      lm.matches[id,"ItoLee_Hemilineage"] = hemibrain.matches$ItoLee_Hemilineage[match(ct,hemibrain.matches$cell.type)]
    }
  }

  # Rename cells
  lm.matches$cell = paste0(lm.matches$cell.type,"#",ave(lm.matches$cell.type,lm.matches$cell.type,FUN= seq.int))
  hemibrain.matches$cell = paste0(hemibrain.matches$cell.type,"#",ave(hemibrain.matches$cell.type,hemibrain.matches$cell.type,FUN= seq.int))

  # Fix hemilineages
  hl = hemibrain.matches$ItoLee_Hemilineage[match(lm.matches$id,hemibrain.matches$lm.match)]
  l = hemibrain.matches$ItoLee_Lineage[match(lm.matches$id,hemibrain.matches$lm.match)]
  lm.matches$ItoLee_Hemilineage[!is.na(hl)] = hl[!is.na(hl)]
  lm.matches$ItoLee_Lineage[!is.na(l)] = l[!is.na(l)]

  # Add cell body fiber infor for FAFB cells
  lm.matches$cellBodyFiber = hemibrain.matches$cellBodyFiber[match(lm.matches$hemibrain.match,hemibrain.matches$bodyid)]

  # Make matching data frame
  hemibrain.matches$dataset = "hemibrain"
  lm.matches$dataset = "lm"
  matched.h = hemibrain.matches[,c("bodyid", "cell.type", "cell", "cellBodyFiber", "ItoLee_Hemilineage",
                                   "LM.match", "LM.match.quality", "FAFB.match", "FAFB.match.quality", "dataset")]
  matched.f = lm.matches[,c("id",  "cell.type",  "cell", "cellBodyFiber", "ItoLee_Hemilineage",
                              "hemibrain.match", "hemibrain.match.quality", "FAFB.match", "FAFB.match.quality","dataset")]
  colnames(matched.h) = colnames(matched.f) = c("id","cell.type", "cell","cellBodyFiber","ItoLee_Hemilineage","match","quality", "FAFB.match", "FAFB.match.quality","dataset")
  matched = rbind(matched.h,matched.f)
  matched$quality[is.na(matched$match)] = "none"
  matched$match[is.na(matched$match)] = "none"

  # Sort out types
  matched$connectivity.type = matched$cell.type
  matched$cell.type = gsub("_.*","",matched$cell.type)
  matched$cell.type[is.na(matched$cell.type)] = "uncertain"
  matched$connectivity.type[is.na(matched$connectivity.type)] = "uncertain"

  # Return
  matched

}


#' Manage hemibrain-FAFB neuron matches
#'
#' @description We can match neurons in the hemibrain data with FAFB neurons (hemibrain->FAFB, hemibrain tab) and
#' FAFB neurons have been matched with their hemibrain counterparts (FAFB->hemibrain, fafb tab). These matches are recorded on a
#' \href{https://docs.google.com/spreadsheets/d/1OSlDtnR3B1LiB5cwI5x5Ql6LkZd8JOS5bBr-HTi0pOw/edit#gid=0}{Google
#'   Sheet} on our hemibrain Google Team Drive, to which you will need access
#' through an authenticated account to view and use this function. We can use these function to manipulate the Google Sheet.
#' This include adding new IDs for matching (\code{hemibrain_matching_add}), transfer matches between the sheet's hemibrain and FAFB tabs (\code{hemibrain_matching_transfers})
#' and add pre-made matches (\code{hemibrain_add_matches}).The function \code{lywire_matching_rewrite} updates a Google sheet with correct flywire IDs and stable
#' locations in flywire space, based on CATMAID skids.
#'
#' @param df a \code{data.frame} of pre-made matches that the user wants to transfer onto the Google Sheet. This will erase any extant matches for the specified neurons. This data frame
#' must have columns: bodyid (i.e. the hemibrain neuron), skid (i.e. the FAFB neuron) and quality (i.e. match quality). Matches have three levels of 'quality', largely dependent on the degree of manual tracing for FAFB neurons
#'  - good (could be the same cell), medium (same cell type) and poor (could be the same or similar cell type).
#' @param ids either hemibrain bodyids or FAFB skids to add to a Google Sheet. You will want to add IDs where they do not already exist, so that
#'  you can use, for example, \code{\link{hemibrain_matching}}.
#' @param dataset the tab to which to add your new information. You are either adding to information related to hemibrain neurons, or FAFB neurons.
#' @param direction the match direction, i.e. hemibrain->FAFB (hemibrain tab) or FAFB->hemibrain (fafb tab). Defaults to updating both.
#' @param selected_file Specifies which Google Sheet to use. Unless you are using a personal Google Sheet, this should be \code{options()$hemibrainr_matching_gsheet}.
#' @param User your initials, so updates can be linked to you. There is a limited number of users, each of whom have been
#'   assigned a number of neurons to match up. In order to add yourself as a user,
#'   simply open this Google Sheet in your browser and add your initials to neurons of your choosing on the rightmost column 'Users'.
#' @param flywire.ids flywire IDs to add to Google sheet if not already present.
#'
#' @param ... arguments passed to methods for, for example, \code{neuprintr::neuprint_get_meta} and \code{elmr::fafb_get_meta}.
#'
#' @return  \code{NULL}. Updates the master Google sheet.
#' @details Currently, information is recorded in a
#'   \href{https://docs.google.com/spreadsheets/d/1OSlDtnR3B1LiB5cwI5x5Ql6LkZd8JOS5bBr-HTi0pOw/edit#gid=0}{Google
#'    Sheet}
#'
#' @examples
#' \donttest{
#' \dontrun{
#'
#' # Add a mising FAFB projection neuron, so we can match it later:
#' hemibrain_matching_add(ids = "16", sheet = "FAFB", User = "ASB")
#'
#'
#' }}
#' @seealso \code{\link{hemibrain_matching}}
#' @rdname hemibrain_add_made_matches
hemibrain_add_made_matches <- function(df,
                                  direction = c("both","hemibrain-FAFB","FAFB-hemibrain"),
                                  User = "flyconnectome",
                                  ...){
  direction = match.arg(direction)
  cnames = c("bodyid","skid","quality")
  if(!all(cnames%in%colnames(df))){
    stop("df must have column names: ", paste(cnames,collapse = ", "))
  }else{
    df = subset(df,
                !is.na(df$skid)
                & !is.na(df$bodyid)
                & df$skid!=""
                & df$bodyid != ""
                & df$quality %in% c("none","good","medium","poor"))
  }

  # Update hemibrain sheet
  if(direction%in%c("both","hemibrain-FAFB")){
    message("Reading hemibrain sheet")
    gs = hemibrain_match_sheet(sheet = "hemibrain")
    missing = setdiff(df$bodyid,gs$bodyid)
    if(length(missing)){
      message("Adding missing hemibrain bodyids")
      hemibrain_matching_add(ids = missing, sheet = "hemibrain", User = User)
      gs = hemibrain_match_sheet(sheet = "hemibrain")
    }
    message("Checking that FAFB matches exist")
    hdf = subset(df, df$bodyid %in% gs$bodyid)
    hdf = hdf[!duplicated(hdf$bodyid),]
    in.fab = suppressWarnings(sapply(hdf$skid,catmaid::catmaid_skids))
    in.fab = sapply(in.fab,function(x) length(x)!=0)
    if(sum(!in.fab)){
      message("Dropping ", sum(!in.fab)," matches as givens skeleton IDs do not exist for FAFB v14")
      hdf = hdf[in.fab,]
    }
    gs[as.character(hdf$bodyid),"FAFB.match"] = hdf$skid
    gs[as.character(hdf$bodyid),"FAFB.match.quality"] = hdf$quality
    gs[as.character(hdf$bodyid),"User"] = User
    message("Adding matches")
    write_matches(gs=gs,
                  ids = as.character(hdf$bodyid),
                  ws="hemibrain",
                  id.field ="bodyid",
                  column = "FAFB.match")
    write_matches(gs=gs,
                  ids = as.character(hdf$bodyid),
                  ws="hemibrain",
                  id.field ="bodyid",
                  column = "FAFB.match.quality")
  }
  if(direction%in%c("both","FAFB-hemibrain")){
    message("Reading FAFB sheet")
    gs = hemibrain_match_sheet(sheet = "FAFB")
    missing = setdiff(df$skid,gs$skid)
    if(length(missing)){
      message("Adding missing FAFB bodyids")
      hemibrain_matching_add(ids = missing, sheet = "FAFB", User = User)
      gs = hemibrain_match_sheet(sheet = "FAFB")
    }
    message("Checking that hemibrain matches exist")
    hdf = subset(df, df$skid %in% gs$skid)
    hdf = hdf[!duplicated(hdf$skid),]
    in.hemi = suppressWarnings(sapply(hdf$skid,function(x) tryCatch(neuprintr::neuprint_ids(x),error = function(e) NULL)))
    in.hemi = sapply(in.hemi,function(x) length(x)!=0)
    if(sum(!in.hemi)){
      message("Dropping ", sum(!in.hemi)," matches as givens skeleton IDs do not exist for FAFB v14")
      hdf = hdf[in.hemi,]
    }
    gs[as.character(hdf$skid),"hemibrain.match"] = hdf$bodyid
    gs[as.character(hdf$skid),"hemibrain.match.quality"] = hdf$quality
    gs[as.character(hdf$skid),"User"] = User
    message("Adding matches")
    write_matches(gs=gs,
                  ids = as.character(hdf$skid),
                  ws="FAFB",
                  id.field ="skid",
                  column = "hemibrain.match")
    write_matches(gs=gs,
                  ids = as.character(hdf$skid),
                  ws="FAFB",
                  id.field ="skid",
                  column = "hemibrain.match.quality")
  }
}

# Get correct GSheet
hemibrain_match_sheet <- function(selected_file = options()$hemibrainr_matching_gsheet,
                                  sheet = c("hemibrain","FAFB","CATMAID","flywire")){
  # Which sheet
  sheet = match.arg(sheet)
  sheet[sheet=="hemibrain"] = "hemibrain"
  # neuron ID name
  if(sheet=="hemibrain"){
    id.field = "bodyid"
  }else if (sheet == "flywire"){
    id.field = "flywire.id"
    sheet = "FAFB"
  }else if (sheet=="CATMAID"){
    sheet = "FAFB"
  }else{
    id.field = "skid"
  }
  # Read sheet
  gs = gsheet_manipulation(FUN = googlesheets4::read_sheet,
                                        ss = selected_file,
                                        sheet = sheet,
                                        return = TRUE)
  gs[[id.field]] = correct_id(gs[[id.field]])
  ids = gs[[id.field]]
  ids[is.na(ids)] = paste0("missing_",1:sum(is.na(ids)))
  ids = paste0(ids,"#",ave(ids,ids,FUN= seq.int))
  ids = gsub("#1$","",ids)
  rownames(gs) = ids
  gs
}


#' @rdname hemibrain_add_made_matches
#' @export
hemibrain_matching_add <- function(ids = NULL,
                                   dataset = c("hemibrain","FAFB","CATMAID","flywire"),
                                   User = dataset,
                                   selected_file  = options()$hemibrainr_matching_gsheet,
                                   ...){
  # Read
  dataset = match.arg(dataset)

  # Have we got the right arguments?
  if(is.null(ids)&dataset!='flywire'){
    stop("Please provide IDs with which you want to update ", dataset, " related information")
  }else if(dataset == "CATMAID"){
    dataset = "FAFB"
  }

  # get gs
  gs = hemibrain_match_sheet(sheet = dataset, selected_file = selected_file)

  # What are we adding?
  if(is.null(ids)&dataset == "flywire"){
    fw.ids = names(flywire_neurons())
    all.ids = unique(gs$flywire.ids)
    add = setdiff(fw.ids, all.ids)
  }
  add = setdiff(ids, rownames(gs))

  # Meta information
  if(dataset=="hemibrain"){
    meta = hemibrain_get_meta(add, ...)
    meta$cell.type = meta$type
  }else if (dataset == "FAFB"){
    meta = elmr::fafb_get_meta(add, ...)
  } else if (dataset == "flywire"){
    add = setdiff(ids, gs$flywire.id)
    meta = flywire_neurons()[as.character(add),]
    if(!nrow(meta)){
      stop("Selected IDs could not be added. They must be among the neurons
         saved on Google drive, see flywire_neurons()")
    }
  }
  if(!length(add)){
    stop("Given IDs already exist in sheet")
  }
  missing = setdiff(colnames(gs),colnames(meta))
  meta = add_blanks(meta, missing)
  meta = meta[,colnames(gs)]
  meta$User = User

  # Add new rows
  sheet[sheet=="hemibrain"] = "hemibrain"
  batches = split(1:nrow(meta), ceiling(seq_along(1:nrow(meta))/500))
  for(i in batches){
    gsheet_manipulation(FUN = googlesheets4::sheet_append,
                                     data = meta[min(i):max(i),],
                                     ss = selected_file,
                                     sheet = sheet)
  }

}

#' @rdname hemibrain_add_made_matches
#' @export
hemibrain_matching_transfers <- function(selected_file = options()$hemibrainr_matching_gsheet){
  #############
  # Transfers #
  #############

  # Read the LM Google Sheet
  lmg = gsheet_manipulation(FUN = googlesheets4::read_sheet,
                            ss = selected_file,
                            sheet = "lm",
                            guess_max = 3000,
                            return = TRUE)
  lmg$id = correct_id(lmg$id)
  rownames(lmg) = lmg$id

  # Read the FAFB Google Sheet
  fg = gsheet_manipulation(FUN = googlesheets4::read_sheet,
                           ss = selected_file,
                           sheet = "FAFB",
                           guess_max = 3000,
                           return = TRUE)
  fg$skid = correct_id(fg$skid)
  rownames(fg) = fg$skid

  # Read the hemibrain Google Sheet
  hg = gsheet_manipulation(FUN = googlesheets4::read_sheet,
                           ss = selected_file,
                           sheet = "hemibrain",
                           guess_max = 3000,
                           return = TRUE)
  hg$bodyid = correct_id(hg$bodyid)
  rownames(hg) = hg$bodyid

  ##############
  # LM -> FAFB #
  ##############
  missing = is.na(fg$LM.match)
  matches = lmg$id[match(fg$skid[missing],lmg$FAFB.match)]
  quality = lmg$FAFB.match.quality[match(fg$skid[missing],lmg$FAFB.match)]
  fg[missing,"LM.match"] = matches
  fg[missing,"LM.match.quality"] = quality
  sorted = which(missing)[!is.na(matches)]
  if(length(sorted)){
    update_gsheet(update = fg[sorted,],
                gs = fg,
                tab = "FAFB",
                match = "lm",
                id = "skid")
  }

  ###################
  # LM -> Hemibrain #
  ###################
  missing = is.na(hg$LM.match)
  matches = lmg$id[match(hg$bodyid[missing],lmg$hemibrain.match)]
  quality = lmg$hemibrain.match.quality[match(hg$bodyid[missing],lmg$hemibrain.match)]
  hg[missing,"LM.match"] = matches
  hg[missing,"LM.match.quality"] = quality
  sorted = which(missing)[!is.na(matches)]
  if(length(sorted)){
    update_gsheet(update = hg[sorted,],
                gs = hg,
                tab = "hemibrain",
                match = "lm",
                id = "bodyid")
  }

  ##############
  # FAFB -> LM #
  ##############
  missing = is.na(lmg$FAFB.match)
  matches = fg$skid[match(lmg$id[missing],fg$LM.match)]
  quality = fg$LM.match.quality[match(lmg$id[missing],fg$LM.match)]
  lmg[missing,"FAFB.match"] = matches
  lmg[missing,"FAFB.match.quality"] = quality
  sorted = which(missing)[!is.na(matches)]
  if(length(sorted)){
    update_gsheet(update = lmg[sorted,],
                  gs = lmg,
                  tab = "lm",
                  match = "FAFB",
                  id = "id")
  }

  #####################
  # FAFB -> Hemibrain #
  #####################
  missing = is.na(hg$FAFB.match)
  matches = fg$skid[match(hg$bodyid[missing],fg$hemibrain.match)]
  quality = fg$hemibrain.match.quality[match(hg$bodyid[missing],fg$hemibrain.match)]
  hg[missing,"FAFB.match"] = matches
  hg[missing,"FAFB.match.quality"] = quality
  sorted = which(missing)[!is.na(matches)]
  if(length(sorted)){
    update_gsheet(update = hg[sorted,],
                gs = hg,
                tab = "hemibrain",
                match = "FAFB",
                id = "bodyid")
  }

  ###################
  # Hemibrain -> LM #
  ###################
  missing = is.na(lmg$hemibrain.match)
  matches = hg$bodyid[match(lmg$id[missing],hg$LM.match)]
  quality = hg$LM.match.quality[match(lmg$id[missing],hg$LM.match)]
  lmg[missing,"hemibrain.match"] = matches
  lmg[missing,"hemibrain.quality"] = quality
  sorted = which(missing)[!is.na(matches)]
  if(length(sorted)){
    update_gsheet(update = lmg[sorted,],
                gs = lmg,
                tab = "lm",
                match = "hemibrain",
                id = "id")
  }

  #####################
  # Hemibrain -> FAFB #
  #####################
  missing = is.na(fg$hemibrain.match)
  matches = hg$bodyid[match(fg$skid[missing],hg$FAFB.match)]
  quality = hg$FAFB.match.quality[match(fg$skid[missing],hg$FAFB.match)]
  fg[missing,"hemibrain.match"] = matches
  fg[missing,"hemibrain.match.quality"] = quality
  sorted = which(missing)[!is.na(matches)]
  if(length(sorted)){
    update_gsheet(update = fg[sorted,],
                gs = fg,
                tab = "FAFB",
                match = "hemibrain",
                id = "skid")
  }
}


# Udate function
update_gsheet <- function(update,
                          gs,
                          selected_file = options()$hemibrainr_matching_gsheet,
                          tab,
                          match = c("hemibrain", "LM", "FAFB", "flywire"),
                          id){
  match = match.arg(match)
  if(match=="flywire"){
    match2 = "FAFB"
  }else{
    match2 = match
  }
  for(row in 1:nrow(update)){
    columns = c(paste0(match,".match"), paste0(match2,".match.quality"))
    r = match(update[row,id],gs[[id]])+1
    if(is.issue(r)){
      next
    }
    for(column in columns){
      letter = LETTERS[match(column,colnames(gs))]
      range = paste0(letter,r)
      gsheet_manipulation(FUN = googlesheets4::range_write,
                          ss = selected_file,
                          range = range,
                          data = as.data.frame(update[row,column]),
                          sheet = tab,
                          col_names = FALSE)
    }
  }

}


#' @rdname hemibrain_add_made_matches
#' @export
fafb_matching_rewrite <- function(selected_file  = options()$hemibrainr_matching_gsheet,
                                   ...){
  matches = hemibrain_matches()
  n = hemibrain_match_sheet(sheet = "FAFB", selected_file = selected_file)
  n1 = elmr::fafb_get_meta("annotation:Lineage_annotated", batch = TRUE, ...)
  n2 = subset(n1, n1$skid %in% n$skid)
  n[match(n2$skid,n$skid),c("skid","ItoLee_Hemilineage", "Hartenstein_Hemilineage", "cell_body_fiber")] = n2[,c("skid","ItoLee_Hemilineage", "Hartenstein_Hemilineage", "cell_body_fiber")]
  ids.missing = setdiff(n1$skid,n$skid)
  n3 = elmr::fafb_get_meta(unique(ids.missing), batch = TRUE, ...)
  n = plyr::rbind.fill(n, n3[,c("skid","ItoLee_Hemilineage", "Hartenstein_Hemilineage", "cell_body_fiber")])
  n = n[!duplicated(n),]
  n$cell.type = matches[as.character(n$skid),"connectivity.type"]
  lskids = as.character(catmaid::catmaid_skids("annotation:side: left", ...))
  n$side = "right"
  n[n$skid%in%lskids,"side"] = "left"
  n$User[is.na(n$User)] = "flyconnectome"
  nblast = tryCatch(hemibrain_nblast('hemibrain-fafb14'), error = function(e) NULL)
  if(!is.null(nblast)){
    nblast = hemibrain_nblast('hemibrain-fafb14')
    nblast.top =nblast[match(n$skid,rownames(nblast)),]
    tops = apply(nblast.top,1,function(r) which.max(r))
    top = colnames(nblast)[unlist(tops)]
    top[!n$skid%in%rownames(nblast)] = NA
    n$nblast.top = top
  }
  n = n[order(n$cell.type),]
  n = n[order(n$ItoLee_Hemilineage),]
  gsheet_manipulation(data = n[0,],
                      ss = selected_file,
                      sheet = "FAFB")
  batches = split(1:nrow(n), ceiling(seq_along(1:nrow(n))/500))
  for(i in batches){
    gsheet_manipulation(FUN = googlesheets4::sheet_append,
                                     data = n[min(i):max(i),],
                                     ss = selected_file,
                                     sheet = "FAFB")
  }
  missing = setdiff(subset(matches,matches$dataset=="hemibrain")$match,subset(matches,matches$dataset=="FAFB")$id)
  missing = unique(missing[!missing%in%c("none","","NA"," ","good","medium","poor","tract")])
  hemibrain_matching_add(ids = missing, sheet="FAFB", ...)
}


#' @rdname hemibrain_add_made_matches
#' @export
hemibrain_matching_rewrite <- function(ids = NULL,
                                       selected_file  = options()$hemibrainr_matching_gsheet,
                                  ...){
  gs = hemibrain_match_sheet(sheet = "hemibrain", selected_file = selected_file)
  if(is.null(ids)){
    ids = hemibrain_neuron_bodyids()
  }
  meta1 = hemibrain_get_meta(unique(ids), ...)
  ids.missing = setdiff(gs$bodyid,meta1$bodyid)
  meta2 = hemibrain_get_meta(unique(ids.missing), ...)
  meta = rbind(meta1,meta2)
  meta$cell.type = meta$type
  meta = meta[,c("bodyid","ItoLee_Hemilineage","Hartenstein_Hemilineage","cellBodyFiber","cell.type","layer","ct.layer")]
  meta$FAFB.match = gs$FAFB.match[match(meta$bodyid,gs$bodyid)]
  meta$FAFB.match.quality = gs$FAFB.match.quality[match(meta$bodyid,gs$bodyid)]
  meta$LM.match = gs$LM.match[match(meta$bodyid,gs$bodyid)]
  meta$LM.match.quality = gs$LM.match.quality[match(meta$bodyid,gs$bodyid)]
  meta$User = gs$User[match(meta$bodyid,gs$bodyid)]
  meta = meta[order(meta$bodyid),]
  meta = meta[order(meta$cell.type),]
  meta = meta[order(meta$ItoLee_Hemilineage),]
  batches = split(1:nrow(meta), ceiling(seq_along(1:nrow(meta))/500))
  gsheet_manipulation(FUN = googlesheets4::write_sheet,
                      data = meta[0,],
                      ss = selected_file,
                      sheet = "hemibrain")
  for(i in batches){
    gsheet_manipulation(FUN = googlesheets4::sheet_append,
                        data = meta[min(i):max(i),],
                        ss = selected_file,
                        sheet = "hemibrain")
  }
}

# hidden
id_selector <- function(gs,
                        ids = NULL,
                        id = c("bodyid","flywire.id","skid","id","flywire.xyz"),
                        overwrite = c("FALSE","mine","mine_empty","TRUE"),
                        quality.field,
                        match.field,
                        initials = NULL,
                        column = NULL,
                        field = NULL){
  id = match.arg(id)
  overwrite = as.character(overwrite)
  overwrite = match.arg(overwrite)
  if(is.null(column)+is.null(field)==1){
    stop("column and field must both be NULL, or both be given")
  }
  # choose possible ids
  if(is.null(ids)|!length(ids)){
    if(!is.null(initials) & overwrite %in% c("FALSE","mine")){
      ids = gs[[id]][gs$User==initials]
    }else if(!is.null(initials) & overwrite == "mine_empty"){
      ids = gs[[id]][gs$User==initials|is.na(gs[[match.field]])|is.na(gs[[quality.field]])]
    }
  }else{
    ids = intersect(ids,gs[[id]])
  }
  ids = ids[!grepl("missing",ids)]
  # further narrow
  selected = subset(gs, gs[[id]] %in% ids)
  if(overwrite %in% c("mine", "mine_empty")){
    doit = subset(selected, is.na(selected[[quality.field]]) | selected[[quality.field]] %in% c("none","n","tract","t","","NA"))
  } else if(overwrite=="TRUE"){
    doit = subset(selected, is.na(selected[[match.field]]) |  is.na(selected[[match.field]])%in%c(""," "))
  }else{
    doit = selected
  }
  # choose wheee we only have certain column values
  if(is.null(column)+is.null(field)==0){
    selected = subset(selected, selected[[column]]%in%field)
  }
  # return
  selected
}



