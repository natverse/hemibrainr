#####################################################################################
################################ FIB-FAFB Matching ##################################
#####################################################################################

#' Match up neurons between the hemibrain and FAFB
#'
#' @description Match up neurons between the hemibrain, FAFB and light level data and save the
#'   result using a
#'   \href{https://docs.google.com/spreadsheets/d/1OSlDtnR3B1LiB5cwI5x5Ql6LkZd8JOS5bBr-HTi0pOw/edit#gid=0}{Google
#'   Sheet} on the hemibrain Google Team Drive operated by the flyconnectome
#'   group at the University of Cambridge. This function
#'   only considers lateral horn neurons. You must have access to the Team Drive
#'   in order to use this function. This function makes use of the Google
#'   Filestream application, which should be installed on your machine. Further,
#'   note that neurons are read from the FAFB CATMAID project, and you must have
#'   log-in details for this project recorded in your .Renviron for this
#'   function to work.
#'
#' @param ids body IDs for hemibrain neurons present in the
#'   \href{https://docs.google.com/spreadsheets/d/1OSlDtnR3B1LiB5cwI5x5Ql6LkZd8JOS5bBr-HTi0pOw/edit#gid=0}{Google
#'    Sheet}, for which the user will attempt to make a match if one has not
#'   been made already. Else, LM neuron IDs in the tab 'lm' when using \code{lm_matching}.
#' @param hemibrain.nblast a FAFB (rows) - hemibrain (columns) normalised
#'   NBLAST matrix. By default this is read from the flyconnectome Team Drive.
#' @param selected_file the Google Sheet database to read and write from. For
#'   now, defaults to a
#'   \href{https://docs.google.com/spreadsheets/d/1OSlDtnR3B1LiB5cwI5x5Ql6LkZd8JOS5bBr-HTi0pOw/edit#gid=0}{Google
#'    Sheet for lateral horn neurons}. No other databases have been prepared.
#' @param batch_size the number of FAFB top matches to read from CATMAID in one
#'   go.
#' @param db Either a neuronlist or the name of a character vector naming a
#'   neuronlist. Defaults to the value of \code{\link{hemibrain_neurons}()}.
#' @param match.type whether to match up FAFB skeletons (\code{"EM"}) or light level
#' skeletons  (\code{"LM"}) from Dolan et al. and Frechter et al. 2019 (eLife), stored in the package \code{lhns} as
#' \code{most.lhns}.
#' @param query a neuronlist of light level neurons to match against. Should correspond to the given NBLAST matrix.
#' Defaults to reading a transformed \code{most.lhns} from the Hemibrain Google Team Drive.
#' @param overwrite logical, whether or not to overwrite matches already made.
#'
#' @details Currently, the
#'   \href{https://docs.google.com/spreadsheets/d/1OSlDtnR3B1LiB5cwI5x5Ql6LkZd8JOS5bBr-HTi0pOw/edit#gid=0}{Google
#'    Sheet} is set up with  limited number of users, each of whom have been
#'   assigned a number of neurons to match up. In order to add yourself as a
#'   user, simply open this Google Sheet in your browser and add your initials
#'   to neurons of your choosing on the rightmost column 'Users'. Once a match
#'   is recorded, the user selects a quality for that match. There can be no
#'   match (n), a poor match (p) an okay match (o) or an exact match (e). As a
#'   rule of thumb, a poor match could be a neuron from a very similar same cell
#'   type or a highly untraced neuron that may be the correct cell type. An okay
#'   match should be a neuron that looks to be from the same morphological cell
#'   type but there may be some discrepancies in its arbour. An exact match is a
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
#' hemibrain_FAFB_matching(hemibrain.nblast = fib.fafb.crossnblast.twigs5.mean.compress)
#' }}
#' @rdname hemibrain_matching
#' @export
#' @seealso \code{\link{hemibrain_adjust_saved_split}}
hemibrain_matching <- function(ids = NULL,
                         hemibrain.nblast = NULL,
                         selected_file = "1OSlDtnR3B1LiB5cwI5x5Ql6LkZd8JOS5bBr-HTi0pOw",
                         batch_size = 10,
                         db=hemibrain_neurons(),
                         match.type = c("FAFB", "LM"),
                         query = NULL,
                         overwrite = FALSE){
  match.type = match.arg(match.type)
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
  if(is.null(hemibrain.nblast) & match.type =="FAFB"){
    matname="fib.fafb.crossnblast.twigs5.mean.compress"
      message("Loading FAFB-FIB NBLAST ", matname,
              " from flyconnectome Google Team Drive using Google Filestream: ")
      load(sprintf("/Volumes/GoogleDrive/Shared drives/flyconnectome/fafbpipeline/%s.rda", matname))
      hemibrain.nblast = get(matname)
      rm("fib.fafb.crossnblast.twigs5.mean.compress")
  }
  if(is.null(hemibrain.nblast) & match.type =="LM"){
    matname="hemibrain.lhns.mean.compressed"
      message("Loading LM-FIB NBLAST ", matname,
              " from hemibrain Google Team Drive using Google Filestream: ")
      load(sprintf("/Volumes/GoogleDrive/Shared\ drives/hemibrain/hemibrain_nblast/%s.rda", matname))
      hemibrain.nblast = get(matname)
      hemibrain.nblast = t(hemibrain.nblast)
      rm("hemibrain.lhns.mean.compressed")
  }
  # Read the Google Sheet
  gs = gsheet_manipulation(FUN = googlesheets4::read_sheet,
                           ss = selected_file,
                           sheet = "lhns",
                           guess_max = 3000,
                           return = TRUE)
  gs$bodyid = correct_id(gs$bodyid)
  rownames(gs) = gs$bodyid
  # Read neuron meta data
  meta = neuprintr::neuprint_get_meta(ids)
  meta = meta[order(meta$type),]
  meta = subset(meta, meta$bodyid%in%ids)
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
  if(match.type=="LM" & is.null(query)){
    fafb = hemibrain_lm_lhns(brainspace = c("JRCFIB2018F"))
  }else if(match.type=="LM"){
    fafb = query
  }
  # How much is done?
  match.field = paste0(match.type,".match")
  quality.field = paste0(match.type,".match.quality")
  done = subset(gs, !is.na(gs[[match.field]]))
  message("Neuron matches: ", nrow(done))
  print(table(gs[[quality.field]]))
  # Choose user
  message("Users: ", paste(unique(gs$User),collapse = " "))
  initials = must_be("Choose a user : ", answers = unique(gs$User))
  say_hello(initials)
  rgl::bg3d("white")
  # choose ids
  if(is.null(ids)){
    ids = gs$bodyid[gs$User==initals]
  }else{
    ids = intersect(ids,gs$bodyid)
  }
  # choose brain
  if(match.type=="FAFB"){
    brain = elmr::FAFB14.surf
  }else{
    brain = hemibrainr::hemibrain_microns.surf
  }
  # Make matches!
  for(n in meta$bodyid){
    # Get bodyid
    n = as.character(n)
    # Remove neurons with matches
    if(!overwrite){
      donotdo = subset(gs, !is.na(gs[[match.field]]) | !bodyid%in%ids)
    }else{
      donotdo = subset(gs, !bodyid%in%ids)
    }
    if(n%in%donotdo$bodyid | !n%in%unlist(dimnames(hemibrain.nblast))){
      next
    }
    # Plot brain
    rgl::clear3d()
    plot3d(brain, alpha = 0.1, col ="grey")
    # Read hemibrain neuron
    if(is.null(db)){
      lhn  = neuprintr::neuprint_read_neurons(n)
    } else {
      lhn = tryCatch(db[as.character(n)], error = function(e) {
        warning("Cannot read neuron: ", n, " from local db; fetching from neuPrint!")
        neuprintr::neuprint_read_neurons(n)
        })
    }
    # Transform hemibrain neuron to FAFB space
    lhn = scale_neurons.neuronlist(lhn, scaling = (8/1000))
    message("Hemibrain body ID: ", lhn[n,"bodyid"])
    message("Hemibrain-assigned cell type : ",lhn[n,"type"])
    # Read top 10 FAFB matches
    if(match.type=="FAFB"){
      lhn = suppressWarnings(nat.templatebrains::xform_brain(lhn, reference = "FAFB14", sample = "JRCFIB2018F"))
      message(sprintf("Reading the top %s FAFB matches",batch_size))
      r = tryCatch(sort(hemibrain.nblast[,as.character(n)],decreasing = TRUE), error = function(e) NULL)
      if(is.null(r)){
        r = sort(hemibrain.nblast[as.character(n),],decreasing = TRUE)
      }
      fafb = catmaid::read.neurons.catmaid(names(r)[1:batch_size])
      j = batch_size
    }else{
      r = sort(hemibrain.nblast[,n],decreasing = TRUE)
      j = length(fafb)
    }
    plot3d(lhn[n], lwd = 2, soma = TRUE, col = "black")
    sel = c("go","for","it")
    k = 1
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
        if(match.type=="FAFB"){
          prog = hemibrain_choice(sprintf("Do you want to read %s more neurons from CATMAID? ", batch_size))
          if(prog){
            k = j
            j = j + batch_size
            fafb = nat::union(fafb, catmaid::read.neurons.catmaid(names(r)[(k+1):j]))
          }
        }
      }
    }
    # Assign match and its quality
    gs[n,match.field] = ifelse(length(sel)==0,'none',sel)
    if(length(sel)){
      rgl::plot3d(fafb[sel],col="blue",lwd=2,soma=TRUE)
      quality = must_be("What is the quality of this match? good(e)/okay(o)/poor(p) ", answers = c("e","o","p"))
    }else{
      quality = "n"
    }
    quality = standardise_quality(quality)
    gs[n,quality.field] = quality
    gs = gs[!duplicated(gs$bodyid),]
    unsaved = c(unsaved, n)
    message(length(unsaved), " unsaved matches")
    print(knitr::kable(gs[unsaved,c("bodyid","type",match.field,quality.field)]))
    p = must_be("Continue (enter) or save (s)? ", answers = c("","s"))
    if(p=="s"){
      plot_inspirobot()
      say_encouragement(initials)
      # Read!
      gs2 = gsheet_manipulation(FUN = googlesheets4::read_sheet,
                               ss = selected_file,
                               sheet = "lhns",
                               guess_max = 3000,
                               return = TRUE)
      gs2$bodyid = correct_id(gs2$bodyid)
      rownames(gs2) = gs2$bodyid
      gs = gs[rownames(gs2),]
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
}

# hidden
## Google sheet write function
write_matches <- function(gs,
                       ids = NULL,
                       id.field = "bodyid",
                       ws = "lhns",
                       selected_file = "1OSlDtnR3B1LiB5cwI5x5Ql6LkZd8JOS5bBr-HTi0pOw",
                       column = colnames(gs)){
  if((!identical(colnames(gs),column) & length(column) > 1)|(sum(column%in%colnames(gs))<1)){
    stop("Column must be one column of the google sheet, or all the columns")
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

# hidden
correct_id <-function(v){
  gsub(" ","",v)
}


#' @rdname hemibrain_matching
#' @export
lm_matching <- function(ids = NULL,
                        hemibrain.nblast = NULL,
                        selected_file = "1OSlDtnR3B1LiB5cwI5x5Ql6LkZd8JOS5bBr-HTi0pOw",
                        batch_size = 50,
                        db=hemibrain_neurons(),
                        query = hemibrain_lm_lhns(brainspace = c("JRCFIB2018F")),
                        overwrite = FALSE){
  # Motivate!
  nat::nopen3d()
  plot_inspirobot()
  unsaved = c()
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
    matname="hemibrain.lhns.mean.compressed"
    if(exists(matname)) {
      message("Using loaded LM-FIB NBLAST: ", matname)
      hemibrain.nblast = get(matname)
    } else {
      message("Loading LM-FIB NBLAST ", matname,
              " from hemibrain Google Team Drive using Google Filestream: ")
      load(sprintf("/Volumes/GoogleDrive/Shared\ drives/hemibrain/hemibrain_nblast/%s.rda", matname))
      hemibrain.nblast = get(matname)
      rm("hemibrain.lhns.mean.compressed")
    }
  }
  # Read the Google Sheet
  gs = gsheet_manipulation(FUN = googlesheets4::read_sheet,
                           ss = selected_file,
                           sheet = "lm",
                           guess_max = 3000,
                           return = TRUE)
  gs$id = correct_id(gs$id)
  rownames(gs) = gs$id
  if(is.null(ids)){
    ids = gs$id
  }else{
    ids = intersect(ids,gs$id)
  }
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
  message("Users: ", paste(unique(gs$User),collapse = " "))
  initials = must_be("Choose a user : ", answers = unique(gs$User))
  say_hello(initials)
  rgl::bg3d("white")
  # choose brain
  brain = hemibrain.surf
  # Make matches!
  for(n in gs$id){
    # Get id
    n = as.character(n)
    end = n==gs$id[length(gs$id)]
    # Remove neurons with matches
    if(!overwrite){
      donotdo = subset(gs, !is.na(gs[[match.field]]) | User != initials | !id%in%ids)
    }else{
      donotdo = subset(gs, User != initials | !id%in%ids)
    }
    if(n%in%donotdo$id | !n%in%names(query)){
      next
    }
    # Plot brain
    rgl::clear3d()
    plot3d(brain, alpha = 0.1, col ="grey")
    # Transform hemibrain neuron to FAFB space
    lhn = query[n]
    lhn = scale_neurons(lhn, scaling = (1000/8))
    plot3d(lhn, lwd = 2, soma = TRUE, col = "black")
    message("ID: ", n)
    message("cell type : ",lhn[n,"cell.type"])
    # Read top 10 FAFB matches
    r = sort(hemibrain.nblast[,n],decreasing = TRUE)
    message(sprintf("Reading the top %s hemibrain hits",batch_size))
    r = sort(hemibrain.nblast[,n],decreasing = TRUE)
    # Read hemibrain neurons
    if(is.null(db)){
      hemi  = neuprintr::neuprint_read_neurons((names(r)[1:batch_size]))
    } else {
      batch = names(r)[1:batch_size]
      batch.in = intersect(batch, names(db))
      hemi = tryCatch(db[match(batch.in,names(db))], error = function(e) NULL)
      if(is.null(hemi)|length(batch.in)!=length(batch)){
        message("Cannot read neuron: ", n, " from local db; fetching from neuPrint!")
        batch.out = setdiff(batch, names(hemi))
        hemi = c(hemi,neuprintr::neuprint_read_neurons(batch.out))
        hemi = hemi[as.character(batch)]
      }
    }
    sel = c("go","for","it")
    k = 1
    j = batch_size
    # Cycle through potential matches
    while(length(sel)>1){
      sel = sel.orig = nat::nlscan(hemi[names(r)[1:j]], col = "red", lwd = 2, soma = TRUE)
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
              hemi2  = neuprintr::neuprint_read_neurons((names(r)[(k+1):j]))
            } else {
              hemi2 = tryCatch(db[(names(r)[(k+1):j])], error = function(e) {
                warning("Cannot read neuron: ", n, " from local db; fetching from neuPrint!")
                neuprintr::neuprint_read_neurons((names(r)[(k+1):j]))
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
      quality = must_be("What is the quality of this match? good(e)/okay(o)/poor(p) ", answers = c("e","o","p"))
    }else{
      quality = "n"
    }
    quality = standardise_quality(quality)
    gs[n,quality.field] = quality
    gs = gs[!duplicated(gs$id),]
    unsaved = c(unsaved, n)
    message(length(unsaved), " unsaved matches")
    print(knitr::kable(gs[unsaved,c("id","type",match.field,quality.field)]))
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
  say_encouragement(initials)
}


#' @rdname hemibrain_matching
#' @export
fafb_matching <- function(ids = NULL,
                        hemibrain.nblast = NULL,
                        selected_file = "1OSlDtnR3B1LiB5cwI5x5Ql6LkZd8JOS5bBr-HTi0pOw",
                        batch_size = 50,
                        db=hemibrain_neurons(),
                        query = NULL,
                        overwrite = FALSE){
  # Motivate!
  nat::nopen3d()
  plot_inspirobot()
  unsaved = c()
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
    matname="fib.fafb.crossnblast.twigs5.mean.compress"
    message("Loading FAFB-FIB NBLAST ", matname,
            " from flyconnectome Google Team Drive using Google Filestream: ")
    load(sprintf("/Volumes/GoogleDrive/Shared drives/flyconnectome/fafbpipeline/%s.rda", matname))
    hemibrain.nblast = get(matname)
    rm("fib.fafb.crossnblast.twigs5.mean.compress")
  }
  # Read the Google Sheet
  gs = gsheet_manipulation(FUN = googlesheets4::read_sheet,
                           ss = selected_file,
                           sheet = "fafb",
                           guess_max = 3000,
                           return = TRUE)
  gs$skid = correct_id(gs$skid)
  rownames(gs) = gs$skid
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
  message("Users: ", paste(unique(gs$User),collapse = " "))
  initials = must_be("Choose a user : ", answers = unique(gs$User))
  say_hello(initials)
  rgl::bg3d("white")
  # choose ids
  if(is.null(ids)){
    ids = gs$skid[gs$User==initals]
  }else{
    ids = intersect(ids,gs$skid)
  }
  # choose brain
  brain = hemibrain.surf
  # Make matches!
  for(n in gs$skid){
    # Get id
    n = as.character(n)
    end = n==gs$skid[length(gs$skid)]
    # Remove neurons with matches
    if(!overwrite){
      donotdo = subset(gs, !is.na(gs[[match.field]]) | !skid%in%ids)
    }else{
      donotdo = subset(gs, !skid%in%ids)
    }
    if(n%in%donotdo$skid){
      next
    }
    # Plot brain
    rgl::clear3d()
    plot3d(brain, alpha = 0.1, col ="grey")
    # Transform hemibrain neuron to FAFB space
    if(!is.null(query)){
      lhn = query[n]
    }else{
      lhn = catmaid::read.neurons.catmaid(n)
    }
    lhn = suppressWarnings(nat.templatebrains::xform_brain(lhn, sample = "FAFB14", reference = "JRCFIB2018F"))
    lhn = scale_neurons(lhn, scaling = (1000/8))
    plot3d(lhn, lwd = 2, soma = TRUE, col = "black")
    message("SKID: ", n)
    message("ItoLee_Hemilineage : ",lhn[n,"ItoLee_Hemilineage"])
    # Read top 10 FIB matches
    r = sort(hemibrain.nblast[n,],decreasing = TRUE)
    message(sprintf("Reading the top %s hemibrain hits",batch_size))
    # Read hemibrain neurons
    if(is.null(db)){
      hemi  = neuprintr::neuprint_read_neurons((names(r)[1:batch_size]))
    } else {
      batch = names(r)[1:batch_size]
      batch.in = intersect(batch, names(db))
      hemi = tryCatch(db[match(batch.in,names(db))], error = function(e) NULL)
      if(is.null(hemi)|length(batch.in)!=length(batch)){
        message("Cannot read neuron: ", n, " from local db; fetching from neuPrint!")
        batch.out = setdiff(batch, names(hemi))
        hemi = c(hemi,neuprintr::neuprint_read_neurons(batch.out))
        hemi = hemi[as.character(batch)]
      }
    }
    sel = c("go","for","it")
    k = 1
    j = batch_size
    # Cycle through potential matches
    while(length(sel)>1){
      sel = sel.orig = nat::nlscan(hemi[names(r)[1:j]], col = "red", lwd = 2, soma = TRUE)
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
            hemi2  = neuprintr::neuprint_read_neurons((names(r)[(k+1):j]))
          } else {
            hemi2 = tryCatch(db[(names(r)[(k+1):j])], error = function(e) {
              warning("Cannot read neuron: ", n, " from local db; fetching from neuPrint!")
              neuprintr::neuprint_read_neurons((names(r)[(k+1):j]))
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
      quality = must_be("What is the quality of this match? good(e)/okay(o)/poor(p) ", answers = c("e","o","p"))
    }else{
      quality = "n"
    }
    quality = standardise_quality(quality)
    gs[n,quality.field] = quality
    gs = gs[!duplicated(gs$skid),]
    unsaved = c(unsaved, n)
    message(length(unsaved), " unsaved matches")
    print(knitr::kable(gs[unsaved,c("skid","ItoLee_Hemilineage",match.field,quality.field)]))
    p = must_be("Continue (enter) or save (s)? ", answers = c("","s"))
    if(p=="s"|end){
      plot_inspirobot()
      say_encouragement(initials)
      # Read!
      gs2 = gsheet_manipulation(FUN = googlesheets4::read_sheet,
                                ss = selected_file,
                                sheet = "fafb",
                                guess_max = 3000,
                                return = TRUE)
      gs2$skid = correct_id(gs2$skid)
      rownames(gs2) = gs2$skid
      gs = gs[rownames(gs2),]
      gs = gs[!duplicated(gs$skid),]
      # Write!
      write_matches(gs=gs,
                    ids = unsaved,
                    id.field = "skid",
                    column = match.field,
                    ws = "fafb")
      write_matches(gs=gs,
                    ids = unsaved,
                    id.field = "skid",
                    column = quality.field,
                    ws = "fafb")
      unsaved = c()
      gs = gs2
      rgl::bg3d("white")
    }
  }
  say_encouragement(initials)
}


#' Retrieve matched up neurons between the hemibrain and FAFB
#'
#' @description Many neurons in the hemibrain data have been matched one to one with FAFB neurons (hemibrain->FAFB) and
#' FAFB neurons have been matched with their hemibrain counterparts (FAFB->hemibrain). These matchs have been done by the
#' Flyconnectome Group at the University of Cambridge, and are recorded on a
#' \href{https://docs.google.com/spreadsheets/d/1OSlDtnR3B1LiB5cwI5x5Ql6LkZd8JOS5bBr-HTi0pOw/edit#gid=0}{Google
#'   Sheet} on our Hemibrain Google Team Drive, to which you will need access
#' through an authenticated account to view and use this function. Matches have three levels of 'quality', largely
#' dependent on the degree of manual tracing for FAFB neurons - good (could be the same cell), medium (same cell type) and poor (could be the same or similar cell type).
#'
#' @param priority whether to use FAFB->hemibrain matches (FAFB) or hemibrain->FAFB matches (hemibrain) in order to ascribe
#' cell type names to FAFB neurons. In both cases, cell type names are attached to hemibrain bodyids, and propogated to their FAFB matches.
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
hemibrain_matches <- function(priority = c("FAFB","hemibrain")){
  priority = match.arg(priority)

  # Get matches
  selected_file = "1OSlDtnR3B1LiB5cwI5x5Ql6LkZd8JOS5bBr-HTi0pOw"
  hemibrain.matches = hemibrainr:::gsheet_manipulation(FUN = googlesheets4::read_sheet,
                                                       ss = selected_file,
                                                       sheet = "lhns",
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

  # Get FAFB matches
  fafb.matches = hemibrainr:::gsheet_manipulation(FUN = googlesheets4::read_sheet,
                                                  ss = selected_file,
                                                  sheet = "fafb",
                                                  return = TRUE)
  fafb.matches$skid = correct_id(fafb.matches$skid)
  fafb.matches = fafb.matches[!duplicated(fafb.matches$skid),]
  fafb.matches = fafb.matches[fafb.matches$skid!="",]
  fafb.matches = subset(fafb.matches, !is.na(fafb.matches$skid))
  fafb.matches = fafb.matches[!duplicated(fafb.matches$skid),]
  fafb.matches$dataset = "FAFB"
  rownames(fafb.matches) = fafb.matches$skid

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
  ntotype = hemibrain.matches$bodyid[is.na(hemibrain.matches$cell.type)]
  meta = neuprintr::neuprint_get_meta(ntotype)
  types = meta$type
  names(types) = meta$bodyid
  hemibrain.matches[names(types),"cell.type"] = types

  # Work out lineages
  for(id in as.character(fafb.matches$skid)){
    if(is.na(id)){
      next
    }
    ct = fafb.matches[id,"cell.type"]
    if(is.na(ct)){
      fafb.matches[id,"cell.type"] = hemibrain.matches$cell.type[match(id,hemibrain.matches$FAFB.match)]
      fafb.matches[id,"ItoLee_Hemilineage"] = hemibrain.matches$ItoLee_Hemilineage[match(id,hemibrain.matches$FAFB.match)]
    }else{
      fafb.matches[id,"ItoLee_Hemilineage"] = hemibrain.matches$ItoLee_Hemilineage[match(ct,hemibrain.matches$cell.type)]
    }
  }

  # Rename cells
  fafb.matches$cell = paste0(fafb.matches$cell.type,"#",ave(fafb.matches$cell.type,fafb.matches$cell.type,FUN= seq.int))
  hemibrain.matches$cell = paste0(hemibrain.matches$cell.type,"#",ave(hemibrain.matches$cell.type,hemibrain.matches$cell.type,FUN= seq.int))

  # Fix hemilineages
  hl = hemibrain.matches$ItoLee_Hemilineage[match(fafb.matches$skid,hemibrain.matches$FAFB.match)]
  l = hemibrain.matches$ItoLee_Lineage[match(fafb.matches$skid,hemibrain.matches$FAFB.match)]
  fafb.matches$ItoLee_Hemilineage[!is.na(hl)] = hl[!is.na(hl)]
  fafb.matches$ItoLee_Lineage[!is.na(l)] = l[!is.na(l)]

  # Make matching data frame
  matched.h = hemibrain.matches[,c("bodyid", "cell.type", "cell", "ItoLee_Hemilineage",
                                   "FAFB.match", "FAFB.match.quality", "LM.match", "LM.match.quality", "dataset")]
  matched.f = fafb.matches[,c("skid",  "cell.type",  "cell", "ItoLee_Hemilineage",
                              "hemibrain.match", "hemibrain.match.quality", "LM.match", "LM.match.quality","dataset")]
  colnames(matched.h) = colnames(matched.f) = c("id","cell.type", "cell","ItoLee_Hemilineage","match","quality", "LM.match", "LM.match.quality","dataset")
  matched = rbind(matched.h,matched.f)
  matched$quality[is.na(matched$match)] = "none"
  matched$match[is.na(matched$match)] = "none"

  # Sort out types
  matched$connectivity.type = matched$cell.type
  matched$cell.type = gsub("[a-z]$","",matched$cell.type)
  matched$cell.type[is.na(matched$cell.type)] = "uncertain"

  # Return
  matched

}


