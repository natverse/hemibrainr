#####################################################################################
################################ FIB-FAFB Matching ##################################
#####################################################################################

#' Match up neurons between the hemibrain and FAFB
#'
#' @description Match up neurons between the hemibrain, FAFB and light level data and save the
#'   result using a
#'   \href{https://docs.google.com/spreadsheets/d/1OSlDtnR3B1LiB5cwI5x5Ql6LkZd8JOS5bBr-HTi0pOw/edit#gid=0}{Google
#'   Sheet} on the hemibrain Google Team Drive operated by the flyconnectome
#'   group at the University of Cambridge. Your aim is to match the light blue neuron to the best red neuron!
#'   You must have access to the Team Drive
#'   in order to use this function. This function makes use of the Google
#'   Filestream application or rclone, which should be installed on your machine. Further,
#'   note that neurons are read from the FAFB CATMAID project when \code{repository=="CATMAID"}, and you must have
#'   login details for this project recorded in your .Renviron for these
#'   functions to work. You are given neurons assigned to your initials on the matching google sheet (the \code{user} column)
#'   though you can look at other's assigned matches using \code{superuser=TRUE}. Your selection is further narrowed
#'   by omitting neurons whose matches have already been made, when using the default \code{overwrite="TRUE"}. However, you
#'   can change this to overwrite matches, and to review made matches.
#'
#' @param ids body IDs for hemibrain neurons present in the
#'   \href{https://docs.google.com/spreadsheets/d/1OSlDtnR3B1LiB5cwI5x5Ql6LkZd8JOS5bBr-HTi0pOw/edit#gid=0}{Google
#'    Sheet}, for which the user will attempt to make a match if one has not
#'   been made already. Else, LM neuron IDs in the tab 'lm' when using \code{lm_matching}.
#' @param hemibrain.nblast a FAFB (rows) - hemibrain (columns) normalised
#'   NBLAST matrix. By default this is read from the flyconnectome Team Drive.
#' @param mirror.nblast a flywire (rows) - flywire-mirrored (columns) normalised
#' NBLAST matrix. By default this is read from the hemibrain Team Drive.
#' @param threshold the minimum normalised NBLAST score between query and target neurons, for the potential match to be show. If set to \code{NULL} then all hits can be shown.
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
#' @param query a neuronlist of neurons for matching. Should correspond to the given NBLAST matrix.
#' Defaults to reading a transformed \code{most.lhns} from the Hemibrain Google Team Drive.
#' @param overwrite Whether or not to overwrite matches already made.
#' The neurons you could possibly be looking at are selected through
#' the arguments: \code{ids}, \code{column}, \code{entry}. If \code{ids} is \code{NULL} you will be given all neurons that have been assigned to your user on the
#' \href{https://docs.google.com/spreadsheets/d/1OSlDtnR3B1LiB5cwI5x5Ql6LkZd8JOS5bBr-HTi0pOw/edit#gid=0}{Google sheet}.
#' If \code{overwrite} is set to \code{"FALSE"}, you will not overwrite any matches that have already been made from among the selected neurons.
#' If \code{"TRUE"} (be careful!) then you overwrite made matches among the selected neurons.
#' If \code{"bad"} then you can overwrite made matches that are 'tract-only' or 'none'.
#' If \code{"review"} then you will only be shown already-matched neurons, which appear in green.
#' If you make no new selection (hit 't' to exit selection mode) then the made selection will persist.
#' If \code{ids} is not \code{NULL}, then the selected neurons will be further sub-setted by their unique id. Note that selection
#' also works with \code{superuser}, the default is to only take neurons allocated to you. In order to change this,
#' you can user \code{superuser=TRUE}.
#' @param column defaults to \code{NULL}, no further subsetting. Else, you can select a column from the Google sheet.
#' Only  neurons with a certain value (\code{entry}) in that column will be chosen for matching.
#' @param entry defaults to \code{NULL}, no further subsetting. Else, it is a value in \code{column}.
#' @param user the initials of the matching 'users', i.e. you, which should also be recorded on the master matching google sheet. However, you can enter a new user
#' for the matches you make in R. If set to \code{NULL} then all users are 'selected'.
#' @param superuser if \code{FALSE} then you will only be given neurons flagged for your user. If \code{TRUE} then
#' you will be given neurons flagged for any user. To select whether or not you want to look at neurons with no match, neurons with a match
#' or either, use the \code{overwrite} argument.
#' @param flywire.good logical, whether or not to only take 'well traced' flywire neurons, as annotated by the Drosophila Connectomics Group. This relies on the status column retrieved
#' by \code{flywire_meta}.
#' @param verbose logical. If \code{TRUE} the pipeline pauses for each neuron that is missing from the NBLAST/data set.
#'
#' @details Currently, the
#'   \href{https://docs.google.com/spreadsheets/d/1OSlDtnR3B1LiB5cwI5x5Ql6LkZd8JOS5bBr-HTi0pOw/edit#gid=0}{Google
#'    Sheet} is set up with  limited number of users, each of whom have been
#'   assigned a number of neurons to match up. In order to add yourself as a
#'   user, simply open this Google Sheet in your browser and add your initials
#'   to neurons of your choosing on the rightmost column 'users'. Once a match
#'   is recorded, the user selects a quality for that match. There can be no
#'   match (none), a poor match (poor) an okay match (okay) or an exact match (good). As a
#'   rule of thumb, a poor match could be a neuron from a very similar same cell
#'   type or a highly untraced neuron that may be the correct cell_type. An okay
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
#'
#' # And example of matching neurons in a flywire tracing sheet
#' sheet = flywire_tracing_sheet("SLPal2_dorsal", regex = TRUE)
#' fafb_matching(ids = unique(sheet$root_id), repository="flywire",
#' overwrite = "bad", user = "AJ", superuser = TRUE)
#'
#' # Look at all poorly made matches
#' fafb_matching(repository="flywire", column = "hemibrain_match_quality",
#' entry = "poor", overwrite = "TRUE", user = "AJ", superuser = TRUE)
#' }}
#' @rdname hemibrain_matching
#' @export
#' @seealso \code{\link{hemibrain_adjust_saved_split}}
hemibrain_matching <- function(ids = NULL,
                         hemibrain.nblast = NULL,
                         threshold = 0,
                         flytable = TRUE,
                         selected_file = options()$hemibrainr_matching_gsheet,
                         batch_size = 50,
                         db=NULL, # brain="FAFB"
                         repository = c("flywire", "CATMAID", "lm"),
                         query = hemibrain_neurons(brain = "FAFB14"), # brain="FAFB"
                         overwrite = c("FALSE","bad","TRUE","review"),
                         column = NULL,
                         entry = NULL,
                         user = NULL,
                         superuser = FALSE,
                         verbose = FALSE){
  repository = match.arg(repository)
  message("Matching hemibrain neurons (blue) to ", repository," neurons (red)")
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
  unsaved = saved = c()
  message("
          #######################Colours#######################
          blue = hemibrain neuron,
          red = potential FAFB matches based on NBLAST score,
          green = a chosen FAFB neuron during scanning,
          dark blue = your selected FAFB match.
          #######################Colours#######################
          ")
  ## Get NBLAST
  if(is.null(hemibrain.nblast) & repository == "CATMAID"){
    hemibrain.nblast = t(hemibrain_nblast("hemibrain-fafb14"))
  }
  if(is.null(hemibrain.nblast) & repository == "flywire"){
    hemibrain.nblast = t(hemibrain_nblast("hemibrain-flywire"))
  }
  if(is.null(hemibrain.nblast) & repository =="lm"){
    hemibrain.lhns.mean.compressed=NULL
    matname="hemibrain.lhns.mean.compressed"
    message("Loading LM-FIB NBLAST ", matname,
            " from hemibrain Google Team Drive using Google Filestream: ")
    load(sprintf("/Volumes/GoogleDrive/Shared\ drives/hemibrain/hemibrain_nblast/%s.rda", matname))
    hemibrain.nblast = get(matname)
    hemibrain.nblast = hemibrain.nblast
    rm("hemibrain.lhns.mean.compressed")
  }
  # Read the Google Sheet or seatable
  if(flytable){
    gs = flytable_matches(dataset="hemibrain")
  }else{
    gs = hemibrain_match_sheet(selected_file = selected_file, sheet = "hemibrain")
  }
  gs$user[is.na(gs$user)] = ""
  # Get hemibrain neurons
  if(missing(query)) {
    query=tryCatch(force(query), error=function(e) {
      message(e)
      message("Unable to use `hemibrain_neurons(brain = 'FAFB14')`. ",
                "I will read neurons from neuPrint, but this will be slower!")
      })
  }else if(is.character(query)) {
    query=tryCatch(get(query), error=function(e) stop("Unable to find neuronlist: ", query))
  }
  if(repository=="flywire" & is.null(db)){
    db=tryCatch(flywire_neurons(), error=function(e) {
      message(e)
      message("Unable to use `flywire_neurons()`. ",
              "I will read neurons from FlyWire, but this will be slower!")
    })
  }else if(repository=="LM"){
    db=tryCatch(lm_lhns(brainspace = c("FAFB14")), error=function(e) {
      message(e)
      message("Unable to use `lm_lhns()`. ",
              "I will read neurons from package lhns, but this will be slower!")
    })
  }else{
    db = NULL
  }
  # fields
  if(repository=="CATMAID"){
    match.field = search.id = "fafb_match"
    chosen.field = "skid"
    quality.field = "fafb_match_quality"
    brain = elmr::FAFB14.surf
  }else if(repository == "flywire"){
    search.id = "root_id"
    match.field = chosen.field = "flywire_xyz"
    quality.field = "fafb_match_quality"
    brain = elmr::FAFB14.surf
  }else{
    match.field = search.id = paste0(repository,".match")
    quality.field = paste0(repository,"_match_quality")
    brain = hemibrainr::hemibrain_microns.surf
  }
  id = "bodyid"
  # How much is done?
  done = subset(gs, !is.na(gs[[match.field]]))
  message("Neuron matches: ", nrow(done))
  print(table(gs[[quality.field]]))
  # Choose user
  initials = choose_user(gs, user = user)
  rgl::bg3d("white")
  # choose ids
  selected = id_selector(gs=gs, ids=ids, id = id, overwrite = overwrite,
                         quality.field = quality.field, match.field = match.field,
                         initials = initials, column = column, entry = entry, superuser = superuser)
  if(!verbose){
    rem = gs[[id]][!gs[[id]]%in%selected[[id]]]
    if(length(rem)){
      message("Removing ", length(rem), " neurons not in NBLAST matrix")
      selected = subset(selected, !selected[[id]]%in%rem)
    }
  }
  if(!nrow(selected)){
    stop("No neurons to match")
  }
  # Make matches!
  match.more = TRUE
  while(match.more){
    match_cycle = neuron_match_scanner(brain = brain,
                                       selected = selected,
                                       id = id,
                                       unsaved = unsaved,
                                       saved = saved,
                                       chosen.field = chosen.field,
                                       nblast = hemibrain.nblast,
                                       threshold = threshold,
                                       batch_size = batch_size,
                                       targets = db,
                                       targets.repository = repository,
                                       query = query,
                                       extra.neurons = NULL,
                                       query.repository = "hemibrain",
                                       extra.repository = "none",
                                       match.field = match.field,
                                       quality.field = quality.field,
                                       soma.size = 4000,
                                       show.columns = c("cell_type","ito_lee_hemilineage","status", match.field, quality.field,"note"),
                                       skip.if.absent = !verbose)
    selected = match_cycle[["selected"]]
    unsaved = match_cycle[["unsaved"]]
    if(length(unsaved)){
      plot_inspirobot()
      say_encouragement(initials)
      if(flytable){ # seatable
        flytable_matches_update(selected, dataset = "hemibrain")
      }else{ # google sheet
        # Read!
        gs2 = hemibrain_match_sheet(selected_file = selected_file, sheet = "hemibrain")
        selected.unsaved = subset(selected, selected[[id]]%in%unsaved)
        gs2[match(selected.unsaved[[id]],gs2[[id]]),match.field]= selected.unsaved[[match.field]]
        gs2[match(selected.unsaved[[id]],gs2[[id]]),quality.field]= selected.unsaved[[quality.field]]
        gs2[match(selected.unsaved[[id]],gs2[[id]]),"root_id"]= selected.unsaved[["root_id"]]
        gs2[match(selected.unsaved[[id]],gs2[[id]]),"note"]= selected.unsaved[["note"]]
        gs2[match(selected.unsaved[[id]],gs2[[id]]),"user"]= initials
        # Write!
        if(!identical(gs,gs2)){
          gsheet_update_cols(
            write.cols = c(match.field,quality.field,"note","user","root_id"),
            gs=gs2,
            selected_sheet = selected_file,
            sheet = "hemibrain",
            Verbose = TRUE)
        }
      }
      saved = c(unsaved, saved)
      unsaved = c()
    }
    match.more = hemibrain_choice("Match more neurons? y/n ")
  }
  say_encouragement(initials)
}

# hidden
## Google sheet write function
write_matches <- function(gs,
                          selected_file,
                       ids = NULL,
                       id.field = "bodyid",
                       ws = "hemibrain",
                       column = colnames(gs)){
  if((!identical(colnames(gs),column) & length(column) > 1)|(sum(column%in%colnames(gs))<1)){
    stop("Column must be one column of the Google sheet, or all the columns")
  }
  rows = (1:nrow(gs))+1
  all.ids = rownames(gs)
  rownames(gs) = rows
  if(!is.null(ids)){
    ids = intersect(ids, all.ids)
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
                        overwrite = c("FALSE","bad","TRUE","review"),
                        column = NULL,
                        entry = NULL,
                        user = NULL,
                        superuser = FALSE){
  # Motivate!
  nat::nopen3d()
  plot_inspirobot()
  message("
          #######################Colours##########################
          blue = LM neuron,
          red = potential hemibrain_matches based on NBLAST score,
          green = a chosen hemibrain neuron during scanning,
          dark blue = your selected hemibrain_match.
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
    q1 = lm_lhns(brainspace = c("JRCFIB2018F"))
    q2 = lm_lhns(cable = "lhins", brainspace = c("JRCFIB2018F"))
    q3 = lm_lhns(cable = "lines", brainspace = c("JRCFIB2018F"))
    query = nat::union(q1, q2)
    query = nat::union(query, q3)
  }
  # Read the Google Sheet
  gs = hemibrain_match_sheet(sheet = "hemibrain", selected_file = "lm")
  gs$user[is.na(gs$user)] = ""
  id = "id"
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
  quality.field = paste0("hemibrain","_match_quality")
  done = subset(gs, !is.na(gs[[match.field]]))
  message("Neuron matches: ", nrow(done), "/", nrow(gs))
  print(table(gs[[quality.field]]))
  # Choose user
  initials = choose_user(gs, user = user)
  rgl::bg3d("white")
  # choose ids
  selected = id_selector(gs=gs, ids=ids, id=id, overwrite = overwrite,
                         quality.field = quality.field, match.field = match.field,
                         initials = initials, column = column, entry = entry, superuser = superuser)
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
      points3d(nat::xyzmatrix(lhn), col = "#1BB6AF")
    }else{
      plot3d(lhn, lwd = 2, soma = TRUE, col = "#1BB6AF")
    }
    message("ID: ", n)
    message("cell_type : ",lhn[n,"cell_type"])
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
      sel = sel.orig = nat::nlscan(cycle, col = "#EE4244", lwd = 2, soma = TRUE)
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
    unsaved = unique(c(unsaved, n))
    message(length(unsaved), " unsaved matches")
    print(knitr::kable(gs[unsaved,c("id","cell_type",match.field,quality.field)]))
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
                    selected_file = selected_file,
                    ws = "lm")
      write_matches(gs=gs,
                    ids = unsaved,
                    id.field = "id",
                    column = quality.field,
                    selected_file = selected_file,
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
                  selected_file = selected_file,
                  ws = "lm")
    write_matches(gs=gs,
                  ids = unsaved,
                  id.field = "id",
                  column = quality.field,
                  selected_file = selected_file,
                  ws = "lm")
  }
  say_encouragement(initials)
}

#' @rdname hemibrain_matching
#' @export
fafb_matching <- function(ids = NULL,
                        repository = c("flywire","CATMAID"),
                        hemibrain.nblast = NULL,
                        threshold = 0,
                        selected_file = options()$hemibrainr_matching_gsheet,
                        batch_size = 20,
                        db=hemibrain_neurons(brain="FAFB14"),
                        query = NULL,
                        overwrite = c("FALSE","bad","TRUE","review"),
                        column = NULL,
                        entry = NULL,
                        user = NULL,
                        superuser = FALSE,
                        flywire.good = FALSE){
  repository = match.arg(repository)
  message("Matching ",repository," FAFB neurons (blue) to hemibrain neurons (red)")
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
  unsaved = saved = c()
  message("
          #######################Colours##########################
          black = FAFB CATMAID neuron you are trying to match,
          blue = flywire neuron you are trying to match,
          light grey = mirrrored flywire neuron,
          red = potential hemibrain_matches based on NBLAST score,
          green = a chosen hemibrain neuron during scanning,
          dark blue = your selected hemibrain_match.
          #######################Colours##########################
          ")
  ## Get NBLAST
  if(is.null(hemibrain.nblast) & repository == "CATMAID"){
    hemibrain.nblast = t(hemibrain_nblast("hemibrain-fafb14"))
  }
  if(is.null(hemibrain.nblast) & repository == "flywire"){
    hemibrain.nblast = hemibrain_nblast("hemibrain-flywire")
  }
  # FlyWire or CATMAID?
  if(repository=="CATMAID"){
    id = "skid"
    brain = elmr::FAFB14.surf
    extra.neurons = tryCatch(flywire_neurons(), error = function(e){
      message(e)
      NULL
    })
  }else if(repository == "flywire"){
    id = "root_id"
    brain = elmr::FAFB14.surf
    if(is.null(ids)){
      ids = tryCatch(names(flywire_neurons()), error = function(e){
        message(e)
        NULL
      })
    }
    if(is.null(query)){
      query = tryCatch(flywire_neurons(), error = function(e){
        message(e)
        NULL
      })
    }
    extra.neurons = tryCatch(flywire_neurons(mirror=TRUE), error = function(e){
      message(e)
      NULL
    })
  }
  chosen.field = "bodyid"
  match.field = search.id = paste0("hemibrain",".match")
  quality.field = paste0("hemibrain","_match_quality")
  # Read the Google Sheet
  gs = hemibrain_match_sheet(selected_file = selected_file, sheet = repository)
  gs$user[is.na(gs$user)] = ""
  # Get hemibrain neurons
  if(missing(db)) {
    db=tryCatch(force(db), error=function(e) {
      warning("Unable to use `hemibrain_neurons()`. ",
           "You must load the hemibrain Google Team Drive")
      NULL
    })
  }
  # How much is done
  done = subset(gs, !is.na(gs[[match.field]]))
  message("Neuron matches: ", nrow(done), "/", nrow(gs))
  print(table(gs[[quality.field]]))
  # Choose user
  initials = choose_user(gs, user = user)
  # choose ids
  selected = id_selector(gs=gs, ids=ids, id=id, overwrite = overwrite,
                         quality.field = quality.field, match.field = match.field,
                         initials = initials, column = column, entry = entry, superuser = superuser,
                         flywire.good = FALSE)
  # Make matches!
  match.more = TRUE
  while(match.more){
    match_cycle = neuron_match_scanner(brain = brain,
                                       selected = selected,
                                       id = id,
                                       unsaved = unsaved,
                                       saved = saved,
                                       chosen.field = chosen.field,
                                       nblast = hemibrain.nblast,
                                       threshold = threshold,
                                       batch_size = batch_size,
                                       targets = db,
                                       targets.repository = "hemibrain",
                                       query = query,
                                       query.repository = repository,
                                       extra.neurons = extra.neurons,
                                       extra.repository = setdiff(c("flywire","CATMAID"),repository),
                                       match.field = match.field,
                                       quality.field = quality.field,
                                       soma.size = 4000,
                                       show.columns = c("cell_type","ito_lee_hemilineage","status", match.field, quality.field,"note"))
    selected = match_cycle[["selected"]]
    unsaved = match_cycle[["unsaved"]]
    selected$user = initials
    if(length(unsaved)){
      plot_inspirobot()
      say_encouragement(initials)
      # Read!
      gs2 = hemibrain_match_sheet(selected_file = selected_file, sheet = repository)
      selected.unsaved = subset(selected, selected[[id]]%in%unsaved)
      for(i in unsaved){
        gs2[gs2[[id]]%in%i,c(match.field,quality.field,"note","user")] = selected.unsaved[match(i,selected.unsaved[[id]]),c(match.field,quality.field,"note","user")]
      }
      # Write!
      if(!identical(gs,gs2)){
        gsheet_update_cols(
          write.cols = c(match.field,quality.field,"note","user"),
          gs=gs2,
          selected_sheet = selected_file,
          sheet = "FAFB",
          Verbose = TRUE)
      }
      saved = c(unsaved, saved)
      unsaved = c()
    }
    match.more = hemibrain_choice("Match more neurons? y/n ")
  }
  say_encouragement(initials)
}

# neuron ID name
get_idfield <- function(repository = c("hemibrain","CATMAID","flywire","LM","lm","FAFB.hemisphere","FAFB"),
                        return = c("id.field","sheet","match.field","match_quality")){
  return = match.arg(return)
  repository = match.arg(repository)
  sheet = repository
  if(repository=="hemibrain"){
    id.field = "bodyid"
  }else if (repository == c("flywire")){
    id.field = "root_id"
    sheet = "FAFB"
  }else if (repository%in%c("CATMAID","FAFB")){
    id.field = "skid"
    sheet = "FAFB"
  }else if (repository=="FAFB"){
    id.field = "skid"
  }else if(repository=="FAFB.hemisphere"){
    id.field = "fafb_hemisphere_match"
  }else{
    id.field = "id"
  }
  if(return=="sheet"){
    sheet
  }else if(return%in%c("match.field","match_quality")){
    if(repository%in%c("CATMAID","FAFB")){
      match.field  = "fafb_match"
      match_quality = "fafb_match_quality"
    }else if(repository%in%c("LM","lm")){
      match.field = "lm_match"
      match_quality = "lm_match_quality"
    }else if(repository=="hemibrain"){
      match.field = "hemibrain_match"
      match_quality = "hemibrain_match_quality"
    }else if(repository=="flywire"){
      match.field = "flywire_xyz"
      match_quality = "fafb_match_quality"
    }else if(repository%in%c("FAFB.hemisphere","hemisphere")){
      match.field = "fafb_hemisphere_match"
      match_quality = "fafb_hemisphere_match_quality"
    }else{
      match.field = "match"
      match_quality = "quality"
    }
    if(return=="match.field"){
      match.field
    }else{
      match_quality
    }
  }else{
    id.field
  }
}

# flag for user
# hidden
matching_user <- function(selected_file,
                          repository,
                          ids,
                          user){
  if(length(user)>1){
    stop("Specify just one user")
  }
  id.field = get_idfield(repository, return = "id.field")
  ws = get_idfield(repository, return = "sheet")
  gs = hemibrain_match_sheet(selected_file = selected_file, sheet = repository)
  gs[gs[[id.field]]%in%ids,"user"] = user
  write_matches(gs=gs,
                ids = ids,
                id.field = id.field,
                column = "user",
                selected_file = selected_file,
                ws = ws)
}


#' Retrieve matched up neurons between the hemibrain and FAFB
#'
#' @description Many neurons in the hemibrain data have been matched one to one with FAFB neurons (hemibrain->FAFB) and
#' FAFB neurons have been matched with their hemibrain counterparts (FAFB->hemibrain). These matches have been done by the
#' Flyconnectome Group at the University of Cambridge, and are recorded on a
#' \href{https://docs.google.com/spreadsheets/d/1OSlDtnR3B1LiB5cwI5x5Ql6LkZd8JOS5bBr-HTi0pOw/edit#gid=0}{Google
#'   Sheet} on our hemibrain Google Team Drive, to which you will need access
#' through an authenticated account to view and use this function. Matches have three levels of 'quality', largely
#' dependent on the degree of manual tracing for FAFB neurons - good (could be the same cell), medium (same cell_type) and poor (could be the same or similar cell_type).
#'
#' @param flytable logical. If \code{TRUE} then match data is pulled from the Drosophila Connectomics Group's seatable databse.Otherwise, we will use a googlesheet.
#' indicated by \code{selected_file}.
#' @param priority whether to use FAFB->hemibrain_matches (FAFB) or hemibrain->FAFB matches (hemibrain) in order to ascribe
#' cell_type names to FAFB neurons. In both cases, cell_type names are attached to hemibrain bodyids, and propagated to their FAFB matches.
#' @param selected_file the Google Sheet database to read and write from. For
#'   now, defaults to this
#'   \href{https://docs.google.com/spreadsheets/d/1OSlDtnR3B1LiB5cwI5x5Ql6LkZd8JOS5bBr-HTi0pOw/edit#gid=0}{Google
#'    Sheet}.
#'
#' @return  a \code{data.frame} which includes neuron's ID (either its CATMAID skeleton ID or neuprint body ID), the data set from which it comes,
#' its putative cell_type and connectivity type, and its match in the other dataset.
#'
#' @details Currently, the
#'   \href{https://docs.google.com/spreadsheets/d/1OSlDtnR3B1LiB5cwI5x5Ql6LkZd8JOS5bBr-HTi0pOw/edit#gid=0}{Google
#'    Sheet} is set up with  limited number of users, each of whom have been
#'   assigned a number of neurons to match up. In order to add yourself as a
#'   user, simply open this Google Sheet in your browser and add your initials
#'   to neurons of your choosing on the rightmost column 'users'.
#'
#' @return a \code{data.frame} where each row is a neuron, either from the hemibrain or FAFB data sets. Each row gives you its matching neuron in the other data set. These matches have been
#' manually assigned using \code{\link{fafb_matching}}, \code{\link{hemibrain_matching}} and \code{\link{LR_matching}}. If you use this information make sure you credit it appropriately.
#' Contact us if unsure:
#' \itemize{
#'
#'   \item{"cell_type"}{ - the neuPrint designated 'type' for the neuron. If \code{dataset} is not \code{"hemibrain"}, then this is based on the hemibrain \code{match}.}
#'
#'   \item{"cell"}{ - the unique cell, which is just \code{cell_type#number}.}
#'
#'   \item{"cellBodyFiber"}{ - the cell body fiber to which this neuron belongs}
#'
#'   \item{"ito_lee_hemilineage"}{ - the hemilineage to which this neuron belongs. Seer \code{\link{hemibrain_hemilineages}}.}
#'
#'   \item{"match"}{ - the ID of the manual match from the other data set. If \code{dataset=="hemibrain"} then this is a \code{root_id} that can be found in \code{flywire_neurons}.If \code{"CATMAID"} or \code{"flywire"} then it is a hemibrain body ID.}
#'
#'   \item{"quality"}{ - the matcher makers qualitative assessment of how good this match is.}
#'
#'   \item{"fafb_hemisphere_match"}{ - the flywire coordinates of a neuron on the opposite hemisphere, which a match maker has designated as this \code{id}'s cognate.}
#'
#'   \item{"fafb_hemisphere_match_quality"}{ - the quality of this match.}
#'
#'   \item{"lm_match"}{ - indicates a light level neuron that is a match for \code{id}. This neuron will be in \code{flycircuit_neurons()} or other light level data.}
#'
#'   \item{"lm_match_quality"}{ - the quality of this match.}
#'
#'   \item{"dataset"}{ - the data set to which \code{id} belongs. Either \code{"hemibrain"}, or one of the two FAFB repositories, \code{"CATMAID"} or \code{"flywire"}.}
#'
#'   \item{"priority"}{ - whether FAFB->hemibrain_matches (\code{"FAFB"}) or hemibrain->FAFB matches (\code{"hemibrain"}) were used in order to ascribe cell_type names to FAFB neurons.
#'   In both cases, cell_type names are attached to hemibrain bodyids, and propagated to their FAFB matches.}
#'
#' }
#' @examples
#' \donttest{
#' \dontrun{
#'
#' # Get matches
#' matched = hemibrain_matches()
#'
#' }}
#' @seealso \code{\link{hemibrain_matching}},\code{\link{fafb_matching}},\code{\link{LR_matching}}
#' @importFrom stats ave
#' @rdname hemibrain_matches
#' @export
hemibrain_matches <- function(flytable = TRUE,
                              priority = c("FAFB","hemibrain"),
                              selected_file = options()$hemibrainr_matching_gsheet){
  priority = match.arg(priority)
  if (flytable){
    matched = flytable_matches()
  }else{
    # Get matches
    hemibrain_matches = gsheet_manipulation(FUN = googlesheets4::read_sheet,
                                            ss = selected_file,
                                            sheet = "hemibrain",
                                            return = TRUE,
                                            Verbose = FALSE)
    hemibrain_matches$bodyid = correct_id(hemibrain_matches$bodyid)
    hemibrain_matches = hemibrain_matches[!duplicated(hemibrain_matches$bodyid),]
    hemibrain_matches = hemibrain_matches[hemibrain_matches$bodyid!="",]
    hemibrain_matches$ito_lee_lineage = gsub("_.*","",hemibrain_matches$ito_lee_hemilineage)
    hemibrain_matches$dataset = "hemibrain"
    hemibrain_matches = subset(hemibrain_matches, !is.na(hemibrain_matches$bodyid))
    hemibrain_matches = hemibrain_matches[!duplicated(hemibrain_matches$bodyid),]
    rownames(hemibrain_matches) = hemibrain_matches$bodyid

    # Get FAFB matches
    fafb_matches = gsheet_manipulation(FUN = googlesheets4::read_sheet,
                                       ss = selected_file,
                                       sheet = "FAFB",
                                       return = TRUE,
                                       Verbose = FALSE)
    fafb_matches = subset(fafb_matches, !is.na(fafb_matches$skid)|!is.na(fafb_matches$root_id))
    fafb_matches = fafb_matches[!(duplicated(fafb_matches$skid)&duplicated(fafb_matches$root_id)),]
    fafb_matches$skid = correct_id(fafb_matches$skid)
    fafb_matches$root_id = correct_id(fafb_matches$root_id)
    fafb_matches$dataset = "FAFB"
    fafb_matches$skid[is.na(fafb_matches$skid)] = "missing"
    rownames(fafb_matches) = paste0(fafb_matches$skid,"#",ave(fafb_matches$skid,fafb_matches$skid,FUN= seq.int))

    # Set unassigned matches to top.nblast
    nblast.tops = c("CATMAID.nblast.top","flywire.nblast.top")
    for(nt in nblast.tops){
      if(!is.null(hemibrain_matches[[nt]])){
        present = !is.na(hemibrain_matches[[nt]])
        missing = is.na(hemibrain_matches$flywire_xyz)
        replace = (present+missing)>1
        if(sum(replace)){
          if(nt=="CATMAID.nblast.top"){
            hemibrain_matches$fafb_match[replace] = hemibrain_matches[[nt]][replace]
            hemibrain_matches$flywire_xyz[replace] = fafb_matches$flywire_xyz[match(hemibrain_matches$fafb_match[replace],fafb_matches$skid)]
            hemibrain_matches$root_id[replace] = fafb_matches$root_id[match(hemibrain_matches$fafb_match[replace],fafb_matches$skid)]
          }else{
            hemibrain_matches$root_id[replace] = hemibrain_matches[[nt]][replace]
            hemibrain_matches$flywire_xyz[replace] = fafb_matches$flywire_xyz[match(fafb_matches$root_id[replace],fafb_matches$root_id)]
            hemibrain_matches$fafb_match[replace] = fafb_matches$skid[match(hemibrain_matches$root_id[replace],fafb_matches$root_id)]
          }
          hemibrain_matches$fafb_match_quality[replace] = "NBLAST"
        }
      }
    }

    # Set unassigned matches to top.nblast
    nblast.tops = c("hemibrain.nblast.top")
    for(nt in nblast.tops){
      if(!is.null(fafb_matches[[nt]])){
        present = !is.na(fafb_matches[[nt]])
        missing = is.na(fafb_matches$hemibrain_match)
        replace = (present+missing)>1
        if(sum(replace)){
          fafb_matches$hemibrain_match[replace] = fafb_matches[[nt]][replace]
          fafb_matches$hemibrain_match_quality[replace] = "NBLAST"
        }
      }
    }

    # Match to other side
    side.match = fafb_matches$fafb_hemisphere_match[match(hemibrain_matches$flywire_xyz,fafb_matches$flywire_xyz)]
    side.quality = fafb_matches$fafb_hemisphere_match_quality[match(hemibrain_matches$flywire_xyz,fafb_matches$flywire_xyz)]
    hemibrain_matches$fafb_hemisphere_match = ifelse(length(side.match),side.match,NA)
    hemibrain_matches$fafb_hemisphere_match_quality = ifelse(length(side.quality),side.quality,NA)

    # Unmatched
    fafb.unmatched = is.na(fafb_matches$hemibrain_match)|fafb_matches$hemibrain_match==""
    fafb_matches$hemibrain_match_quality[fafb.unmatched] = "none"
    fafb_matches$hemibrain_match[fafb.unmatched] = "none"
    hemibrain.unmatched = is.na(hemibrain_matches$fafb_match)|hemibrain_matches$fafb_match==""
    hemibrain_matches$fafb_match_quality[hemibrain.unmatched] = "none"
    hemibrain_matches$fafb_match[hemibrain.unmatched] = "none"
    side.unmatched = is.na(hemibrain_matches$fafb_hemisphere_match)|hemibrain_matches$fafb_hemisphere_match==""
    hemibrain_matches$fafb_hemisphere_match_quality[side.unmatched] = "none"
    hemibrain_matches$fafb_hemisphere_match[side.unmatched] = "none"
    side.unmatched = is.na(fafb_matches$fafb_hemisphere_match)|fafb_matches$fafb_hemisphere_match==""
    fafb_matches$fafb_hemisphere_match_quality[side.unmatched] = "none"
    fafb_matches$fafb_hemisphere_match[side.unmatched] = "none"

    # Add in neuprint types
    meta = neuprintr::neuprint_get_meta(hemibrain_matches$bodyid)
    types = meta$type
    names(types) = meta$bodyid
    hemibrain_matches[,"cell_type"] = types[as.character(hemibrain_matches$bodyid)]

    # Address FAFB cell_types
    fafb_matches$cell_type = NA
    for(q in c("good","medium","poor")){
      if(is.na(q)){
        next
      }
      if(priority=="FAFB"){
        # FAFB -> Hemibrain
        inq = as.character(subset(fafb_matches, fafb_matches$hemibrain_match_quality==q)$flywire_xyz)
        inq = inq[!is.na(inq)]
        isna = is.na(fafb_matches[inq,]$cell_type)
        cts = hemibrain_matches$cell_type[match(fafb_matches[inq,]$hemibrain_match,hemibrain_matches$bodyid)]
        fafb_matches[inq,]$cell_type[isna] = cts[isna]
        # Hemibrain -> FAFB
        inq = as.character(subset(hemibrain_matches, hemibrain_matches$fafb_match_quality==q)$flywire_xyz)
        inq = inq[!is.na(inq)]
        isna = is.na(fafb_matches[inq,]$cell_type)
        fafb_matches[inq,]$cell_type[isna] = hemibrain_matches$cell_type[match(fafb_matches[inq,]$flywire_xyz[isna],hemibrain_matches$flywire_xyz)]
      }else{
        # Hemibrain -> FAFB
        inq = as.character(subset(hemibrain_matches, hemibrain_matches$fafb_match_quality==q)$flywire_xyz)
        inq = inq[!is.na(inq)]
        isna = is.na(fafb_matches[inq,]$cell_type)
        fafb_matches[inq,]$cell_type[isna] = hemibrain_matches$cell_type[match(fafb_matches[inq,]$flywire_xyz[isna],hemibrain_matches$flywire_xyz)]
        # FAFB -> Hemibrain
        inq = as.character(subset(fafb_matches, fafb_matches$hemibrain_match_quality==q)$flywire_xyz)
        inq = inq[!is.na(inq)]
        isna = is.na(fafb_matches[inq,]$cell_type)
        cts = hemibrain_matches$cell_type[match(fafb_matches[inq,]$hemibrain_match,hemibrain_matches$bodyid)]
        fafb_matches[inq,]$cell_type[isna] = cts[isna]
      }
    }

    # Add in neuprint types
    ntotype = fafb_matches$hemibrain_match[is.na(fafb_matches$cell_type)]
    ntotype = ntotype[ntotype!="none" & !is.na(ntotype)]
    meta = neuprintr::neuprint_get_meta(ntotype)
    types = meta$type
    names(types) = meta$bodyid
    types = types[!is.na(types)]
    fafb_matches[match(names(types),fafb_matches$hemibrain_match),"cell_type"] = types

    # Work out lineages
    for(repo in c("CATMAID","flywire")){
      id.field = get_idfield(repository = repo, return = "match.field")
      match.field = get_idfield(repository = repo, return = "match.field")
      ids = id_okay(fafb_matches[[id.field]])
      for(id in ids){
        idx = match(id,fafb_matches[[id.field]])
        idx2 = match(id,hemibrain_matches[[match.field]])
        if(!is.na(idx2)){
          ct = fafb_matches[idx,"cell_type"]
          if(is.na(ct)){
            fafb_matches[idx,"cell_type"] = hemibrain_matches$cell_type[idx2]
            if(is.na(fafb_matches[idx,"cell_type"])){
              fafb_matches[idx,"cell_type"] = "uncertain"
            }else{
              fafb_matches[idx,"ito_lee_hemilineage"] = hemibrain_matches$ito_lee_hemilineage[idx2]
            }
          }else{
            fafb_matches[idx,"ito_lee_hemilineage"] = hemibrain_matches$ito_lee_hemilineage[match(ct,hemibrain_matches$cell_type)]
          }
        }
      }
    }

    # Add cell body fiber infor for FAFB cells
    cbf.match = hemibrain_matches$cellBodyFiber[match(fafb_matches$hemibrain_match,hemibrain_matches$bodyid)]
    fafb_matches$cellBodyFiber = ifelse(length(cbf.match),cbf.match,NA)

    # Rename cells
    fafb_matches = fafb_matches[order(fafb_matches$hemibrain_match,decreasing = TRUE),]
    hemibrain_matches = hemibrain_matches[order(hemibrain_matches$bodyid,decreasing = TRUE),]
    fafb_matches$cell = paste0(fafb_matches$cell_type,"#",ave(fafb_matches$cell_type,fafb_matches$cell_type,FUN= seq.int))
    hemibrain_matches$cell = paste0(hemibrain_matches$cell_type,"#",ave(hemibrain_matches$cell_type,hemibrain_matches$cell_type,FUN= seq.int))

    # Separate flywire and CATMAID matches
    catmaid.matches = subset(fafb_matches, !is.na(fafb_matches$skid)& !fafb_matches$skid%in%c(""," ","NA","0","error","none"))
    flywire.matches = subset(fafb_matches, !is.na(fafb_matches$root_id)& !fafb_matches$root_id%in%c(""," ","NA","0","error","none"))
    hemibrain_matches$dataset = "hemibrain"
    flywire.matches$dataset = "flywire"
    catmaid.matches$dataset = "CATMAID"
    hemibrain_matches$match.dataset = "flywire"
    flywire.matches$match.dataset = "hemibrain"
    catmaid.matches$match.dataset = "hemibrain"
    hemibrain_matches.catmaid = hemibrain_matches
    hemibrain_matches.catmaid$match.dataset = "CATMAID"

    # Make matching data frame
    matched.h = hemibrain_matches[,c("bodyid", "cell_type", "cell", "cellBodyFiber", "ito_lee_hemilineage",
                                     "root_id", "fafb_match_quality",
                                     "fafb_hemisphere_match", "fafb_hemisphere_match_quality",
                                     "lm_match", "lm_match_quality", "dataset", "match.dataset")]
    matched.hc = hemibrain_matches.catmaid[,c("bodyid", "cell_type", "cell", "cellBodyFiber", "ito_lee_hemilineage",
                                              "fafb_match", "fafb_match_quality",
                                              "fafb_hemisphere_match", "fafb_hemisphere_match_quality",
                                              "lm_match", "lm_match_quality", "dataset", "match.dataset")]
    matched.c = catmaid.matches[,c("skid",  "cell_type",  "cell", "cellBodyFiber", "ito_lee_hemilineage",
                                   "hemibrain_match", "hemibrain_match_quality",
                                   "fafb_hemisphere_match", "fafb_hemisphere_match_quality",
                                   "lm_match", "lm_match_quality","dataset", "match.dataset")]
    matched.f = flywire.matches[,c("root_id",  "cell_type",  "cell", "cellBodyFiber", "ito_lee_hemilineage",
                                   "hemibrain_match", "hemibrain_match_quality",
                                   "fafb_hemisphere_match", "fafb_hemisphere_match_quality",
                                   "lm_match", "lm_match_quality","dataset", "match.dataset")]
    colnames(matched.h) = colnames(matched.hc) = colnames(matched.c) = colnames(matched.f) = c("id","cell_type", "cell","cellBodyFiber","ito_lee_hemilineage",
                                                                                               "match","quality",
                                                                                               "fafb_hemisphere_match", "fafb_hemisphere_match_quality",
                                                                                               "lm_match", "lm_match_quality","dataset", "match.dataset")
    matched = rbind(matched.h,matched.hc,matched.f,matched.c)
    matched$quality[is.na(matched$match)] = "none"
    matched$match[is.na(matched$match)] = "none"
    matched$ito_lee_lineage = gsub("_.*","",matched$ito_lee_hemilineage)

    # Sort out types
    matched$connectivity.type = matched$cell_type
    matched$cell_type = gsub("_[a-z]{1}$","",matched$cell_type)
    matched$cell_type[is.na(matched$cell_type)] = "uncertain"
    matched$connectivity.type[is.na(matched$connectivity.type)] = "uncertain"
    matched$priority = priority
    rownames(matched) = 1:nrow(matched)


  }

  # Return
  matched
}

#' @rdname hemibrain_matches
#' @export
lm_matches <- function(priority = c("hemibrain","lm"), selected_file = options()$hemibrainr_matching_gsheet){
  priority = match.arg(priority)

  # Get matches
  hemibrain_matches = gsheet_manipulation(FUN = googlesheets4::read_sheet,
                                          ss = selected_file,
                                          sheet = "hemibrain",
                                          return = TRUE,
                                          Verbose = FALSE)
  hemibrain_matches$bodyid = correct_id(hemibrain_matches$bodyid)
  hemibrain_matches = hemibrain_matches[!duplicated(hemibrain_matches$bodyid),]
  hemibrain_matches = hemibrain_matches[hemibrain_matches$bodyid!="",]
  hemibrain_matches$ito_lee_lineage = gsub("_.*","",hemibrain_matches$ito_lee_hemilineage)
  hemibrain_matches$dataset = "hemibrain"
  hemibrain_matches = subset(hemibrain_matches, !is.na(hemibrain_matches$bodyid))
  hemibrain_matches = hemibrain_matches[!duplicated(hemibrain_matches$bodyid),]
  rownames(hemibrain_matches) = hemibrain_matches$bodyid
  hemibrain_matches$cell_type = hemibrain_matches$connectivity.type

  # Get lm matches
  lm_matches = gsheet_manipulation(FUN = googlesheets4::read_sheet,
                                     ss = selected_file,
                                     sheet = "lm",
                                     return = TRUE,
                                   Verbose = FALSE)
  lm_matches$id = correct_id(lm_matches$id)
  lm_matches = lm_matches[!duplicated(lm_matches$id),]
  lm_matches = lm_matches[lm_matches$id!="",]
  lm_matches = subset(lm_matches, !is.na(lm_matches$id))
  lm_matches = lm_matches[!duplicated(lm_matches$id),]
  lm_matches$dataset = "lm"
  rownames(lm_matches) = lm_matches$id

  # Unmatched
  unmatched = is.na(lm_matches$hemibrain_match)|lm_matches$hemibrain_match==""
  lm_matches$hemibrain_match_quality[unmatched] = "none"
  lm_matches$hemibrain_match[unmatched] = "none"

  # Add in neuprint types
  meta = neuprintr::neuprint_get_meta(hemibrain_matches$bodyid)
  types = meta$type
  names(types) = meta$bodyid
  hemibrain_matches[,"cell_type"] = types[as.character(hemibrain_matches$bodyid)]

  # Address lm cell_types
  lm_matches$cell_type = NA
  for(q in c("good","medium","poor")){
    if(is.na(q)){
      next
    }
    if(priority=="lm"){
      # lm -> Hemibrain
      inq = as.character(subset(lm_matches, lm_matches$hemibrain_match_quality==q)$id)
      inq = inq[!is.na(inq)]
      isna = is.na(lm_matches[inq,]$cell_type)
      cts = hemibrain_matches$cell_type[match(lm_matches[inq,]$hemibrain_match,hemibrain_matches$bodyid)]
      lm_matches[inq,]$cell_type[isna] = cts[isna]
      # Hemibrain -> lm
      inq = as.character(subset(hemibrain_matches, hemibrain_matches$lm_match_quality==q)$lm_match)
      inq = inq[!is.na(inq)]
      isna = is.na(lm_matches[inq,]$cell_type)
      lm_matches[inq,]$cell_type[isna] = hemibrain_matches$cell_type[match(lm_matches[inq,]$id[isna],hemibrain_matches$lm_match)]
    }else{
      # Hemibrain -> lm
      inq = as.character(subset(hemibrain_matches, hemibrain_matches$lm_match_quality==q)$lm_match)
      inq = inq[!is.na(inq)]
      isna = is.na(lm_matches[inq,]$cell_type)
      lm_matches[inq,]$cell_type[isna] = hemibrain_matches$cell_type[match(lm_matches[inq,]$id[isna],hemibrain_matches$lm_match)]
      # lm -> Hemibrain
      inq = as.character(subset(lm_matches, lm_matches$hemibrain_match_quality==q)$id)
      inq = inq[!is.na(inq)]
      isna = is.na(lm_matches[inq,]$cell_type)
      cts = hemibrain_matches$cell_type[match(lm_matches[inq,]$hemibrain_match,hemibrain_matches$bodyid)]
      lm_matches[inq,]$cell_type[isna] = cts[isna]
    }
  }

  # Add in neuprint types
  ntotype = lm_matches$hemibrain_match[is.na(lm_matches$cell_type)]
  ntotype = ntotype[ntotype!="none" & !is.na(ntotype)]
  meta = neuprintr::neuprint_get_meta(ntotype)
  types = meta$type
  names(types) = meta$bodyid
  types = types[!is.na(types)]
  lm_matches[match(names(types),lm_matches$hemibrain_match),"cell_type"] = types

  # Work out lineages
  for(id in as.character(lm_matches$id)){
    if(is.na(id)){
      next
    }
    ct = lm_matches[id,"cell_type"]
    if(is.na(ct)){
      lm_matches[id,"cell_type"] = hemibrain_matches$cell_type[match(id,hemibrain_matches$lm_match)]
      lm_matches[id,"ito_lee_hemilineage"] = hemibrain_matches$ito_lee_hemilineage[match(id,hemibrain_matches$lm_match)]
    }else{
      lm_matches[id,"ito_lee_hemilineage"] = hemibrain_matches$ito_lee_hemilineage[match(ct,hemibrain_matches$cell_type)]
    }
  }

  # Rename cells
  lm_matches$cell = paste0(lm_matches$cell_type,"#",ave(lm_matches$cell_type,lm_matches$cell_type,FUN= seq.int))
  hemibrain_matches$cell = paste0(hemibrain_matches$cell_type,"#",ave(hemibrain_matches$cell_type,hemibrain_matches$cell_type,FUN= seq.int))

  # Fix hemilineages
  hl = hemibrain_matches$ito_lee_hemilineage[match(lm_matches$id,hemibrain_matches$lm_match)]
  l = hemibrain_matches$ito_lee_lineage[match(lm_matches$id,hemibrain_matches$lm_match)]
  lm_matches$ito_lee_hemilineage[!is.na(hl)] = hl[!is.na(hl)]
  lm_matches$ito_lee_lineage[!is.na(l)] = l[!is.na(l)]

  # Add cell body fiber infor for FAFB cells
  lm_matches$cellBodyFiber = hemibrain_matches$cellBodyFiber[match(lm_matches$hemibrain_match,hemibrain_matches$bodyid)]

  # Make matching data frame
  hemibrain_matches$dataset = "hemibrain"
  lm_matches$dataset = "lm"
  matched.h = hemibrain_matches[,c("bodyid", "cell_type", "cell", "cellBodyFiber", "ito_lee_hemilineage",
                                   "lm_match", "lm_match_quality", "fafb_match", "fafb_match_quality", "dataset")]
  matched.f = lm_matches[,c("id",  "cell_type",  "cell", "cellBodyFiber", "ito_lee_hemilineage",
                              "hemibrain_match", "hemibrain_match_quality", "fafb_match", "fafb_match_quality","dataset")]
  colnames(matched.h) = colnames(matched.f) = c("id","cell_type", "cell","cellBodyFiber","ito_lee_hemilineage","match","quality", "fafb_match", "fafb_match_quality","dataset")
  matched = rbind(matched.h,matched.f)
  matched$quality[is.na(matched$match)] = "none"
  matched$match[is.na(matched$match)] = "none"

  # Sort out types
  matched$connectivity.type = matched$cell_type
  matched$cell_type = gsub("_.*","",matched$cell_type)
  matched$cell_type[is.na(matched$cell_type)] = "uncertain"
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
#'  - good (could be the same cell), medium (same cell_type) and poor (could be the same or similar cell_type).
#' @param ids either hemibrain bodyids or FAFB skids to add to a Google Sheet. You will want to add IDs where they do not already exist, so that
#'  you can use, for example, \code{\link{hemibrain_matching}}.
#' @param dataset the tab to which to add your new information. You are either adding to information related to hemibrain neurons, or FAFB neurons.
#' @param direction the match direction, i.e. hemibrain->FAFB (hemibrain tab) or FAFB->hemibrain (fafb tab). Defaults to updating both.
#' @param selected_file Specifies which Google Sheet to use. Unless you are using a personal Google Sheet, this should be \code{options()$hemibrainr_matching_gsheet}.
#' @param user your initials, so updates can be linked to you. There is a limited number of users, each of whom have been
#'   assigned a number of neurons to match up. In order to add yourself as a user,
#'   simply open this Google Sheet in your browser and add your initials to neurons of your choosing on the rightmost column 'users'.
#' @param flywire_ids flywire IDs to add to Google sheet if not already present.
#' @param flycircuit.ids flycircuit IDs to add to Google sheet if not already present.
#' @param meta meta data for the given flycircuit IDs.
#' @param top.nblast logical. Whether or not to also give the top NBLAST match for each entry.
#' @param catmaid.update logical. Whether or not to update \code{flywire_xyz} and \code{root_id} columns, based on e CATMAID neuron specified by a \code{skid} column.
#' @param reorder logical. Whether or not to re-write the sheet so that it is ordered by hemilineage.
#' @param nblast if \code{top.nblast} is \code{TRUE} this nblast matrix is used to update the column \code{top.nblast}. If set to \code{NULL} defaults to using \code{hemibrain_nblast}. Columns should be hemibrain neurons, and rows the other data set.
#' @param nblast.hemibrain.catmaid if \code{top.nblast} is \code{TRUE} this nblast matrix is used to update the column \code{top.nblast}. If set to \code{NULL} defaults to using \code{hemibrain_nblast}. Columns should be hemibrain neurons, and rows CATMAID neurons.
#' @param nblast.hemibrain.flywire if \code{top.nblast} is \code{TRUE} this nblast matrix is used to update the column \code{top.nblast}. If set to \code{NULL} defaults to using \code{hemibrain_nblast}. Columns should be hemibrain neurons, and rows flywire neurons.
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
#' hemibrain_matching_add(ids = "16", dataset = "FAFB", user = "ASB")
#'
#'
#' }}
#' @seealso \code{\link{hemibrain_matching}}, \code{\link{fafb_matching}}
#' @rdname hemibrain_add_made_matches
hemibrain_add_made_matches <- function(df,
                                  direction = c("both","hemibrain-FAFB","FAFB-hemibrain","hemisphere"),
                                  user = "flyconnectome",
                                  selected_file = options()$hemibrainr_matching_gsheet,
                                  ...){
  direction = match.arg(direction)
  if(direction %in% c("both","hemibrain-FAFB","FAFB-hemibrain")){
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
        hemibrain_matching_add(ids = missing, dataset = "hemibrain", user = user, selected_file = selected_file)
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
      gs[as.character(hdf$bodyid),"fafb_match"] = hdf$skid
      gs[as.character(hdf$bodyid),"fafb_match_quality"] = hdf$quality
      gs[as.character(hdf$bodyid),"user"] = user
      message("Adding matches")
      write_matches(gs=gs,
                    ids = as.character(hdf$bodyid),
                    ws="hemibrain",
                    id.field ="bodyid",
                    selected_file = selected_file,
                    column = "fafb_match")
      write_matches(gs=gs,
                    ids = as.character(hdf$bodyid),
                    ws="hemibrain",
                    id.field ="bodyid",
                    selected_file = selected_file,
                    column = "fafb_match_quality")
    }
    if(direction%in%c("both","FAFB-hemibrain")){
      message("Reading FAFB sheet")
      gs = hemibrain_match_sheet(sheet = "FAFB")
      missing = setdiff(df$skid,gs$skid)
      if(length(missing)){
        message("Adding missing FAFB bodyids")
        hemibrain_matching_add(ids = missing, dataset = "FAFB", user = user, selected_file = selected_file)
        gs = hemibrain_match_sheet(sheet = "FAFB")
      }
      message("Checking that hemibrain_matches exist")
      hdf = subset(df, df$skid %in% gs$skid)
      hdf = hdf[!duplicated(hdf$skid),]
      in.hemi = suppressWarnings(sapply(hdf$skid,function(x) tryCatch(neuprintr::neuprint_ids(x),error = function(e) NULL)))
      in.hemi = sapply(in.hemi,function(x) length(x)!=0)
      if(sum(!in.hemi)){
        message("Dropping ", sum(!in.hemi)," matches as givens skeleton IDs do not exist for FAFB v14")
        hdf = hdf[in.hemi,]
      }
      gs[as.character(hdf$skid),"hemibrain_match"] = hdf$bodyid
      gs[as.character(hdf$skid),"hemibrain_match_quality"] = hdf$quality
      gs[as.character(hdf$skid),"user"] = user
      message("Adding matches")
      write_matches(gs=gs,
                    ids = as.character(hdf$skid),
                    ws="FAFB",
                    id.field ="skid",
                    selected_file = selected_file,
                    column = "hemibrain_match")
      write_matches(gs=gs,
                    ids = as.character(hdf$skid),
                    ws="FAFB",
                    id.field ="skid",
                    selected_file = selected_file,
                    column = "hemibrain_match_quality")
    }
  }else{
    cnames = c("root_id","fafb_hemisphere_match","fafb_hemisphere_match_quality")
    if(!all(cnames%in%colnames(df))){
      stop("df must have column names: ", paste(cnames,collapse = ", "))
    }
    message("Adding matches")
    df = subset(df, !is.na(df$root_id) & !is.na(df$fafb_hemisphere_match) & !is.na(df$fafb_hemisphere_match_quality))
    df$root_id = correct_id(df$root_id)
    df$fafb_hemisphere_match_quality = correct_id(df$fafb_hemisphere_match_quality)
    df$root_id = unlist(sapply(df$root_id, function(x) tryCatch(fafbseg::flywire_latestid(x), error = function(e) x)))
    gs = hemibrain_match_sheet(sheet = "flywire")
    missing = setdiff(df$root_id,gs$root_id)
    if(length(missing)){
      message("Adding missing flywire neurons")
      meta = flywire_meta()
      hemibrain_matching_add(ids = missing,  dataset="flywire", meta = meta, selected_file = selected_file) # meta
      gs = hemibrain_match_sheet(sheet = "flywire")
    }
    gs.good = subset(gs, !is.na(gs$fafb_hemisphere_match_quality) & gs$fafb_hemisphere_match_quality!="none")
    df = subset(df, !df$root_id %in% gs.good$root_id)
    gs.new = gs
    gs.new[,"fafb_hemisphere_match"] = df[match(gs$root_id,df$root_id),"fafb_hemisphere_match"]
    gs.new[,"fafb_hemisphere_match_quality"] = df[match(gs$root_id,df$root_id),"fafb_hemisphere_match_quality"]
    gs.new.na = is.na(gs.new[,"fafb_hemisphere_match"])
    gs.old.filled = !is.na(gs[,"fafb_hemisphere_match_quality"])
    gs.new[gs.new.na,"fafb_hemisphere_match"] = gs[gs.new.na,"fafb_hemisphere_match"]
    gs.new[gs.new.na,"fafb_hemisphere_match_quality"] = gs[gs.new.na,"fafb_hemisphere_match_quality"]
    gs.new[gs.old.filled,"fafb_hemisphere_match"] = gs[gs.old.filled,"fafb_hemisphere_match"]
    gs.new[gs.old.filled,"fafb_hemisphere_match_quality"] = gs[gs.old.filled,"fafb_hemisphere_match_quality"]
    write_matches(gs=gs.new,
                  ids = as.character(df$root_id),
                  ws="FAFB",
                  id.field ="root_id",
                  selected_file = selected_file,
                  column = "fafb_hemisphere_match")
    write_matches(gs=gs.new,
                  ids = as.character(df$root_id),
                  ws="FAFB",
                  id.field ="root_id",
                  selected_file = selected_file,
                  column = "fafb_hemisphere_match_quality")
  }
}

# Get correct GSheet
hemibrain_match_sheet <- function(selected_file = options()$hemibrainr_matching_gsheet,
                                  sheet = c("hemibrain","FAFB","CATMAID","flywire","lm")){
  # Which sheet
  sheet = match.arg(sheet)

  # neuron ID name
  id.field = get_idfield(repository = sheet, return = "id.field")
  ws = get_idfield(repository = sheet, return = "sheet")

  # Read sheet
  gs = gsheet_manipulation(FUN = googlesheets4::read_sheet,
                                        ss = selected_file,
                                        sheet = ws,
                                        return = TRUE,
                           Verbose = FALSE)
  if(!id.field%in%colnames(gs)){
    stop(id.field," missing from column names of googlesheet: ", selected_file)
  }
  if(nrow(gs)){
    gs[[id.field]] = correct_id(gs[[id.field]])
    ids = gs[[id.field]]
    ids[is.na(ids)] = paste0("missing_",1:sum(is.na(ids)))
    ids = paste0(ids,"#",ave(ids,ids,FUN= seq.int))
    ids = gsub("#1$","",ids)
    rownames(gs) = ids
  }
  gs
}

#' @rdname hemibrain_add_made_matches
#' @export
hemibrain_matching_add <- function(ids = NULL,
                                   meta = NULL,
                                   dataset = c("hemibrain","FAFB","CATMAID","flywire","lm"),
                                   user = dataset,
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
    all.ids = unique(gs$flywire_ids)
    add = setdiff(fw.ids, all.ids)
  }
  add = setdiff(ids, rownames(gs))
  if(length(add)){
    # Meta information
    if(is.null(meta)){
      if(dataset=="hemibrain"){
        meta = hemibrain_get_meta(add, ...)
        meta$cell_type = meta$type
      }else if (dataset == "FAFB"){
        meta = elmr::fafb_get_meta(add, ...)
      } else if (dataset == "flywire"){
        add = setdiff(ids, gs$root_id)
        if(is.null(meta)){
          meta = flywire_neurons()[as.character(add),]
        }else{
          meta = meta[meta$root_id%in%add,]
        }
        meta = meta[!meta$flywire_xyz%in%gs$flywire_xyz,]
        if(!nrow(meta)){
          stop("Selected IDs could not be added. They must be among the neurons
         saved on Google drive, see flywire_neurons()")
        }
      }
    }else{
      id = if(dataset%in%c("FAFB","CATMAID")){
        "skid"
      }else if(dataset=="flywire"){
        "root_id"
      }else if(dataset =="hemibrain"){
        "bodyid"
      }else{
        "id"
      }
      add = setdiff(ids, gs[[id]])
      meta = meta[meta[[id]]%in%add,]
    }
    if(!nrow(meta)){
      warning("Given IDs already exist in sheet")
    }else{
      missing = setdiff(colnames(gs),colnames(meta))
      meta = add_blanks(meta, missing)
      meta = meta[,colnames(gs)]
      meta$user = user

      # Add new rows
      dataset[dataset=="flywire"] = "FAFB"
      batches = split(1:nrow(meta), ceiling(seq_along(1:nrow(meta))/500))
      for(i in batches){
        gsheet_manipulation(FUN = googlesheets4::sheet_append,
                            data = meta[min(i):max(i),],
                            ss = selected_file,
                            sheet = dataset)
      }
    }
  }else{
    warning("Given IDs already exist in sheet")
  }
}

#' @rdname hemibrain_add_made_matches
#' @export
hemibrain_matching_transfers <- function(selected_file = options()$hemibrainr_matching_gsheet){
  #############
  # Transfers #
  #############
  for(repo1 in c("CATMAID","flywire","hemibrain","FAFB.hemisphere")){
    for(repo2 in c("CATMAID","flywire","hemibrain","FAFB.hemisphere")){
      ws.1 = get_idfield(repo1, return = "sheet")
      ws.2 = get_idfield(repo2, return = "sheet")
      if(repo1==repo2|ws.1==ws.2){
        next
      }else{
        gs.1 = hemibrain_match_sheet(selected_file = selected_file, sheet = ws.1)
        gs.2 = hemibrain_match_sheet(selected_file = selected_file, sheet = ws.2)
        used.cols = colnames(gs.1)
        match.field.2 = get_idfield(repo2, return = "match.field")
        match_quality.2 = get_idfield(repo2, return = "match_quality")
        id.field.1 = get_idfield(repo1, return = "id.field")
        match.field.1 = get_idfield(repo1, return = "match.field")
        match_quality.1 = get_idfield(repo1, return = "match_quality")
        id.field.2 = get_idfield(repo2, return = "id.field")
        missing = is.na(gs.1[[match.field.2]])| (gs.1[[match.field.2]] %in% c("none"," ",""))
        matches = gs.2[[id.field.2]][match(gs.1[[id.field.1]][missing],gs.2[[match.field.1]])]
        if(length(matches)&sum(missing)){
          quality = gs.2[[match_quality.1]][match(gs.1[[id.field.1]][missing],gs.2[[match.field.1]])]
          matches[is.na(quality)] = NA
          quality[is.na(matches)] = NA
          gs.1[missing,match.field.2] = matches
          gs.1[missing,match_quality.2] = quality
          write.cols = intersect(c(match.field.2,match_quality.2),used.cols)
          message("Updating ",ws.1," with ", sum(!is.na(matches))," ",match.field.2,"es for column: ", id.field.1)
          if(length(quality)&sum(missing)){
            gsheet_update_cols(
              write.cols = write.cols,
              gs=gs.1[,used.cols],
              selected_sheet = selected_file,
              sheet = ws.1)
          }
        }
      }
    }
  }
}

# Udate function
update_gsheet <- function(update,
                          gs,
                          selected_file,
                          tab,
                          match = c("hemibrain", "lm", "FAFB", "flywire"),
                          id){
  match = match.arg(match)
  if(match=="flywire"){
    match2 = "FAFB"
  }else{
    match2 = match
  }
  for(row in 1:nrow(update)){
    columns = c(paste0(match,".match"), paste0(match2,"_match_quality"))
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
                          data = as.data.frame(update[row,column], stringsAsFactors = FALSE),
                          sheet = tab,
                          col_names = FALSE,
                          Verbose = TRUE)
    }
  }

}

#' @rdname hemibrain_add_made_matches
#' @export
fafb_matching_rewrite <- function(selected_file  = options()$hemibrainr_matching_gsheet,
                                  top.nblast = FALSE,
                                  reorder = FALSE,
                                  nblast = NULL,
                                   ...){
  n = hemibrain_match_sheet(sheet = "FAFB", selected_file = selected_file)
  n1 = elmr::fafb_get_meta("annotation:Lineage_annotated", batch = 10, ...)
  n2 = subset(n1, n1$skid %in% n$skid)
  n[match(n2$skid,n$skid),c("skid","ito_lee_hemilineage", "hartenstein_hemilineage", "cell_body_fiber")] = n2[,c("skid","ito_lee_hemilineage", "hartenstein_hemilineage", "cell_body_fiber")]
  ids.missing = as.character(setdiff(n1$skid,n$skid))
  ids.missing = id_okay(ids.missing)
  if(length(ids.missing)){
    n3 = elmr::fafb_get_meta(unique(ids.missing), batch = 10, ...)
    n = plyr::rbind.fill(n, n3[,c("skid","ito_lee_hemilineage", "hartenstein_hemilineage", "cell_body_fiber")])
    reorder = TRUE
  }
  matches = tryCatch(hemibrain_matches(selected_file=selected_file), error = function(e) NULL)
  matches = subset(matches, ! matches$quality %in% "NBLAST")
  if(!is.null(matches)){
    n$cell_type = matches[match(n$skid, matches$id),"connectivity.type"]
  }else{
    n$cell_type = NULL
  }
  lskids = as.character(catmaid::catmaid_skids("annotation:side: left", ...))
  n$side[is.na(n$side)] = "right"
  n[n$skid%in%lskids,"side"] = "left"
  if(top.nblast){
    if(is.null(nblast)){
      nblast = tryCatch(hemibrain_nblast('hemibrain-fafb14'), error = function(e) NULL)
    }
    if(!is.null(nblast)){
      nblast.top =nblast[match(n$skid,rownames(nblast)),]
      tops = apply(nblast.top,1,function(r) which.max(r))
      top = colnames(nblast)[unlist(tops)]
      top[!n$skid%in%rownames(nblast)] = NA
      n$nblast.top = top
    }
  }

  # Write to google sheet
  if(reorder){
    n$skid = correct_id(n$skid)
    gsheet_reorder(gs=n,tab="FAFB",selected_sheet=selected_file,field = "flywire_xyz", remove.duplicates = FALSE)
  }else{
    # Update
    write.cols = intersect(c("connectivity.type","cell_type","nblast.top","side",
                             "ito_lee_hemilineage","hartenstein_hemilineage","cell_body_fiber"),
                           colnames(n))
    gsheet_update_cols(
      write.cols = write.cols,
      gs=n,
      selected_sheet = selected_file,
      sheet = "FAFB",
      Verbose = FALSE)
  }

  # Add any missing data
  if(!is.null(matches)){
    missing = setdiff(subset(matches,matches$dataset=="hemibrain" & matches$match.dataset == "CATMAID")$match,
                      subset(matches,matches$dataset=="CATMAID" & matches$match.dataset == "hemibrain")$id)
    missing = setdiff(missing, n$skid)
    missing = missing[!grepl("missing|none|NA|good|medium|poor|tract|,|;|)",missing)]
    missing = id_okay(missing)
    if(length(missing)){
      tryCatch(hemibrain_matching_add(ids = missing, dataset="FAFB", selected_file = selected_file, ...),
               error = function(e) warning(e))
    }
  }
}

# hidden
id_okay <- function(x, zero = TRUE, logical = FALSE){
  x = x[!is.na(x)]
  x = x[!x%in%c(""," ","NA","none","missing","error")]
  if(zero){
    x = x[x!="0"]
    x = x[x!=0]
  }
  x = as.character(correct_id(x))
  x
}

#' @rdname hemibrain_add_made_matches
#' @export
hemibrain_matching_rewrite <- function(ids = NULL,
                                       selected_file  = options()$hemibrainr_matching_gsheet,
                                       top.nblast = FALSE,
                                       meta = NULL,
                                       nblast.hemibrain.catmaid = NULL,
                                       nblast.hemibrain.flywire = NULL,
                                  ...){
  gs = hemibrain_match_sheet(sheet = "hemibrain", selected_file = selected_file)
  if(is.null(ids)){
    ids = tryCatch(hemibrain_neuron_bodyids(), error  = function(e) unique(hemibrainr::hemibrain_metrics$bodyid))
  }
  meta1 = hemibrain_get_meta(unique(ids), ...)
  ids.missing = setdiff(gs$bodyid,meta1$bodyid)
  ids.missing = ids.missing[!is.na(ids.missing)]
  ids.missing = id_okay(ids.missing)
  if(length(ids.missing)){
    meta2 = hemibrain_get_meta(unique(ids.missing), ...)
    meta = rbind(meta1,meta2)
  }else{
    meta = meta1
  }
  meta$cell_type = meta$type
  chosen.cols = c("bodyid","ito_lee_hemilineage","hartenstein_hemilineage","cellBodyFiber","cell_type","layer","ct.layer","flywire_xyz","root_id","top.nblast")
  meta = meta[,intersect(colnames(meta),chosen.cols)]
  meta$fafb_match = gs$fafb_match[match(meta$bodyid,gs$bodyid)]
  meta$fafb_match_quality = gs$fafb_match_quality[match(meta$bodyid,gs$bodyid)]
  meta$lm_match = gs$lm_match[match(meta$bodyid,gs$bodyid)]
  meta$lm_match_quality = gs$lm_match_quality[match(meta$bodyid,gs$bodyid)]
  meta$user = gs$user[match(meta$bodyid,gs$bodyid)]
  meta = meta[order(meta$bodyid),]
  meta = meta[order(meta$cell_type),]
  meta = meta[order(meta$ito_lee_hemilineage),]
  if(top.nblast){
    nblast.hemibrain.catmaid = tryCatch(hemibrain_nblast('hemibrain-fafb14'), error = function(e) NULL)
    if(!is.null(nblast.hemibrain.catmaid)){
      nblast.top = nblast.hemibrain.catmaid[,match(meta$bodyid,colnames(nblast.hemibrain.catmaid))]
      tops = apply(nblast.top,1,function(r) which.max(r))
      top = rownames(nblast.hemibrain.catmaid)[unlist(tops)]
      top[!meta$bodyid%in%colnames(nblast.hemibrain.catmaid)] = NA
      meta$catmaid.nblast.top = top
    }
    nblast.hemibrain.flywire = tryCatch(hemibrain_nblast('hemibrain-flywire'), error = function(e) NULL)
    if(!is.null(nblast.hemibrain.flywire) & !is.null(meta)){
      nblast.top =nblast.hemibrain.flywire[,match(meta$bodyid,colnames(nblast.hemibrain.flywire))]
      tops = apply(nblast.top,1,function(r) which.max(r))
      top = rownames(nblast.hemibrain.flywire)[unlist(tops)]
      top[!meta$bodyid%in%colnames(nblast.hemibrain.flywire)] = NA
      meta$flywire.nblast.top = meta[match(unlist(top),meta$root_id),"flywire_xyz"]
    }
  }
  if(!identical(gs,meta)){
    write.cols = unique(colnames(meta),chosen.cols)
    for(wc in write.cols){
      gs[[wc]] = meta[match(gs$bodyid,meta$bodyid),wc]
    }
    gsheet_update_cols(
      write.cols = write.cols,
      gs=gs,
      selected_sheet = selected_file,
      sheet = "hemibrain")
  }
  if(nrow(meta)>nrow(gs)){
    missed.cols = setdiff(colnames(gs),colnames(meta))
    if(length(missed.cols)){
      for(mc in missed.cols){
        meta[[mc]] = NA
      }
    }
    meta = meta[,colnames(gs)]
    newrows = subset(meta, !meta$bodyid%in%gs$bodyid)
    gsheet_manipulation(FUN = googlesheets4::sheet_append,
                                     data = newrows,
                                     ss = selected_file,
                                     sheet = "hemibrain")
  }
  # batches = split(1:nrow(meta), ceiling(seq_along(1:nrow(meta))/500))
  # gsheet_manipulation(FUN = googlesheets4::write_sheet,
  #                     data = meta[0,],
  #                     ss = selected_file,
  #                     sheet = "hemibrain")
  # for(i in batches){
  #   gsheet_manipulation(FUN = googlesheets4::sheet_append,
  #                       data = meta[min(i):max(i),],
  #                       ss = selected_file,
  #                       sheet = "hemibrain")
  #}
}

#' @rdname hemibrain_add_made_matches
#' @export
flycircuit_matching_rewrite <- function(flycircuit.ids = names(flycircuit_neurons()),
                                        meta = flycircuit_neurons[,],
                                        selected_file  = options()$hemibrainr_matching_gsheet){

  # Get the LM matching Google sheet
  gs = hemibrain_match_sheet(sheet = "lm", selected_file = selected_file)
  ids = unique(gs$id)

  # Get meta
  shared.cols = intersect(colnames(gs),colnames(fc.meta))
  fc.meta = meta[,shared.cols]
  fc.meta = subset(fc.meta, fc.meta$id %in% flycircuit.ids)

  # gs merge
  if(nrow(gs)){
    fc.meta = subset(fc.meta, !fc.meta$id %in% gs$id)
    gs = plyr::rbind.fill(gs, fc.meta)
  }else{
    gs = fc.meta
  }

  # Update
  rownames(gs) = NULL
  googlesheets4::write_sheet(gs[0,],
                             ss = selected_file,
                             sheet = "lm")
  batches = split(1:nrow(gs), ceiling(seq_along(1:nrow(gs))/500))
  for(i in batches){
    gsheet_manipulation(FUN = googlesheets4::sheet_append,
                        data = gs[min(i):max(i),],
                        ss = selected_file,
                        sheet = "lm")
  }
}

# hidden
id_selector <- function(gs,
                        ids = NULL,
                        id = c("bodyid","root_id","skid","id","flywire_xyz","fafb_match","hemibrain_match","lm_match"),
                        overwrite = c("FALSE","bad","TRUE","review"),
                        quality.field,
                        match.field,
                        initials = NULL,
                        column = NULL,
                        entry = NULL,
                        superuser = FALSE,
                        flywire.good = FALSE){
  id = match.arg(id)
  ids = id_okay(ids)
  overwrite = match.arg(overwrite)
  if(is.null(column)+is.null(entry)==1){
    stop("column and entry must both be NULL, or both be given")
  }
  if("duplicated" %in% colnames(gs)){
    gs$duplicated = as.logical(gs$duplicated)
    gs = subset(gs, !gs$duplicated)
  }
  # choose possible ids
  id.len = ifelse(length(ids),length(ids),"all")
  if((overwrite %in% c("TRUE") & is.null(ids))){
    doit = gs
  }else if(overwrite=="bad"){
    doit = subset(gs, ( (is.na(gs[[match.field]]) | gs[[match.field]] == 0)
                  & is.na(gs[[quality.field]]))
                  | gs[[quality.field]] %in% c("none","n","tract","t",""," ","NA","0")
                  | gs[[match.field]] %in% c("none","n","tract","t",""," ","NA","0"))
    if(!nrow(doit)){
      stop("No valid 'bad' neurons for matching could be found. Try a different set of arguments")
    }
  }else if(overwrite=="FALSE"){
    doit = subset(gs, ( (is.na(gs[[match.field]]) | gs[[match.field]] == 0) & is.na(gs[[quality.field]])) )
  }else{ # review
    doit = subset(gs, ( (!is.na(gs[[match.field]]) | gs[[match.field]] != 0) | !is.na(gs[[quality.field]])) )
  }
  if(!superuser){
    doit = subset(doit, doit$user==initials)
    users = paste0("for matches flagged for user ", initials)
  }else{
    users = "for any match"
  }
  able = if(overwrite%in%c("bad","TRUE")){
    "enabled"
  }else{
    "disabled"
  }
  if(flywire.good & id%in%c("root_id","flywire_xyz")){
    message("Only taking flywire neurons of completion status adequate/complete")
    fw.meta = flywire_meta(sql=FALSE)
    good.ids = subset(fw.meta, grepl("adequate|complete",fw.meta$status))
    doit = subset(doit, doit[[id]]%in%good.ids$root_id)
    if(!nrow(doit)){
      stop("No valid neurons for matching could be found adequate/complete in our flywire meta data (flywire_meta()). Try a different set of arguments, e.g. flywire.good = FALSE")
    }
  }
  if(!is.null(column)){
    if(!is.null(entry)){
      doit = subset(doit, doit[[column]]%in%entry)
      if(!nrow(doit)){
        stop("No valid neurons for matching could be found with column: ", column, " with entry: ", entry," . Try a different set of arguments.")
      }
    }
  }
  message(sprintf("Reviewing made matches for %s %ss for user %s, with overwrite %s %s",
                  id.len,
                  id,
                  initials,
                  able,
                  users
                  ))
  if(is.null(ids)||!length(ids)){
    ids = doit[[id]]
  }else{
    ids = intersect(ids,doit[[id]])
    doit = subset(doit, doit[[id]]%in%ids)
  }
  ids = ids[!grepl("missing",ids)]
  ids = ids[!is.na(ids)]
  if(!nrow(doit)){
    stop(sprintf("No matches to be made with the given IDs and parameters -
id = '%s'; overwrite = '%s'; quality.field = '%s'; match.field = '%s'; initials = '%s'; column = '%s'; entry = '%s'",
                 id,
                 overwrite,
                 quality.field,
                 match.field,
                 initials,
                 nullToNA(column),
                 nullToNA(entry)))
  }else{
    message(sprintf("Out of %s possible %ss there are %s valid ones considering arguments user + superuser + overwrite + column + entry",id.len,id,length(ids)))
  }
  # order doit
  user.order = unique(c(initials,sort(doit$user)))
  if(length(user.order)>1){
    message(sprintf("We will go through matches in this user order: ", initials, paste(user.order,collapse=", ")))
    print(table(doit$user)[user.order])
  }
  doit = doit[order(match(doit$user,user.order)),]
  # return
  doit
}

# hidden
choose_user <- function(gs, user = NULL){
  users = sort(unique(gs$user))
  message("users: ", paste(users,collapse = " "))
  if(is.null(user)||!user%in%users){
    initials = must_be("Enter for new user. Or else, choose a user : ", answers = c(users,""))
  }else{
    initials = user
  }
  if(initials%in%c(""," ")){
    initials = readline(prompt="Please give your initials : ")
  }
  say_hello(initials)
  initials
}
