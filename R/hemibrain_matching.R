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
#' If \code{overwrite} is set to \code{FALSE}, you will not overwrite any matches that have already been made from among the selected neurons.
#' If \code{'mine'} you will re-examine and overwrite any matches you have made out of the selected neurons.
#' With \code{'mine_empty'} the same happens, but you will also retain neurons that have no match, and have been assigned to any other user.
#' If \code{FALSE} (be careful!) then you overwrite made matches among the selected neurons.
#' @param column defaults to \code{NULL}, no further subsetting. Else, you can select a column from the Google sheet.
#' Only  neurons with a certain value (\code{entry}) in that column will be chosen for matching.
#' @param entry defaults to \code{NULL}, no further subsetting. Else, it is a value in \code{column}.
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
                         threshold = 0,
                         selected_file = options()$hemibrainr_matching_gsheet,
                         batch_size = 50,
                         db=NULL, # brain="FAFB"
                         repository = c("flywire", "CATMAID", "lm"),
                         query = hemibrain_neurons(brain = "FAFB14"), # brain="FAFB"
                         overwrite = c("FALSE","mine","mine_empty","TRUE"),
                         column = NULL,
                         entry = NULL){
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
  # Read the Google Sheet
  gs = hemibrain_match_sheet(sheet = "hemibrain", selected_file = selected_file)
  gs$User[is.na(gs$User)] = ""
  # Get hemibrain neurons
  if(missing(query)) {
    query=tryCatch(force(query), error=function(e) {
        message("Unable to use `hemibrain_neurons(brain = 'FAFB14')`. ",
                "I will read neurons from neuPrint, but this will be slower!")
      })
  }else if(is.character(query)) {
    query=tryCatch(get(query), error=function(e) stop("Unable to find neuronlist: ", query))
  }
  if(repository=="flywire" & is.null(db)){
    db=tryCatch(flywire_neurons(), error=function(e) {
      message("Unable to use `flywire_neurons()`. ",
              "I will read neurons from FlyWire, but this will be slower!")
    })
  }else if(repository=="LM"){
    db=tryCatch(lm_lhns(brainspace = c("FAFB14")), error=function(e) {
      message("Unable to use `lm_lhns()`. ",
              "I will read neurons from package lhns, but this will be slower!")
    })
  }else{
    db = NULL
  }
  # fields
  if(repository=="CATMAID"){
    match.field = search.id = "FAFB.match"
    chosen.field = "skid"
    quality.field = "FAFB.match.quality"
    brain = elmr::FAFB14.surf
  }else if(repository == "flywire"){
    search.id = "flywire.id"
    match.field = chosen.field = "flywire.xyz"
    quality.field = "FAFB.match.quality"
    brain = elmr::FAFB14.surf
  }else{
    match.field = search.id = paste0(repository,".match")
    quality.field = paste0(repository,".match.quality")
    brain = hemibrainr::hemibrain_microns.surf
  }
  id = "bodyid"
  # How much is done?
  done = subset(gs, !is.na(gs[[match.field]]))
  message("Neuron matches: ", nrow(done))
  print(table(gs[[quality.field]]))
  # Choose user
  initials = choose_user(gs)
  rgl::bg3d("white")
  # choose ids
  selected = id_selector(gs=gs, ids=ids, id = id, overwrite = overwrite,
                         quality.field = quality.field, match.field = match.field,
                         initials = initials, column = column, entry = entry)
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
                                       soma.size = 400,
                                       show.columns = c("cell.type","ItoLee_Hemilineage","note"))
    selected = match_cycle[["selected"]]
    unsaved = match_cycle[["unsaved"]]
    if(length(unsaved)){
      plot_inspirobot()
      say_encouragement(initials)
      # Read!
      gs2 = hemibrain_match_sheet(selected_file = selected_file, sheet = "hemibrain")
      selected.unsaved = subset(selected, selected[[id]]%in%unsaved)
      gs2[match(selected.unsaved[[id]],gs2[[id]]),match.field]= selected.unsaved[[match.field]]
      gs2[match(selected.unsaved[[id]],gs2[[id]]),quality.field]= selected.unsaved[[quality.field]]
      gs2[match(selected.unsaved[[id]],gs2[[id]]),"flywire.id"]= selected.unsaved[["flywire.id"]]
      gs2[match(selected.unsaved[[id]],gs2[[id]]),"note"]= selected.unsaved[["note"]]
      gs2[match(selected.unsaved[[id]],gs2[[id]]),"User"]= initials
      # Write!
      write_matches(gs=gs2,
                    ids = unsaved,
                    id.field = id,
                    column = match.field,
                    selected_file = selected_file,
                    ws = "hemibrain")
      write_matches(gs=gs2,
                    ids = unsaved,
                    id.field = id,
                    column = quality.field,
                    selected_file = selected_file,
                    ws = "hemibrain")
      write_matches(gs=gs2,
                    ids = unsaved,
                    id.field = id,
                    column = "note",
                    selected_file = selected_file,
                    ws = "hemibrain")
      write_matches(gs=gs2,
                    ids = unsaved,
                    id.field = id,
                    column = "User",
                    selected_file = selected_file,
                    ws = "hemibrain")
      write_matches(gs=gs2,
                    ids = unsaved,
                    id.field = id,
                    column = "flywire.id",
                    selected_file = selected_file,
                    ws = "hemibrain")
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
                        overwrite = c("FALSE","mine","mine_empty","TRUE"),
                        column = NULL,
                        entry = NULL){
  # Motivate!
  nat::nopen3d()
  plot_inspirobot()
  message("
          #######################Colours##########################
          blue = LM neuron,
          red = potential hemibrain matches based on NBLAST score,
          green = a chosen hemibrain neuron during scanning,
          dark blue = your selected hemibrain match.
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
  gs$User[is.na(gs$User)] = ""
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
  quality.field = paste0("hemibrain",".match.quality")
  done = subset(gs, !is.na(gs[[match.field]]))
  message("Neuron matches: ", nrow(done), "/", nrow(gs))
  print(table(gs[[quality.field]]))
  # Choose user
  initials = choose_user(gs)
  rgl::bg3d("white")
  # choose ids
  selected = id_selector(gs=gs, ids=ids, id=id, overwrite = overwrite,
                         quality.field = quality.field, match.field = match.field,
                         initials = initials, column = column, entry = entry)
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
                        overwrite = c("FALSE","mine","mine_empty","TRUE"),
                        column = NULL,
                        entry = NULL){
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
          red = potential hemibrain matches based on NBLAST score,
          green = a chosen hemibrain neuron during scanning,
          dark blue = your selected hemibrain match.
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
    extra.neurons = tryCatch(flywire_neurons(), error = function(e) NULL)
  }else if(repository == "flywire"){
    id = "flywire.id"
    brain = elmr::FAFB14.surf
    if(is.null(ids)){
      ids = tryCatch(names(flywire_neurons()), error = function(e) NULL)
    }
    if(is.null(query)){
      query = tryCatch(flywire_neurons(), error = function(e) NULL)
    }
    extra.neurons = tryCatch(flywire_neurons(mirror=TRUE), error = function(e) NULL)
  }
  chosen.field = "bodyid"
  match.field = search.id = paste0("hemibrain",".match")
  quality.field = paste0("hemibrain",".match.quality")
  # Read the Google Sheet
  gs = hemibrain_match_sheet(selected_file = selected_file, sheet = repository)
  gs$User[is.na(gs$User)] = ""
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
  initials = choose_user(gs)
  # choose ids
  selected = id_selector(gs=gs, ids=ids, id=id, overwrite = overwrite,
                         quality.field = quality.field, match.field = match.field,
                         initials = initials, column = column, entry = entry)
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
                                       soma.size = 400,
                                       show.columns = c("cell.type","ItoLee_Hemilineage","note"))
    selected = match_cycle[["selected"]]
    unsaved = match_cycle[["unsaved"]]
    if(length(unsaved)){
      plot_inspirobot()
      say_encouragement(initials)
      # Read!
      gs2 = hemibrain_match_sheet(selected_file = selected_file, sheet = "hemibrain")
      selected.unsaved = subset(selected, selected[[id]]%in%unsaved)
      gs2[match(selected.unsaved[[id]],gs2[[id]]),match.field]= selected.unsaved[[match.field]]
      gs2[match(selected.unsaved[[id]],gs2[[id]]),quality.field]= selected.unsaved[[quality.field]]
      gs2[match(selected.unsaved[[id]],gs2[[id]]),"note"]= selected.unsaved[["note"]]
      gs2[match(selected.unsaved[[id]],gs2[[id]]),"User"]= initials
      # Write!
      write_matches(gs=gs2,
                    ids = unsaved,
                    id.field = id,
                    column = match.field,
                    selected_file = selected_file,
                    ws = "FAFB")
      write_matches(gs=gs2,
                    ids = unsaved,
                    id.field = id,
                    column = quality.field,
                    selected_file = selected_file,
                    ws = "FAFB")
      write_matches(gs=gs2,
                    ids = unsaved,
                    id.field = id,
                    column = "note",
                    selected_file = selected_file,
                    ws = "FAFB")
      write_matches(gs=gs2,
                    ids = unsaved,
                    id.field = id,
                    column = "User",
                    selected_file = selected_file,
                    ws = "FAFB")
      saved = c(unsaved, saved)
      unsaved = c()
    }
    match.more = hemibrain_choice("Match more neurons? y/n ")
  }
  say_encouragement(initials)
}

# neuron ID name
get_idfield <- function(repository = c("hemibrain","CATMAID","flywire","LM","lm","FAFB.hemisphere","FAFB"),
                        return = c("id.field","sheet","match.field","match.quality")){
  return = match.arg(return)
  repository = match.arg(repository)
  sheet = repository
  if(repository=="hemibrain"){
    id.field = "bodyid"
  }else if (repository == c("flywire")){
    id.field = "flywire.id"
    sheet = "FAFB"
  }else if (repository%in%c("CATMAID","FAFB")){
    id.field = "skid"
    sheet = "FAFB"
  }else if (repository=="FAFB"){
    id.field = "skid"
  }else if(repository=="FAFB.hemisphere"){
    id.field = "FAFB.hemisphere.match"
  }else{
    id.field = "id"
  }
  if(return=="sheet"){
    sheet
  }else if(return%in%c("match.field","match.quality")){
    if(repository%in%c("CATMAID","FAFB")){
      match.field  = "FAFB.match"
      match.quality = "FAFB.match.quality"
    }else if(repository%in%c("LM","lm")){
      match.field = "LM.match"
      match.quality = "LM.match.quality"
    }else if(repository=="hemibrain"){
      match.field = "hemibrain.match"
      match.quality = "hemibrain.match.quality"
    }else if(repository=="flywire"){
      match.field = "flywire.xyz"
      match.quality = "FAFB.match.quality"
    }else if(repository%in%c("FAFB.hemisphere","hemisphere")){
      match.field = "FAFB.hemisphere.match"
      match.quality = "FAFB.hemisphere.match.quality"
    }else{
      match.field = "match"
      match.quality = "quality"
    }
    if(return=="match.field"){
      match.field
    }else{
      match.quality
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
                          User){
  if(length(User)>1){
    stop("Specify just one user")
  }
  id.field = get_idfield(repository, return = "id.field")
  ws = get_idfield(repository, return = "sheet")
  gs = hemibrain_match_sheet(selected_file = selected_file, sheet = repository)
  gs[gs[[id.field]]%in%ids,"User"] = User
  write_matches(gs=gs,
                ids = ids,
                id.field = id.field,
                column = "User",
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
#' dependent on the degree of manual tracing for FAFB neurons - good (could be the same cell), medium (same cell type) and poor (could be the same or similar cell type).
#'
#' @param priority whether to use FAFB->hemibrain matches (FAFB) or hemibrain->FAFB matches (hemibrain) in order to ascribe
#' cell type names to FAFB neurons. In both cases, cell type names are attached to hemibrain bodyids, and propagated to their FAFB matches.
#' @param selected_file the Google Sheet database to read and write from. For
#'   now, defaults to this
#'   \href{https://docs.google.com/spreadsheets/d/1OSlDtnR3B1LiB5cwI5x5Ql6LkZd8JOS5bBr-HTi0pOw/edit#gid=0}{Google
#'    Sheet}.
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
#' @return a \code{data.frame} where each row is a neuron, either from the hemibrain or FAFB data sets. Each row gives you its matching neuron in the other data set. These matches have been
#' manually assigned using \code{\link{fafb_matching}}, \code{\link{hemibrain_matching}} and \code{\link{LR_matching}}. If you use this information make sure you credit it appropriately.
#' Contact us if unsure:
#' \itemize{
#'
#'   \item{"cell.type"}{ - the neuPrint designated 'type' for the neuron. If \code{dataset} is not \code{"hemibrain"}, then this is based on the hemibrain \code{match}.}
#'
#'   \item{"cell"}{ - the unique cell, which is just \code{cell.type#number}.}
#'
#'   \item{"cellBodyFiber"}{ - the cell body fiber to which this neuron belongs}
#'
#'   \item{"ItoLee_Hemilineage"}{ - the hemilineage to which this neuron belongs. Seer \code{\link{hemibrain_hemilineages}}.}
#'
#'   \item{"match"}{ - the ID of the manual match from the other data set. If \code{dataset=="hemibrain"} then this is a \code{flywire.id} that can be found in \code{flywire_neurons}.If \code{"CATMAID"} or \code{"flywire"} then it is a hemibrain body ID.}
#'
#'   \item{"quality"}{ - the matcher makers qualitative assessment of how good this match is.}
#'
#'   \item{"FAFB.hemisphere.match"}{ - the flywire coordinates of a neuron on the opposite hemisphere, which a match maker has designated as this \code{id}'s cognate.}
#'
#'   \item{"FAFB.hemisphere.match.quality"}{ - the quality of this match.}
#'
#'   \item{"LM.match"}{ - indicates a light level neuron that is a match for \code{id}. This neuron will be in \code{flycircuit_neurons()} or other light level data.}
#'
#'   \item{"LM.match.quality"}{ - the quality of this match.}
#'
#'   \item{"dataset"}{ - the data set to which \code{id} belongs. Either \code{"hemibrain"}, or one of the two FAFB repositories, \code{"CATMAID"} or \code{"flywire"}.}
#'
#'   \item{"priority"}{ - whether FAFB->hemibrain matches (\code{"FAFB"}) or hemibrain->FAFB matches (\code{"hemibrain"}) were used in order to ascribe cell type names to FAFB neurons.
#'   In both cases, cell type names are attached to hemibrain bodyids, and propagated to their FAFB matches.}
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
hemibrain_matches <- function(priority = c("FAFB","hemibrain"),
                              selected_file = options()$hemibrainr_matching_gsheet){
  priority = match.arg(priority)

  # Get matches
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

  # Get FAFB matches
  fafb.matches = gsheet_manipulation(FUN = googlesheets4::read_sheet,
                                                  ss = selected_file,
                                                  sheet = "FAFB",
                                                  return = TRUE)
  fafb.matches = subset(fafb.matches, !is.na(fafb.matches$skid)|!is.na(fafb.matches$flywire.id))
  fafb.matches = fafb.matches[!(duplicated(fafb.matches$skid)&duplicated(fafb.matches$flywire.id)),]
  fafb.matches$skid = correct_id(fafb.matches$skid)
  fafb.matches$flywire.id = correct_id(fafb.matches$flywire.id)
  fafb.matches$dataset = "FAFB"
  fafb.matches$skid[is.na(fafb.matches$skid)] = "missing"
  rownames(fafb.matches) = paste0(fafb.matches$skid,"#",ave(fafb.matches$skid,fafb.matches$skid,FUN= seq.int))

  # Set unassigned matches to top.nblast
  nblast.tops = c("CATMAID.nblast.top","flywire.nblast.top")
  for(nt in nblast.tops){
    if(!is.null(hemibrain.matches[[nt]])){
      present = !is.na(hemibrain.matches[[nt]])
      missing = is.na(hemibrain.matches$flywire.xyz)
      replace = (present+missing)>1
      if(sum(replace)){
        if(nt=="CATMAID.nblast.top"){
          hemibrain.matches$FAFB.match[replace] = hemibrain.matches[[nt]][replace]
          hemibrain.matches$flywire.xyz[replace] = fafb.matches$flywire.xyz[match(hemibrain.matches$FAFB.match[replace],fafb.matches$skid)]
          hemibrain.matches$flywire.id[replace] = fafb.matches$flywire.id[match(hemibrain.matches$FAFB.match[replace],fafb.matches$skid)]
        }else{
          hemibrain.matches$flywire.id[replace] = hemibrain.matches[[nt]][replace]
          hemibrain.matches$flywire.xyz[replace] = fafb.matches$flywire.xyz[match(fafb.matches$flywire.id[replace],fafb.matches$flywire.id)]
          hemibrain.matches$FAFB.match[replace] = fafb.matches$skid[match(hemibrain.matches$flywire.id[replace],fafb.matches$flywire.id)]
        }
        hemibrain.matches$FAFB.match.quality[replace] = "NBLAST"
      }
      }
  }

  # Set unassigned matches to top.nblast
  nblast.tops = c("hemibrain.nblast.top")
  for(nt in nblast.tops){
    if(!is.null(fafb.matches[[nt]])){
      present = !is.na(fafb.matches[[nt]])
      missing = is.na(fafb.matches$hemibrain.match)
      replace = (present+missing)>1
      if(sum(replace)){
        fafb.matches$hemibrain.match[replace] = fafb.matches[[nt]][replace]
        fafb.matches$hemibrain.match.quality[replace] = "NBLAST"
      }
    }
  }

  # Match to other side
  side.match = fafb.matches$FAFB.hemisphere.match[match(hemibrain.matches$flywire.xyz,fafb.matches$flywire.xyz)]
  side.quality = fafb.matches$FAFB.hemisphere.match.quality[match(hemibrain.matches$flywire.xyz,fafb.matches$flywire.xyz)]
  hemibrain.matches$FAFB.hemisphere.match = ifelse(length(side.match),side.match,NA)
  hemibrain.matches$FAFB.hemisphere.match.quality = ifelse(length(side.quality),side.quality,NA)

  # Unmatched
  fafb.unmatched = is.na(fafb.matches$hemibrain.match)|fafb.matches$hemibrain.match==""
  fafb.matches$hemibrain.match.quality[fafb.unmatched] = "none"
  fafb.matches$hemibrain.match[fafb.unmatched] = "none"
  hemibrain.unmatched = is.na(hemibrain.matches$FAFB.match)|hemibrain.matches$FAFB.match==""
  hemibrain.matches$FAFB.match.quality[hemibrain.unmatched] = "none"
  hemibrain.matches$FAFB.match[hemibrain.unmatched] = "none"
  side.unmatched = is.na(hemibrain.matches$FAFB.hemisphere.match)|hemibrain.matches$FAFB.hemisphere.match==""
  hemibrain.matches$FAFB.hemisphere.match.quality[side.unmatched] = "none"
  hemibrain.matches$FAFB.hemisphere.match[side.unmatched] = "none"
  side.unmatched = is.na(fafb.matches$FAFB.hemisphere.match)|fafb.matches$FAFB.hemisphere.match==""
  fafb.matches$FAFB.hemisphere.match.quality[side.unmatched] = "none"
  fafb.matches$FAFB.hemisphere.match[side.unmatched] = "none"

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
      inq = as.character(subset(fafb.matches, fafb.matches$hemibrain.match.quality==q)$flywire.xyz)
      inq = inq[!is.na(inq)]
      isna = is.na(fafb.matches[inq,]$cell.type)
      cts = hemibrain.matches$cell.type[match(fafb.matches[inq,]$hemibrain.match,hemibrain.matches$bodyid)]
      fafb.matches[inq,]$cell.type[isna] = cts[isna]
      # Hemibrain -> FAFB
      inq = as.character(subset(hemibrain.matches, hemibrain.matches$FAFB.match.quality==q)$flywire.xyz)
      inq = inq[!is.na(inq)]
      isna = is.na(fafb.matches[inq,]$cell.type)
      fafb.matches[inq,]$cell.type[isna] = hemibrain.matches$cell.type[match(fafb.matches[inq,]$flywire.xyz[isna],hemibrain.matches$flywire.xyz)]
    }else{
      # Hemibrain -> FAFB
      inq = as.character(subset(hemibrain.matches, hemibrain.matches$FAFB.match.quality==q)$flywire.xyz)
      inq = inq[!is.na(inq)]
      isna = is.na(fafb.matches[inq,]$cell.type)
      fafb.matches[inq,]$cell.type[isna] = hemibrain.matches$cell.type[match(fafb.matches[inq,]$flywire.xyz[isna],hemibrain.matches$flywire.xyz)]
      # FAFB -> Hemibrain
      inq = as.character(subset(fafb.matches, fafb.matches$hemibrain.match.quality==q)$flywire.xyz)
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
  for(repo in c("CATMAID","flywire")){
    id.field = get_idfield(repository = repo, return = "match.field")
    match.field = get_idfield(repository = repo, return = "match.field")
    ids = id_okay(fafb.matches[[id.field]])
    for(id in ids){
      idx = match(id,fafb.matches[[id.field]])
      idx2 = match(id,hemibrain.matches[[match.field]])
      if(!is.na(idx2)){
        ct = fafb.matches[idx,"cell.type"]
        if(is.na(ct)){
          fafb.matches[idx,"cell.type"] = hemibrain.matches$cell.type[idx2]
          if(is.na(fafb.matches[idx,"cell.type"])){
            fafb.matches[idx,"cell.type"] = "uncertain"
          }else{
            fafb.matches[idx,"ItoLee_Hemilineage"] = hemibrain.matches$ItoLee_Hemilineage[idx2]
          }
        }else{
          fafb.matches[idx,"ItoLee_Hemilineage"] = hemibrain.matches$ItoLee_Hemilineage[match(ct,hemibrain.matches$cell.type)]
        }
      }
    }
  }

  # Add cell body fiber infor for FAFB cells
  cbf.match = hemibrain.matches$cellBodyFiber[match(fafb.matches$hemibrain.match,hemibrain.matches$bodyid)]
  fafb.matches$cellBodyFiber = ifelse(length(cbf.match),cbf.match,NA)

  # Rename cells
  fafb.matches = fafb.matches[order(fafb.matches$hemibrain.match,decreasing = TRUE),]
  hemibrain.matches = hemibrain.matches[order(hemibrain.matches$bodyid,decreasing = TRUE),]
  fafb.matches$cell = paste0(fafb.matches$cell.type,"#",ave(fafb.matches$cell.type,fafb.matches$cell.type,FUN= seq.int))
  hemibrain.matches$cell = paste0(hemibrain.matches$cell.type,"#",ave(hemibrain.matches$cell.type,hemibrain.matches$cell.type,FUN= seq.int))

  # Separate flywire and CATMAID matches
  catmaid.matches = subset(fafb.matches, !is.na(fafb.matches$skid)& !fafb.matches$skid%in%c(""," ","NA","0","error","none"))
  flywire.matches = subset(fafb.matches, !is.na(fafb.matches$flywire.id)& !fafb.matches$flywire.id%in%c(""," ","NA","0","error","none"))
  hemibrain.matches$dataset = "hemibrain"
  flywire.matches$dataset = "flywire"
  catmaid.matches$dataset = "CATMAID"
  hemibrain.matches$match.dataset = "flywire"
  flywire.matches$match.dataset = "hemibrain"
  catmaid.matches$match.dataset = "hemibrain"
  hemibrain.matches.catmaid = hemibrain.matches
  hemibrain.matches.catmaid$match.dataset = "CATMAID"

  # Make matching data frame
  matched.h = hemibrain.matches[,c("bodyid", "cell.type", "cell", "cellBodyFiber", "ItoLee_Hemilineage",
                                   "flywire.id", "FAFB.match.quality",
                                   "FAFB.hemisphere.match", "FAFB.hemisphere.match.quality",
                                   "LM.match", "LM.match.quality", "dataset", "match.dataset")]
  matched.hc = hemibrain.matches.catmaid[,c("bodyid", "cell.type", "cell", "cellBodyFiber", "ItoLee_Hemilineage",
                                   "FAFB.match", "FAFB.match.quality",
                                   "FAFB.hemisphere.match", "FAFB.hemisphere.match.quality",
                                   "LM.match", "LM.match.quality", "dataset", "match.dataset")]
  matched.c = catmaid.matches[,c("skid",  "cell.type",  "cell", "cellBodyFiber", "ItoLee_Hemilineage",
                              "hemibrain.match", "hemibrain.match.quality",
                              "FAFB.hemisphere.match", "FAFB.hemisphere.match.quality",
                              "LM.match", "LM.match.quality","dataset", "match.dataset")]
  matched.f = flywire.matches[,c("flywire.id",  "cell.type",  "cell", "cellBodyFiber", "ItoLee_Hemilineage",
                              "hemibrain.match", "hemibrain.match.quality",
                              "FAFB.hemisphere.match", "FAFB.hemisphere.match.quality",
                              "LM.match", "LM.match.quality","dataset", "match.dataset")]
  colnames(matched.h) = colnames(matched.hc) = colnames(matched.c) = colnames(matched.f) = c("id","cell.type", "cell","cellBodyFiber","ItoLee_Hemilineage",
                                                "match","quality",
                                                "FAFB.hemisphere.match", "FAFB.hemisphere.match.quality",
                                                "LM.match", "LM.match.quality","dataset", "match.dataset")
  matched = rbind(matched.h,matched.hc,matched.f,matched.c)
  matched$quality[is.na(matched$match)] = "none"
  matched$match[is.na(matched$match)] = "none"
  matched$ItoLee_Lineage = gsub("_.*","",matched$ItoLee_Hemilineage)

  # Sort out types
  matched$connectivity.type = matched$cell.type
  matched$cell.type = gsub("_[a-z]{1}$","",matched$cell.type)
  matched$cell.type[is.na(matched$cell.type)] = "uncertain"
  matched$connectivity.type[is.na(matched$connectivity.type)] = "uncertain"
  matched$priority = priority
  rownames(matched) = 1:nrow(matched)

  # Return
  matched
}

#' @rdname hemibrain_matches
#' @export
lm_matches <- function(priority = c("hemibrain","lm"), selected_file = options()$hemibrainr_matching_gsheet){
  priority = match.arg(priority)

  # Get matches
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
#' @param flycircuit.ids flycircuit IDs to add to Google sheet if not already present.
#' @param meta meta data for the given flycircuit IDs.
#' @param top.nblast logical. Whether or not to also give the top NBLAST match for each entry.
#' @param catmaid.update logical. Whether or not to update \code{flywire.xyz} and \code{flywire.id} columns, based on e CATMAID neuron specified by a \code{skid} column.
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
#' hemibrain_matching_add(ids = "16", dataset = "FAFB", User = "ASB")
#'
#'
#' }}
#' @seealso \code{\link{hemibrain_matching}}, \code{\link{fafb_matching}}
#' @rdname hemibrain_add_made_matches
hemibrain_add_made_matches <- function(df,
                                  direction = c("both","hemibrain-FAFB","FAFB-hemibrain"),
                                  User = "flyconnectome",
                                  selected_file = options()$hemibrainr_matching_gsheet,
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
      hemibrain_matching_add(ids = missing, dataset = "hemibrain", User = User, selected_file = selected_file)
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
                  selected_file = selected_file,
                  column = "FAFB.match")
    write_matches(gs=gs,
                  ids = as.character(hdf$bodyid),
                  ws="hemibrain",
                  id.field ="bodyid",
                  selected_file = selected_file,
                  column = "FAFB.match.quality")
  }
  if(direction%in%c("both","FAFB-hemibrain")){
    message("Reading FAFB sheet")
    gs = hemibrain_match_sheet(sheet = "FAFB")
    missing = setdiff(df$skid,gs$skid)
    if(length(missing)){
      message("Adding missing FAFB bodyids")
      hemibrain_matching_add(ids = missing, dataset = "FAFB", User = User, selected_file = selected_file)
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
                  selected_file = selected_file,
                  column = "hemibrain.match")
    write_matches(gs=gs,
                  ids = as.character(hdf$skid),
                  ws="FAFB",
                  id.field ="skid",
                  selected_file = selected_file,
                  column = "hemibrain.match.quality")
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
                                        return = TRUE)
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
  if(length(add)){
    # Meta information
    if(is.null(meta)){
      if(dataset=="hemibrain"){
        meta = hemibrain_get_meta(add, ...)
        meta$cell.type = meta$type
      }else if (dataset == "FAFB"){
        meta = elmr::fafb_get_meta(add, ...)
      } else if (dataset == "flywire"){
        add = setdiff(ids, gs$flywire.id)
        if(is.null(meta)){
          meta = flywire_neurons()[as.character(add),]
        }else{
          meta = meta[meta$flywire.id%in%add,]
        }
        meta = meta[!meta$flywire.xyz%in%gs$flywire.xyz,]
        if(!nrow(meta)){
          stop("Selected IDs could not be added. They must be among the neurons
         saved on Google drive, see flywire_neurons()")
        }
      }
    }else{
      id = if(dataset%in%c("FAFB","CATMAID")){
        "skid"
      }else if(dataset=="flywire"){
        "flywire.id"
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
      meta$User = User

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
  for(repo1 in c("CATMAID","flywire","hemibrain","lm")){
    for(repo2 in c("CATMAID","flywire","hemibrain","lm")){
      ws.1 = get_idfield(repo1, return = "sheet")
      ws.2 = get_idfield(repo2, return = "sheet")
      if(repo1==repo2|ws.1==ws.2){
        next
      }else{
        gs.1 = hemibrain_match_sheet(selected_file = selected_file, sheet = ws.1)
        gs.2 = hemibrain_match_sheet(selected_file = selected_file, sheet = ws.2)
        used.cols = colnames(gs.1)
        match.field.2 = get_idfield(repo2, return = "match.field")
        match.quality.2 = get_idfield(repo2, return = "match.quality")
        id.field.1 = get_idfield(repo1, return = "id.field")
        match.field.1 = get_idfield(repo1, return = "match.field")
        match.quality.1 = get_idfield(repo1, return = "match.quality")
        id.field.2 = get_idfield(repo2, return = "id.field")
        missing = is.na(gs.1[[match.field.2]])| (gs.1[[match.field.2]] %in% c("none"," ",""))
        matches = gs.2[[id.field.2]][match(gs.1[[id.field.1]][missing],gs.2[[match.field.1]])]
        if(length(matches)&sum(missing)){
          quality = gs.2[[match.quality.1]][match(gs.1[[id.field.1]][missing],gs.2[[match.field.1]])]
          matches[is.na(quality)] = NA
          quality[is.na(matches)] = NA
          gs.1[missing,match.field.2] = matches
          gs.1[missing,match.quality.2] = quality
          write.cols = intersect(c(match.field.2,match.quality.2),used.cols)
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
                          data = as.data.frame(update[row,column], stringsAsFactors = FALSE),
                          sheet = tab,
                          col_names = FALSE)
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
  n[match(n2$skid,n$skid),c("skid","ItoLee_Hemilineage", "Hartenstein_Hemilineage", "cell_body_fiber")] = n2[,c("skid","ItoLee_Hemilineage", "Hartenstein_Hemilineage", "cell_body_fiber")]
  ids.missing = as.character(setdiff(n1$skid,n$skid))
  ids.missing = id_okay(ids.missing)
  if(length(ids.missing)){
    n3 = elmr::fafb_get_meta(unique(ids.missing), batch = 10, ...)
    n = plyr::rbind.fill(n, n3[,c("skid","ItoLee_Hemilineage", "Hartenstein_Hemilineage", "cell_body_fiber")])
    reorder = TRUE
  }
  matches = tryCatch(hemibrain_matches(selected_file=selected_file), error = function(e) NULL)
  matches = subset(matches, ! matches$quality %in% "NBLAST")
  if(!is.null(matches)){
    n$cell.type = matches[match(n$skid, matches$id),"connectivity.type"]
  }else{
    n$cell.type = NULL
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
    n = n[order(n$cell.type),]
    n = n[order(n$ItoLee_Hemilineage),]
    n = n[!duplicated(n),]
    gsheet_manipulation(FUN = googlesheets4::write_sheet,
                        data = n[0,],
                        ss = selected_file,
                        sheet = "FAFB")
    batches = split(1:nrow(n), ceiling(seq_along(1:nrow(n))/500))
    for(i in batches){
      gsheet_manipulation(FUN = googlesheets4::sheet_append,
                          data = n[min(i):max(i),],
                          ss = selected_file,
                          sheet = "FAFB")
    }
  }else{
    # Update
    write.cols = intersect(c("connectivity.type","cell.type","nblast.top","side",
                             "ItoLee_Hemilineage","Hartenstein_Hemilineage","cell_body_fiber"),
                           colnames(n))
    gsheet_update_cols(
      write.cols = write.cols,
      gs=n,
      selected_sheet = selected_file,
      sheet = "FAFB")
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
  meta$cell.type = meta$type
  chosen.cols = c("bodyid","ItoLee_Hemilineage","Hartenstein_Hemilineage","cellBodyFiber","cell.type","layer","ct.layer","flywire.xyz","flywire.id","top.nblast")
  meta = meta[,intersect(colnames(meta),chosen.cols)]
  meta$FAFB.match = gs$FAFB.match[match(meta$bodyid,gs$bodyid)]
  meta$FAFB.match.quality = gs$FAFB.match.quality[match(meta$bodyid,gs$bodyid)]
  meta$LM.match = gs$LM.match[match(meta$bodyid,gs$bodyid)]
  meta$LM.match.quality = gs$LM.match.quality[match(meta$bodyid,gs$bodyid)]
  meta$User = gs$User[match(meta$bodyid,gs$bodyid)]
  meta = meta[order(meta$bodyid),]
  meta = meta[order(meta$cell.type),]
  meta = meta[order(meta$ItoLee_Hemilineage),]
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
      meta$flywire.nblast.top = meta[match(unlist(top),meta$flywire.id),"flywire.xyz"]
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
                        id = c("bodyid","flywire.id","skid","id","flywire.xyz","FAFB.match","hemibrain.match","LM.match"),
                        overwrite = c("FALSE","mine","mine_empty","TRUE"),
                        quality.field,
                        match.field,
                        initials = NULL,
                        column = NULL,
                        entry = NULL){
  id = match.arg(id)
  ids = id_okay(ids)
  overwrite = match.arg(overwrite)
  if(is.null(column)+is.null(entry)==1){
    stop("column and entry must both be NULL, or both be given")
  }
  # choose possible ids
  if(overwrite %in% c("TRUE") | is.null(initials)){
    doit = gs
  }else if(overwrite=="mine"){
    doit = subset(gs, gs$User == initials)
  } else if(overwrite=="mine_empty"){
    doit = subset(gs, gs$User == initials
                  | (is.na(gs[[quality.field]])
                  | is.na(gs[[quality.field]])
                  | gs[[quality.field]] %in% c("none","n","tract","t",""," ","NA")
                  | gs[[match.field]] %in% c("none","n","tract","t",""," ","NA")))
  }else if(overwrite=="FALSE"){
    doit = subset(gs, gs$User == initials & (is.na(gs[[match.field]]) | is.na(gs[[quality.field]])) )
  }
  if(is.null(ids)||!length(ids)){
    ids = doit[[id]]
  }else{
    ids = intersect(ids,doit[[id]])
  }
  ids = ids[!grepl("missing",ids)]
  ids = ids[!is.na(ids)]
  # further narrow
  selected = subset(doit, doit[[id]] %in% ids)
  if(!is.null(column)+!is.null(entry)==2){
    selected = subset(selected, selected[[column]]%in%entry)
  }
  if(!nrow(selected)){
    stop(sprintf("No matches to be made with the given IDs and parameters -
id = '%s'; overwrite = '%s'; quality.field = '%s'; match.field = '%s'; initials = '%s'; column = '%s'; entry = '%s'",
                 id,
                 overwrite,
                 quality.field,
                 match.field,
                 initials,
                 nullToNA(column),
                 nullToNA(entry)))
  }
  # return
  selected
}

# hidden
choose_user <- function(gs){
  message("Users: ", paste(sort(unique(gs$User)),collapse = " "))
  initials = must_be("Enter for new user. Or else, choose a user : ", answers = c(sort(unique(gs$User)),""))
  if(initials==""){
    initials = readline(prompt="Please give your initials : ")
  }
  say_hello(initials)
  initials
}



