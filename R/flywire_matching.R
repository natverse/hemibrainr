# FlyWire neuron matching

#' @rdname hemibrain_add_made_matches
#' @export
flywire_matching_rewrite <- function(flywire.ids = names(flywire_neurons()),
                                     selected_file  = options()$hemibrainr_matching_gsheet,
                                     ...){
  # Get the FAFB matching Google sheet
  gs = hemibrain_match_sheet(sheet = "FAFB", selected_file = selected_file)
  skids = unique(gs$skid)

  # Get FAFBv14 coordinates
  cats = nat::neuronlist()
  batches = split(1:length(skids), round(seq(from = 1, to = 100, length.out = length(skids))))
  all.ids = c()
  for(i in 1:10){
    # Read CATMAID neurons
    message("Batch:", i, "/100")
    cat = catmaid::read.neurons.catmaid(skids[batches[[i]]], OmitFailures = TRUE)
    cats = union(cats,cat)

    # Get xyz for primary branch points
    simp = nat::nlapply(cat,nat::simplify_neuron,n=1, .parallel = TRUE, OmitFailures = TRUE)
    branchpoints = sapply(simp, function(y) nat::xyzmatrix(y)[ifelse(length(nat::branchpoints(y)),nat::branchpoints(y),max(nat::endpoints(y))),])
    branchpoints = t(branchpoints)
    FAFB.xyz = apply(branchpoints, 1, paste, collapse = ";")

    # Get FlyWire voxel coordinates
    branchpoints.flywire = nat.templatebrains::xform_brain(branchpoints, reference = "FlyWire", sample = "FAFB14", .parallel = TRUE, verbose = TRUE)
    rownames(branchpoints.flywire) = rownames(branchpoints)
    branchpoints.flywire.raw = scale(branchpoints.flywire, scale = c(4, 4, 40), center = FALSE)
    fw.ids = fafbseg::flywire_xyz2id(branchpoints.flywire.raw, rawcoords = TRUE)
    fw.ids[fw.ids=="0"] = NA
    flywire.xyz = apply(branchpoints.flywire.raw, 1, paste, collapse = ";")

    # Add
    gs[rownames(branchpoints),]$FAFB.xyz = FAFB.xyz
    gs[rownames(branchpoints),]$flywire.xyz = flywire.xyz
    gs[rownames(branchpoints),]$flywire.id = fw.ids
    all.ids=unique(c(all.ids,fw.ids))
  }

  # Update
  rownames(gs) = NULL
  googlesheets4::write_sheet(gs[0,],
                             ss = selected_file,
                             sheet = "FAFB")
  batches = split(1:nrow(gs), ceiling(seq_along(1:nrow(gs))/500))
  for(i in batches){
    gsheet_manipulation(FUN = googlesheets4::sheet_append,
                        data = gs[min(i):max(i),],
                        ss = selected_file,
                        sheet = "FAFB")
  }

  # Add flywire match to the hemibrain and LM sheets
  ## Read the FAFB Google Sheet
  fg = hemibrain_match_sheet(sheet = "FAFB", selected_file = selected_file)

  ## Read the LM Google Sheet
  lmg = hemibrain_match_sheet(sheet = "lm", selected_file = selected_file)
  orig = lmg$flywire.xyz
  lmg$flywire.xyz = fg$flywire.xyz[match(lmg$bodyid,fg$hemibrain.match)]
  different = paste(orig)!=paste(lmg$flywire.xyz)
  lmg$flywire.xyz = fg$flywire.xyz[match(lmg$id,fg$LM.match)]
  update_gsheet(update = lmg[different,],
                  gs = lmg,
                  tab = "lm",
                  match = "flywire",
                  id = "id")

  ## Read the hemibrain Google Sheet
  hg = hemibrain_match_sheet(sheet = "hemibrain", selected_file = selected_file)
  orig = hg$flywire.xyz
  hg$flywire.xyz = fg$flywire.xyz[match(hg$bodyid,fg$hemibrain.match)]
  different = paste(orig)!=paste(hg$flywire.xyz)
  update_gsheet(update = hg[different,],
                  gs = hg,
                  tab = "hemibrain",
                  match = "flywire",
                  id = "bodyid")

  # Add missing flywire information
  missing = setdiff(flywire.ids, all.ids)
  hemibrain_matching_add(ids = missing, dataset="flywire", ...)
}

#' @rdname hemibrain_matching
#' @export
LR_matching <- function(ids = NULL,
                        mirror.nblast = NULL,
                        selected_file = options()$hemibrainr_matching_gsheet,
                        batch_size = 50,
                        db = flywire_neurons(),
                        query = flywire_neurons(mirror=TRUE),
                        overwrite = c("FALSE","mine","mine_empty","TRUE"),
                        column = NULL,
                        field = NULL){
  # Packages
  if(!requireNamespace("elmr", quietly = TRUE)) {
    stop("Please install elmr using:\n", call. = FALSE,
         "remotes::install_github('natverse/elmr')")
  }
  if(!requireNamespace("fafbseg", quietly = TRUE)) {
    stop("Please install fafbseg using:\n", call. = FALSE,
         "remotes::install_github('natverse/fafbseg')")
  }
  # Motivate!
  nat::nopen3d()
  plot_inspirobot()
  unsaved = c()
  message("
          #######################Colours##########################
          black = FAFB CATMAID neuron,
          dark grey = flywire neuron,
          blue = mirrrored flywire neuron you are trying to match,
          red = potential hemibrain matches based on NBLAST score,
          green = a chosen hemibrain neuron during scanning,
          dark blue = your selected hemibrain match.
          #######################Colours##########################
          ")
  ## Get NBLAST
  if(is.null(mirror.nblast)){
    message("Loading flywire NBLAST from flyconnectome Google Team Drive using Google Filestream: ")
    message(paste0(options()$Gdrive_hemibrain_data,"hemibrain_nblast/flywire.mirror.mean.rda"))
    mirror.nblast = hemibrain_nblast("flywire-mirror")
  }
  # Read the Google Sheet
  gs = hemibrain_match_sheet(selected_file = selected_file, sheet = "flywire")
  id = "flywire.id"
  # Get neuron data repo
  if(missing(db)) {
    db=tryCatch(force(db), error=function(e) {
      stop("Unable to use `flywire_neurons()`. ",
           "You must load the hemibrain Google Team Drive")
    })
  }
  if(missing(query)) {
    db=tryCatch(force(query), error=function(e) {
      stop("Unable to use `flywire_neurons(mirror=TRUE)`. ",
           "You must load the hemibrain Google Team Drive")
    })
  }
  # How much is done?
  match.field = paste0("FAFB.hemisphere",".match")
  quality.field = paste0("FAFB.hemisphere",".match.quality")
  done = subset(gs, !is.na(gs[[match.field]]))
  message("Neuron matches: ", nrow(done), "/", nrow(gs))
  print(table(gs[[quality.field]]))
  # choose user
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
  brain = elmr::FAFB.surf
  # Make matches!
  for(n in unique(selected[[id]])){
    # Get id
    n = as.character(n)
    end = n==selected[[id]][length(selected[[id]])]
    # Plot brain
    rgl::clear3d()
    hemibrain_view()
    rgl::bg3d("white")
    plot3d(brain, alpha = 0.1, col ="grey")
    # Transform hemibrain neuron to FAFB space
    fw.m = tryCatch(query[n], error = function(e){
      message("Could not read mirrored neuron: ", n)
      NULL
    })
    fw.n = tryCatch(db[n], error = function(e){
      message("Could not read native neuron: ", n)
      NULL
    })
    sk = gs[n,]$skid[1]
    if(!is.na(sk)){
      lhn = tryCatch(catmaid::read.neurons.catmaid(sk, OmitFailures = TRUE), error = function(e) NULL)
    }else{
      lhn = NULL
    }
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
    if(!is.null(lhn)){plot3d(lhn, lwd = 2, soma = TRUE, col = "black")}
    if(!is.null(fw.n)) {plot3d(fw.n, lwd = 3, soma = TRUE, col = "grey50")}
    if(!is.null(fw.m)) {plot3d(fw.m, lwd = 3, soma = TRUE, col = "#1BB6AF")}
    message("ID: ", n)
    message("ItoLee_Hemilineage : ",lhn[n,"ItoLee_Hemilineage"])
    # Read top 10 FIB matches
    r = tryCatch(sort(mirror.nblast[n,],decreasing = TRUE), error = function(e) NULL)
    if(is.null(r)){
      message(n, " not in NBLAST matrix, skipping ...")
      next
    }
    message(sprintf("Reading the top %s flywire hits",batch_size))
    # Read hemibrain neurons
    if(is.null(db)){
      native  = fafbseg::skeletor((names(r)[1:batch_size]), mesh3d = FALSE, clean = FALSE)
    } else {
      batch = names(r)[1:batch_size]
      batch.in = intersect(batch, names(db))
      native = tryCatch(db[match(batch.in,names(db))], error = function(e) NULL)
      if(is.null(native)|length(batch.in)!=length(batch)){
        # message("Cannot read neurons from local db; fetching from flywire!")
        # if(!requireNamespace("fafbseg", quietly = TRUE)) {
        #   stop("Please install fafbseg using:\n", call. = FALSE,
        #        "remotes::install_github('natverse/fafbseg')")
        # }
        # batch.out = setdiff(batch, names(native))
        # native2 =(tryCatch(fafbseg::skeletor(batch.out, mesh3d = TRUE, clean = FALSE), error=function(e) NULL) )
        # if(!is.null(native2)){
        #   native = nat::union(native, native2)
        #   native = native[as.character(batch)]
        # }
        message("Dropping ",length(batch)-length(batch.in) ," neuron missing from db" )
        batch = intersect(batch, batch.in)
      }
    }
    sel = c("go","for","it")
    k = 1
    j = batch_size
    # Cycle through potential matches
    while(length(sel)>1){
      sel = sel.orig = tryCatch(nat::nlscan(native[names(r)[1:j]], col = "#EE4244", lwd = 3, soma = TRUE),
                                error = function(e) NULL)
      if(is.null(sel)){
        next
      }
      if(length(sel)>1){
        message("Note: You selected more than one neuron")
      }
      if(length(sel) > 0){
        rgl::plot3d(native[sel], lwd = 2, soma = TRUE)
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
            native2  = fafbseg::skeletor((names(r)[(k+1):j]), mesh3d = TRUE, clean = FALSE)
          } else {
            native2 = tryCatch(db[(names(r)[(k+1):j])], error = function(e) {
              warning("Cannot read neuron: ", n, " from local db; fetching from flywire!")
              fafbseg::skeletor((names(r)[(k+1):j]), mesh3d = TRUE, clean = FALSE)
            })
          }
          native = nat::union(native, native2)
        }
      }else{
        while(length(sel)>1){
          message("Choose single best match: ")
          sel = nat::nlscan(native[as.character(sel.orig)], col = "orange", lwd = 2, soma = TRUE)
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
    if(length(sel)){
      sel = as.character(sel)
      if(!is.na(db[sel,"flywire.xyz"])){
        hit = db[sel,"flywire.xyz"]
      }else{
        fixed = flywire_basics(db[sel])
        hit = fixed[,"flywire.xyz"]
      }
    }else{
      hit = "none"
    }
    gs[gs[[id]]%in%n,match.field] = hit
    if(length(sel)){
      rgl::plot3d(native[sel],col="blue",lwd=2,soma=TRUE)
      quality = must_be("What is the quality of this match? good(e)/okay(o)/poor(p)/tract-only(t) ", answers = c("e","o","p","t"))
    }else{
      quality = "n"
    }
    quality = standardise_quality(quality)
    gs[gs[[id]]%in%n,quality.field] = quality
    unsaved = unique(c(unsaved, n))
    message(length(unsaved), " unsaved matches")
    print(knitr::kable(gs[unsaved,c(id,"ItoLee_Hemilineage",match.field,quality.field)]))
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

