# FlyWire neuron matching

#' @rdname hemibrain_add_made_matches
#' @export
flywire_matching_gsheet_rewrite <- function(flywire_ids = names(flywire_neurons()),
                                     meta = flywire_neurons()[,],
                                     catmaid.update = TRUE,
                                     selected_file  = options()$hemibrainr_matching_gsheet, # 1_RXfVRw2nVjzk6yXOKiOr_JKwqk2dGSVG_Xe7_v8roU
                                     reorder = FALSE,
                                     top.nblast = FALSE,
                                     nblast = NULL,
                                     ...){
  # Get the FAFB matching Google sheet
  gs = hemibrain_match_sheet(sheet = "FAFB", selected_file = selected_file)
  skids = as.character(unique(gs$skid))
  skids = id_okay(skids)

  # Get FAFBv14 coordinates
  if(length(skids) & catmaid.update){
    cats = nat::neuronlist()
    batches = split(1:length(skids), round(seq(from = 1, to = 100, length.out = length(skids))))
    all.ids = c()
    for(i in 1:length(batches)){
      # Read CATMAID neurons
      message("Batch:", i, "/10")
      search = skids[batches[[i]]]
      cat = tryCatch(catmaid::read.neurons.catmaid(search, OmitFailures = TRUE), error = function(e) NULL)
      if(!length(cat)){
        next
      }
      cats = nat::union(cats,cat)

      # Get xyz for primary branch points
      simp = nat::nlapply(cat,nat::simplify_neuron,n=1, .parallel = TRUE, OmitFailures = TRUE)
      branchpoints = sapply(simp, function(y) nat::xyzmatrix(y)[ifelse(length(nat::branchpoints(y)),nat::branchpoints(y),max(nat::endpoints(y))),])
      branchpoints = t(branchpoints)
      fafb_xyz = apply(branchpoints, 1, paste_coords)

      # Get FlyWire voxel coordinates
      branchpoints.flywire = nat.templatebrains::xform_brain(branchpoints, reference = "FlyWire", sample = "FAFB14", .parallel = TRUE, verbose = TRUE)
      rownames(branchpoints.flywire) = rownames(branchpoints)
      branchpoints.flywire.raw = scale(branchpoints.flywire, scale = c(4, 4, 40), center = FALSE)
      fw.ids = fafbseg::flywire_xyz2id(branchpoints.flywire.raw, rawcoords = TRUE)
      flywire_xyz = apply(branchpoints.flywire.raw, 1, paste_coords)

      # Add
      indices = match(names(fafb_xyz),gs$skid)
      if(length(indices)){
        gs[indices,]$fafb_xyz = fafb_xyz
        gs[indices,]$flywire_xyz = flywire_xyz
        gs[indices,]$root_id = fw.ids
      }
    }
  }

  # Top NBLAST
  if(top.nblast){
    if(is.null(nblast)){
      nblast = tryCatch(hemibrain_nblast('hemibrain-flywire'), error = function(e) NULL)
    }
    if(!is.null(nblast)){
      nblast.top = nblast[match(gs$root_id,rownames(nblast)),]
      tops = apply(nblast.top,1,function(r) which.max(r))
      top = colnames(nblast)[unlist(tops)]
      top[!gs$root_id%in%rownames(nblast)] = NA
      gs$hemibrain.nblast.top = top
    }
  }

  # flywire_svids
  svids = fafbseg::flywire_xyz2id(nat::xyzmatrix(gs$flywire_xyz), root=FALSE, rawcoords = TRUE)
  gs[,]$flywire_svid = svids

  # Update fafb_xyz column
  empty = is.na(gs$fafb_xyz) & ! is.na(gs$flywire_xyz)
  if(sum(empty)){
    fafb_xyz = meta[gs[empty,"root_id"],"fafb_xyz"]
    gs[empty,"fafb_xyz"] = fafb_xyz
  }
  empty = is.na(gs$fafb_xyz) & ! is.na(gs$flywire_xyz)
  if(sum(empty)){
    points.raw = nat::xyzmatrix(gs[empty,"flywire_xyz"])
    points.nm = scale(points.raw, scale = c(4, 4, 40), center = FALSE)
    points.fafb = nat.templatebrains::xform_brain(points.nm, sample = "FlyWire", reference = "FAFB14", .parallel = TRUE, verbose = TRUE)
    fafb_xyz = apply(points.fafb, 1, paste_coords)
    gs[empty,"fafb_xyz"] = fafb_xyz
  }

  # Update side information
  if(!is.null(meta$side)){
    # Update side
    sides = meta[match(gs$root_id,meta$root_id),"side"]
    sides[is.na(sides)] = gs$side[is.na(sides)]
    gs$side = sides
  }

  # Update
  write.cols = intersect(c("fafb_xyz", "flywire_xyz", "root_id", "flywire_svid", "side", "nblast.top"),colnames(gs))
  gsheet_update_cols(
      write.cols = write.cols,
      gs=gs,
      selected_sheet = selected_file,
      sheet = "FAFB")

  # Figure out duplicate entries
  fg = hemibrain_match_sheet(sheet = "FAFB", selected_file = selected_file)
  fg$index = 1:nrow(fg)+1
  removals = data.frame()
  for(set in c('skid',"flywire_xyz","flywire_svid")){
    dupes = unique(fg[[set]][duplicated(fg[[set]])])
    dupes = id_okay(dupes)
    for(dupe in dupes){
      many = fg[[set]] == dupe
      many[is.na(many)] = FALSE
      sub = fg[many,]
      skd = unique(sub$skid)
      skd = id_okay(skd)
      if(length(skd)>1){
        next
      }
      best = which.max(apply(sub, 1, function(r) sum(!is.na(r[c("hemibrain_match", "hemibrain_match_quality",
                                                                "LM.match", "LM.match.quality",
                                                                "fafb_hemisphere_match", "fafb_hemisphere_match.quality")]))))
      remove = sub[-best,]
      removals = rbind(removals, remove)
    }
  }
  if(reorder){
    if(nrow(removals)){
      n = fg[!fg$index%in%removals$index,]
    }
    n$index = NULL
    n = n[!is.na(n$flywire_xyz)|!is.na(n$skid),]
    gsheet_reorder(gs=n,tab="FAFB",selected_sheet=selected_file,field = "flywire_xyz", remove.duplicates = FALSE)
  }else if (nrow(removals)){
    for(r in sort(removals$index,decreasing = TRUE)){
      range.del = googlesheets4::cell_rows(r)
      message("Removing a row for: ", dupe)
      gsheet_manipulation(FUN = googlesheets4::range_delete,
                          ss = selected_file,
                          range = range.del,
                          sheet = "FAFB")
    }
  }

  # Add missing flywire information
  fg = hemibrain_match_sheet(sheet = "FAFB", selected_file = selected_file)
  all.ids = correct_id(unique(fg$root_id))
  missing = setdiff(flywire_ids, all.ids)
  if(length(missing)){
    hemibrain_matching_add(ids = missing, meta = meta, dataset="flywire", selected_file = selected_file, ...)
  }

  # ## Read the LM Google Sheet
  # lmg = hemibrain_match_sheet(sheet = "lm", selected_file = selected_file)
  # if(nrow(lmg)){
  #   lmg$flywire_xyz = fg$flywire_xyz[match(lmg$id,fg$LM.match)]
  #   gsheet_update_cols(
  #     write.cols = "flywire_xyz",
  #     gs=lmg,
  #     selected_sheet = selected_file,
  #     sheet = "lm")
  # }
  #
  # ## Read the hemibrain Google Sheet
  # hg = hemibrain_match_sheet(sheet = "hemibrain", selected_file = selected_file)
  # if(nrow(hg)){
  #   hg$flywire_xyz = fg$flywire_xyz[match(hg$bodyid,fg$hemibrain_match)]
  #   gsheet_update_cols(
  #     write.cols = "flywire_xyz",
  #     gs=hg,
  #     selected_sheet = selected_file,
  #     sheet = "hemibrain")
  # }

}

#' @rdname hemibrain_matching
#' @export
LR_matching <- function(ids = NULL,
                        threshold = 0,
                        mirror.nblast = NULL,
                        flytable = TRUE,
                        selected_file = options()$hemibrainr_matching_gsheet,
                        batch_size = 50,
                        db = flywire_neurons(),
                        query = flywire_neurons(mirror=TRUE),
                        overwrite = c("FALSE","bad","TRUE","review"),
                        column = NULL,
                        entry = NULL,
                        user = NULL,
                        superuser = FALSE,
                        flywire.good = FALSE,
                        verbose = FALSE){
  message("Matching mirrored flywire neurons (blue) to non-mirrored flywire neurons (red)")
  mirror = TRUE
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
  unsaved = saved = c()
  message("
          #######################Colours##########################
          black = FAFB CATMAID neuron,
          dark grey = flywire neuron,
          blue = mirrrored flywire neuron you are trying to match,
          red = potential hemibrain_matches based on NBLAST score,
          green = a chosen hemibrain neuron during scanning,
          dark blue = your selected hemibrain_match.
          #######################Colours##########################
          ")
  ## Get NBLAST
  if(is.null(mirror.nblast)){
    message("Loading flywire NBLAST from flyconnectome Google Team Drive using Google Filestream: ")
    message(file.path(options()$remote_connectome_data,"hemibrain_nblast/flywire.mirror.mean.rda"))
    mirror.nblast = hemibrain_nblast("flywire-mirror")
  }
  # Read the Google Sheet or seatable
  if(flytable){
    gs = flytable_matches(dataset="flywire")
  }else{
    gs = hemibrain_match_sheet(selected_file = selected_file, sheet = "flywire")
  }
  id = "root_id"
  # Get neuron data repo
  if(missing(db)) {
    db=tryCatch(force(db), error=function(e) {
      stop("Unable to use `flywire_neurons()`. ",
           "You must load the hemibrain Google Team Drive")
    })
  }
  if(missing(query)) {
    query=tryCatch(force(query), error=function(e) {
      stop("Unable to use `flywire_neurons(mirror=TRUE)`. ",
           "You must load the hemibrain Google Team Drive")
    })
  }
  # How much is done?
  match.field = "hemisphere_match_xyz"
  quality.field = "hemisphere_match_quality"
  done = subset(gs, !is.na(gs[[match.field]]) & !gs[[match.field]]%in%c("none","unknown","NA","na"))
  message("Neuron matches: ", nrow(done), "/", nrow(gs))
  print(table(gs[[quality.field]]))
  # choose user
  initials = choose_user(gs, user = user)
  # choose ids
  selected = id_selector(gs=gs, ids=ids, id=id, overwrite = overwrite,
                         quality.field = quality.field, match.field = match.field,
                         initials = initials, column = column, entry = entry, superuser = superuser,
                         flywire.good = flywire.good)
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
  # choose brain
  brain = elmr::FAFB.surf
  # Make matches!
  match.more = TRUE
  while(match.more){
    match_cycle = neuron_match_scanner(brain = brain,
                                       selected = selected,
                                       id = id,
                                       unsaved = unsaved,
                                       saved = saved,
                                       chosen.field = "flywire_xyz",
                                       nblast = mirror.nblast,
                                       threshold = threshold,
                                       batch_size = batch_size,
                                       targets = db,
                                       targets.repository = "flywire",
                                       query = query,
                                       extra.neurons = db,
                                       query.repository = "flywire",
                                       extra.repository = "CATMAID",
                                       match.field = match.field,
                                       quality.field = quality.field,
                                       soma.size = 4000,
                                       show.columns = c("cell_type","ito_lee_hemilineage","status", match.field, quality.field,"note"),
                                       skip.if.absent = !verbose,
                                       mirror.query = mirror)
    if(is.null(match_cycle)){
      next
    }
    selected = match_cycle[["selected"]]
    unsaved = match_cycle[["unsaved"]]
    selected$user = initials
    if(length(unsaved)){
      plot_inspirobot()
      say_encouragement(initials)
      if(flytable){ # seatable
        flytable_matches_update(selected)
      }else{ # google sheet
        # Read!
        gs2 = hemibrain_match_sheet(selected_file = selected_file, sheet = "flywire")
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
      }
      saved = c(unsaved, saved)
      unsaved = c()
    }
    match.more = hemibrain_choice("Match more neurons? ")
  }
  say_encouragement(initials)
}

# Hidden
neuron_match_scanner <- function(brain,
                                 selected,
                                 id,
                                 chosen.field,
                                 nblast,
                                 threshold,
                                 batch_size,
                                 targets,
                                 targets.repository = c("CATMAID","flywire","hemibrain","lm"),
                                 query,
                                 query.repository = c("CATMAID","flywire","hemibrain","lm"),
                                 extra.neurons = NULL,
                                 extra.repository = c("none","CATMAID","flywire","hemibrain","lm"),
                                 match.field,
                                 quality.field,
                                 unsaved = c(),
                                 saved = c(),
                                 soma.size = 4000,
                                 show.columns = c("cell_type","ito_lee_hemilineage","status",  match.field, quality.field,"note"),
                                 skip.if.absent = TRUE,
                                 mirror.query = FALSE){
  targets.repository = match.arg(targets.repository)
  extra.repository = match.arg(extra.repository)
  query.repository = match.arg(query.repository)
  check = unique(setdiff(unique(selected[[id]]),c(unsaved,saved)))
  check = id_okay(check)
  message("We'll look at ", length(check)," ", query.repository, " neurons sequentially.")
  for(i in 1:length(check)){
    # Get id
    message("neuron query ",i+length(saved),"/",length(check)+length(saved))
    n = as.character(check[i])
    end = n==selected[[id]][length(selected[[id]])]
    # Read top 10  matches
    r = tryCatch(sort(nblast[n,],decreasing = TRUE), error = function(e) NULL)
    if(is.null(r)||!ncol(r)){
      message(n, " not in NBLAST matrix, skipping ...")
      if(!skip.if.absent){
        progress = readline(prompt = "Press any key to continue ")
      }
      next
    }
    if(!is.null(threshold)){
      res = names(r)[r>threshold]
      if(!length(res)&!skip.if.absent){
        message(" no normalised NBLAST score greater or equal to ", threshold," for neuron ", n," ...")
        rgl::clear3d()
        rgl::rgl.viewpoint(userMatrix = structure(c(0.990777730941772, 0.049733679741621,
                                                    -0.126039981842041, 0, 0.060652956366539, -0.994590044021606,
                                                    0.084330290555954, 0, -0.121164083480835, -0.091197244822979,
                                                    -0.988434314727783, 0, 0, 0, 0, 1), .Dim = c(4L, 4L)), zoom = 0.644609212875366) # FAFB14 view
        rgl::bg3d("white")
        plot3d(brain, alpha = 0.1, col ="grey")
        query.n = get_match_neuron(query = query, n = n, query.repository = query.repository, skip.if.absent = skip.if.absent, mirror = mirror.query)
        if(!is.null(query.n)&&length(query.n)){plot3d(query.n, lwd = 3, soma = soma.size, col = "#1BB6AF")}
        progress = readline(prompt = "This neuron will be skipped. Press any key to continue ")
        next
      }else if(!length(res)){
        message(" no normalised NBLAST score greater or equal to ", threshold," for neuron ", n," ...")
        next
      }
    }
    batch.size = ifelse(length(res)>=batch_size,batch_size, length(res))
    # Plot brain
    rgl::clear3d()
    rgl::rgl.viewpoint(userMatrix = structure(c(0.990777730941772, 0.049733679741621,
                             -0.126039981842041, 0, 0.060652956366539, -0.994590044021606,
                             0.084330290555954, 0, -0.121164083480835, -0.091197244822979,
                             -0.988434314727783, 0, 0, 0, 0, 1), .Dim = c(4L, 4L)), zoom = 0.644609212875366) # FAFB14 view
    rgl::bg3d("white")
    plot3d(brain, alpha = 0.1, col ="grey")
    # Get data
    query.n = get_match_neuron(query = query, n = n, query.repository = query.repository, mirror = mirror.query)
    if(is.null(query.n)||!length(query.n)){
      message("Could not find query neuron: ", n)
      return(NULL)
    }
    # Other neurons to always plot
    if(extra.repository!="none"){
      if(extra.repository=="CATMAID" & "skid" %in% colnames(selected)){
        sk = selected[n,]$skid[1]
        extra.n = try(get_match_neuron(query = NULL, n = sk, query.repository = extra.repository), silent = FALSE)
      }else{
        fw.id = selected[n,]$root_id[1]
        extra.n = try(get_match_neuron(query = flywire_neurons(), n = fw.id, query.repository = extra.repository), silent = FALSE)
      }
    }else{
      extra.n = NULL
    }
    if(!is.null(extra.neurons)){
      another.n = tryCatch(extra.neurons[n], error = function(e){
        message(e)
        NULL
      })
      if(is.null(another.n)){
        fw.id = selected[n,]$root_id[1]
        another.n = get_match_neuron(query = flywire_neurons(), n = fw.id, query.repository = query.repository)
      }
    }else{
      another.n = NULL
    }
    ### Get old matches
    if(match.field=="flywire_xyz"){
      old.match = selected[selected[[id]]%in%n,"root_id"][1]
    }else{
      old.match = selected[selected[[id]]%in%n,match.field][1]
    }
    old.quality = selected[selected[[id]]%in%n,quality.field][1]
    if(!is.na(old.match)&&!old.match%in%c(""," ")){
      match.n =  tryCatch(targets[old.match], error = function(e){
        message("Could not immediately load match neuron: ", old.match)
        message(e)
        NULL})
    }else{
      match.n = NULL
    }
    ### Plot in 3D
    if(!length(extra.n)){
      extra.n=NULL
    }
    if(!length(another.n)){
      extra.n=NULL
    }
    if(!is.null(query.n)){plot3d(query.n, lwd = 3, soma = soma.size, col = "#1BB6AF")}
    if(!is.null(extra.n)&&length(extra.n)){plot3d(extra.n, lwd = 2, soma = soma.size, col = "black")}
    if(!is.null(another.n)&&length(another.n)){plot3d(another.n, lwd = 3, soma = soma.size, col = "grey50")}
    if(!is.null(match.n)&&length(match.n)){plot3d(match.n, lwd = 2, soma = soma.size, col = "#348E53")}
    message("ID: ", n)
    show.columns = intersect(show.columns,colnames(query.n[,]))
    display = c()
    for(sc in show.columns){
      if(sc%in%colnames(selected)){
        display = c(display,paste0(sc," : ", selected[match(n,selected[[id]]),sc]))
      }else{
        display = c(display,paste0(sc," : ", query.n[n,sc]))
      }
    }
    message(paste(display,collapse="  |  "))
    # Read database neurons
    message(sprintf("Reading the top %s %s hits",batch.size, targets.repository))
    batch = res[1:batch.size]
    if(is.null(targets)){
      if(targets.repository=="flywire"){
        fafbseg::choose_segmentation("flywire")
        native  = fafbseg::skeletor(batch, mesh3d = FALSE, clean = FALSE)
        native[,"root_id"] = names(native)
      }else if (targets.repository == "hemibrain"){
        native  = neuprintr::neuprint_read_neurons(batch, all_segments = TRUE, heal = FALSE)
        native = scale_neurons.neuronlist(native, scaling = (8/1000))
        native = suppressWarnings(nat.templatebrains::xform_brain(native, reference = "FAFB14", sample = "JRCFIB2018F"))
      }else if (targets.repository == "CATMAID"){
        native  = catmaid::read.neurons.catmaid(batch,  OmitFailures = TRUE)
      }
    } else {
      batch.in = intersect(batch, names(targets))
      native = tryCatch(targets[match(batch.in,names(targets))], error = function(e){
        try(file.remove(paste0(attributes(targets)$db@datafile,"___LOCK")), silent = TRUE)
        message(e)
        NULL
      })
      if(is.null(native)|length(batch.in)!=length(batch)){
        message("Dropping ",length(batch)-length(batch.in) ," neuron missing from targets" )
        batch = intersect(batch, batch.in)
      }
    }
    sel = c("go","for","it")
    k = 1
    j = batch.size
    # Cycle through potential matches
    while(length(sel)>1){
      plot.order = match(res[1:j],names(native))
      plot.order = plot.order[!is.na(plot.order)]
      sel = sel.orig = tryCatch(nat::nlscan(native[plot.order], col = "#EE4244", lwd = 3, soma = soma.size),
                                error = function(e){
                                  message(e)
                                  NULL
                                })
      if(!is.issue(old.match)){
        sel = union(old.match,sel)
        sel = sel[!sel%in%c("none","NA","na","unknown")]
      }
      if(length(sel)>1){
        message("Note: You selected more than one neuron")
      }
      sel = sel[sel!=0]
      if(length(sel) > 0){
        missing = setdiff(sel,names(native))
        if(length(missing)&&missing!="0"){
          missed = get_match_neuron(query = NULL, n = missing, query.repository = targets.repository)
          if(is.null(missed)){
            message("Could not fetch selection: ", paste0(missing,collapse=", "))
            sel = NULL
          }else{
            names(missed) = missing
            native = nat::union(missed, native)
          }
        }
      }
      if(length(sel)){
        rgl::plot3d(native[sel], lwd = 2, soma = soma.size, col = hemibrain_colour_ramp(length(sel)))
        prog = hemibrain_choice(sprintf("You selected %s neurons. Are you happy with that? ",length(sel)))
        if(length(sel)>0){
          nat::npop3d()
        }
        if(!prog){
          sel = c("go","for","it")
          if(batch.size < length(r)){
            prog = hemibrain_choice(sprintf("Do you want to read %s more neurons? ", batch.size))
            if(prog){
              k = j
              j = j + batch_size
              if(!is.null(targets)){
                native2 = tryCatch(targets[(res[(k+1):j])], error = function(e) {
                  message("Cannot read neuron: ", n, " from local targets, fetching from remote!")
                  message(e)
                  NULL
                })
              }
              if(is.null(targets)|is.null(native2)){
                if(targets.repository=="flywire"){
                  fafbseg::choose_segmentation("flywire")
                  native2  = fafbseg::skeletor((res[1:batch.size]), mesh3d = FALSE, clean = FALSE)
                }else if (targets.repository == "hemibrain"){
                  native2  = neuprintr::neuprint_read_neurons((res[1:batch.size]), all_segments = TRUE, heal = FALSE)
                  native = scale_neurons.neuronlist(native2, scaling = (8/1000))
                  native2 = suppressWarnings(nat.templatebrains::xform_brain(native2, reference = "FAFB14", sample = "JRCFIB2018F"))
                }else if (targets.repository == "CATMAID"){
                  native2  = catmaid::read.neurons.catmaid((res[1:batch.size]), .progress = 'text', OmitFailures = TRUE)
                }
              }
              native = nat::union(native, native2)
            }
          }
        }else{
          while(length(sel)>1){
            message("Choose single best match: ")
            plot.order = match(sel.orig,names(native))
            plot.order = plot.order[!is.na(plot.order)]
            sel = nat::nlscan(native[plot.order], col = hemibrain_bright_colours["orange"], lwd = 2, soma = TRUE)
            message(sprintf("You selected %s neurons", length(sel)))
            if(!length(sel)){
              noselection = hemibrain_choice("You selected no neurons. Are you happy with that? ")
              if(!noselection){
                sel = sel.orig
              }
            }
          }
        }
      }else{
        sel = NULL
      }
    }
    # Assign match and its quality
    if(length(sel)){
      sel = as.character(sel)
      if(!is.na(native[sel,chosen.field])){
        hit = as.character(native[sel,chosen.field])
      }else if(targets.repository=="flywire"){
        theone = native[sel]
        if("mesh3d"%in%class(theone[[1]])){
          theone = fafbseg::skeletor(sel)
        }
        fixed = flywire_basics(theone)
        hit = as.character(fixed[,chosen.field])
      }else{
        hit = names(targets[sel])
      }
    }else{
      hit = selected[selected[[id]]%in%n,match.field]
      hit[is.issue(hit)] = "none"
      quality = old.quality
      quality[is.issue(quality)] = "none"
    }
    message("You chose: ", hit)
    selected[selected[[id]]%in%n,match.field] = hit
    if(match.field%in%c("fafb_hemisphere_match","flywire.match","FAFB.match","flywire_xyz")){
      try({selected[selected[[id]]%in%n,"flywire.match.id"] = fafbseg::flywire_xyz2id(selected[selected[[id]]%in%n,match.field], rawcoords = TRUE)}, silent = FALSE)
      show.columns = unique(c(show.columns,"flywire.match.id"))
    }
    if(length(sel)){
      rgl::plot3d(native[sel],col= hemibrain_bright_colours["navy"],lwd=2,soma=TRUE)
      quality = must_be("What is the quality of this match? good(e)/okay(o)/poor(p)/tract-only(t) ", answers = c("e","o","p","t"))
    }else{
      quality = "n"
    }
    quality = standardise_quality(quality)
    selected[selected[[id]]%in%n,quality.field] = quality
    message("At quality: ", quality)
    # Make a note?
    orig.note = selected[selected[[id]]%in%n,'note']
    if(!is.issue(orig.note)){
      message("This neuron has note: ", orig.note)
    }
    make.note = hemibrain_choice("Would you like to record a note? y/n ")
    while(make.note){
      note = readline(prompt = "Your note on this match/these neurons:  ")
      message(note)
      note[note%in%c(" ","","NA")] = NA
      selected[selected[[id]]%in%n,'note'] = note
      make.note = !hemibrain_choice("Happy with this note? y/n ")
      show.columns = unique(c(show.columns,"note"))
    }
    unsaved = unique(c(unsaved, n))
    message(length(unsaved), " unsaved matches")
    show.columns = intersect(show.columns,colnames(selected))
    print(knitr::kable(selected[selected[[id]]%in%unsaved,c(id,show.columns,match.field,quality.field)]))
    p = must_be("Continue (enter) or save (s)? ", answers = c("","s"))
    if(p=="s"|end){
      break
    }
  }
  list(selected = selected, unsaved = unsaved)
}

# hidden
get_match_neuron <- function(query = NULL, n, query.repository, skip.if.absent = FALSE, mirror = FALSE){
  if(!is.null(query)){
    query.n = tryCatch(query[n], error = function(e){
      message("Could not immediately load query neuron: ", n)
      try(file.remove(paste0(attributes(query)$db@datafile,"___LOCK")), silent = TRUE)
      message(e)
      NULL
    })
  }else{
    query.n = NULL
  }
  if(is.null(query.n)||!length(query.n)&&!skip.if.absent){ # in FAFB14 space.
    query.n = tryCatch({
      message("Neuron not found on Google drive, attempting to read ...")
      if(!requireNamespace("fafbseg", quietly = TRUE)) {
        stop("Please install fafbseg using:\n", call. = FALSE,
             "remotes::install_github('natverse/fafbseg')")
      }
      if(query.repository == "flywire"){
        fafbseg::choose_segmentation("flywire")
        # query.n = fafbseg::skeletor(n)
        query.n = fafbseg::read_cloudvolume_meshes(n)
      }else if(query.repository == "hemibrain"){
        query.n  = neuprintr::neuprint_read_neurons(n, all_segments = TRUE, heal = FALSE)
        query.n = scale_neurons.neuronlist(query.n, scaling = (8/1000))
        query.n = suppressWarnings(nat.templatebrains::xform_brain(query.n, reference = "FAFB14", sample = "JRCFIB2018F"))
      }else if (query.repository=="CATMAID"){
        query.n = catmaid::read.neurons.catmaid(n)
      }else{
        NULL
      }
      if(mirror){
        if(!requireNamespace("elmr", quietly = TRUE)) {
          stop("Please install fafbseg using:\n", call. = FALSE,
               "remotes::install_github('natverse/elmr')")
        }
        if(!requireNamespace("nat.templatebrains", quietly = TRUE)) {
          stop("Please install fafbseg using:\n", call. = FALSE,
               "remotes::install_github('natverse/nat.templatebrains')")
        }
        if(!requireNamespace("nat.flybrains", quietly = TRUE)) {
          stop("Please install fafbseg using:\n", call. = FALSE,
               "remotes::install_github('natverse/nat.flybrains')")
        }
        t = java_xform_brain(query.n, reference = "JRC2018F", sample = "FAFB14", .parallel = FALSE, verbose = TRUE, OmitFailures = TRUE, progress.rjava=TRUE)
        m = nat.templatebrains::mirror_brain(x = t, brain = nat.flybrains::JRC2018F, .parallel = TRUE, OmitFailures = FALSE, transform = "flip")
        query.n = java_xform_brain(m, reference = "FAFB14", sample = "JRC2018F",.parallel = TRUE, verbose = TRUE, OmitFailures = TRUE, progress.rjava=TRUE)
      }
    }, error = function(e) {NULL})
  }
  query.n
}

