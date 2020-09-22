# Flywire neuron matching

#' @rdname hemibrain_add_made_matches
#' @export
flywire_matching_rewrite <- function(flywire.ids = names(flywire_neurons()),
                                     selected_file  = "1OSlDtnR3B1LiB5cwI5x5Ql6LkZd8JOS5bBr-HTi0pOw",
                                     ...){
  # Get the FAFB matching google sheet
  gs = gsheet_manipulation(FUN = googlesheets4::read_sheet,
                           ss = selected_file,
                           sheet = "FAFB",
                           guess_max = 3000,
                           return = TRUE)
  skids = unique(gs$skid)
  skids = as.integer(skids)
  skids = skids[!is.na(skids)]
  rownames(gs) = correct_id(gs$skid)

  # Get FAFBv14 coordinates
  cats = nat::neuronlist()
  batches = split(1:length(skids), round(seq(from = 1, to = 100, length.out = length(skids))))
  all.ids = c()
  for(i in 1:10){
    # Read CATMAID neurons
    message("Batch:", i, "/100")
    cat = catmaid::read.neurons.catmaid(skids[batches[[i]]], OmitFailures = TRUE)
    cats = c(cats,cat)

    # Get xyz for root points
    roots = sapply(cats, function(y) nat::xyzmatrix(y)[nat::rootpoints(y),])
    roots = t(roots)
    FAFB.xyz = apply(roots, 1, paste, collapse = ",")

    # Get Flywire voxel coordinates
    roots.flywire = nat.templatebrains::xform_brain(roots, reference = "FlyWire", sample = "FAFB14", .parallel = TRUE, verbose = TRUE)
    rownames(roots.flywire) = rownames(roots)
    roots.flywire.raw = scale(roots.flywire, scale = c(4, 4, 40), center = FALSE)
    fw.ids = fafbseg::flywire_xyz2id(roots.flywire.raw, rawcoords = TRUE)
    fw.ids[fw.ids=="0"] = NA
    flywire.xyz = apply(roots.flywire.raw, 1, paste, collapse = ",")

    # Add
    g = gs[rownames(roots),c("FAFB.xyz","flywire.xyz","flywire.id","skid")]
    g[,]$FAFB.xyz = FAFB.xyz
    g[,]$flywire.xyz = flywire.xyz
    g[,]$flywire.id = fw.ids
    rownames(g) = NULL

    # Update
    if(nrow(g)){
      for(row in 1:nrow(g)){
        columns = c("FAFB.xyz","flywire.xyz","flywire.id")
        r = match(g[row,"skid"],gs[["skid"]])+1
        if(is.issue(r)){
          next
        }
        for(column in columns){
          letter = LETTERS[match(column,colnames(gs))]
          range = paste0(letter,r)
          gsheet_manipulation(FUN = googlesheets4::range_write,
                              ss = selected_file,
                              range = range,
                              data = as.data.frame(g[row,column]),
                              sheet = "FAFB",
                              col_names = FALSE)
        }
      }
    }
    all.ids=unique(c(all.ids,fw.ids))
  }

  # Add flywire match to the hemibrain and LM sheets
  ## Read the FAFB Google Sheet
  fg = hemibrain_match_sheet(sheet = "FAFB", selected_file = selected_file)

  ## Read the LM Google Sheet
  lmg = hemibrain_match_sheet(sheet = "lm", selected_file = selected_file)
  orig = lmg$flywire.match
  lmg$flywire.match = fg$flywire.xyz[match(lmg$bodyid,fg$hemibrain.match)]
  different = paste(orig)!=paste(lmg$flywire.match)
  lmg$flywire.match = fg$flywire.xyz[match(lmg$id,fg$LM.match)]
  update_gsheet(update = lmg[different,],
                  gs = lmg,
                  tab = "LM",
                  match = "flywire",
                  id = "id")

  ## Read the hemibrain Google Sheet
  hg = hemibrain_match_sheet(sheet = "hemibrain", selected_file = selected_file)
  orig = hg$flywire.match
  hg$flywire.match = fg$flywire.xyz[match(hg$bodyid,fg$hemibrain.match)]
  different = paste(orig)!=paste(hg$flywire.match)
  update_gsheet(update = hg[different,],
                  gs = hg,
                  tab = "hemibrain",
                  match = "flywire",
                  id = "bodyid")

  # Add missing flywire information
  missing = setdiff(flywire.ids, all.ids)
  hemibrain_matching_add(ids = missing, dataset="flywire", ...)
}


