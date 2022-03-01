##### Flywire-Hemibrain matches

# hidden
flytable_meta <- function(){
  ft=fafbseg::flytable_query("select _id, root_id, supervoxel_id, proofread, status, pos_x, pos_y, pos_z, nucleus_id, side, ito_lee_hemilineage, hartenstein_hemilineage, cell_class, cell_type, root_duplicated from info")
  ft$flywire_xyz = apply(ft[,c("pos_x", "pos_y", "pos_z")], 1, paste_coords)
  ft = as.data.frame(ft)
  ft
}

# hidden
flytable_matches <- function(dataset=c("flywire","hemibrain")){
  dataset = match.arg(dataset)
  ft = get_matching_flytable(dataset=dataset)
  ft = as.data.frame(ft)
  ft = ft[!duplicated(ft$id),]
  rownames(ft) = ft$id
  ft
}

# hidden
get_matching_flytable <- function(dataset=c("flywire","hemibrain")){
  dataset = match.arg(dataset)
  tab = "flywire_right" #switch(dataset, flywire = {"flywire_matching"}, hemibrain = {"hemibrain_matches"})
  ft=fafbseg::flytable_query(sprintf("select _id, root_id, proofread, status, flywire_xyz, cell_type,
cell_class, ito_lee_hemilineage, top_nt, hemisphere_match, hemisphere_match_xyz, hemisphere_match_quality, hemisphere_match_top_nt,
hemibrain_match, hemibrain_match_quality, hemibrain_match_top_nt, hemisphere_match_nblast, hemibrain_match_nblast, hemisphere_ngl_scene,
hemibrain_ngl_scene, user from %s", tab), base = "flywire_matching")
  if(dataset=="flywire"){
    ft$id = ft$root_id
  }else if(dataset=="hemibrain"){
    ft$id = ft$hemibrain_match
  }
  ft
}

# hidden
flytable_matches_update <- function(df, dataset=c("flywire","hemibrain"), DryRun = FALSE){
  ft=flytable_matches(dataset=dataset)
  ft$id = NULL
  tab = "flywire_right" #switch(dataset, flywire = {"matching"}, hemibrain = {"hemibrain_matches"})
  toupdate = setdiff(df, ft)
  if(DryRun)
    toupdate
  else {
    message(sprintf("updating %d matches to %s in %s", nrow(toupdate),  dataset, tab))
    fafbseg::flytable_update_rows(toupdate, table = tab)
    invisible(toupdate)
  }
}

# hidden
flytable_matches_update <- function(df = NULL, update_roots=TRUE, DryRun = FALSE, tab = "flywire_right"){

  # Get master data
  if(is.null(df)){
    df = flytable_meta()
  }
  df[is.na(df)] = "none"

  # Update IDs
  if(update_roots) {
    ft=fafbseg::flytable_query(sprintf("select _id, root_id, flywire_xyz, hemisphere_match, hemisphere_match_xyz from %s", tab))
    ft$root_id=fafbseg::flywire_updateids(ft$root_id, xyz = ft$flywire_xyz, rawcoords = TRUE)
    hemisphere.good = ! ft$hemisphere_match_xyz %in% c("NA","none","") & !is.na(ft$hemisphere_match_xyz)
    ft$hemisphere_match[hemisphere.good] = fafbseg::flywire_updateids(ft$hemisphere_match[hemisphere.good], xyz = ft$hemisphere_match_xyz[hemisphere.good], rawcoords = TRUE)
    df$root_id=fafbseg::flywire_updateids(df$root_id, svids = df$supervoxel_id)
    fafbseg::flytable_update_rows(ft, table = tab)
  }

  # Add missing rows
  if(update_roots){
    ft=fafbseg::flytable_query(sprintf("select _id, root_id, flywire_xyz from %s", tab))
    df.r = subset(df, !is.na(df$side) & df$side %in% c("right","center"))
    new = subset(df.r, ! df.r$root_id %in% ft$root_id)
    new = subset(new, ! new$flywire_xyz %in% ft$flywire_xyz)
    new$id = new$`_id` = NULL
    for(cl in colnames(ft)){
      if(!cl %in% colnames(new)){
        new[[cl]] = "none"
      }
    }
    fafbseg::flytable_append_rows(new, "flywire_right")
  }

  # Update match top_nt
  ft=fafbseg::flytable_query(sprintf("select _id, root_id, proofread, status, flywire_xyz, cell_type, cell_class, ito_lee_hemilineage, top_nt, hemisphere_match_top_nt from %s", tab))
  df$hemisphere_match = ft$hemisphere_match[match(df$root_id,ft$root_id)]
  df$hemisphere_match_top_nt = df$top_nt[match(df$hemisphere_match,df$root_id)]
  same.cols = intersect(colnames(ft), colnames(df))
  toupdate = dplyr::anti_join(ft[,same.cols], df[,same.cols], by = "root_id")

  # Write!
  if(DryRun)
    toupdate
  else if(nrow(toupdate)) {
    message(sprintf("updating %d rows for %s in %s", nrow(toupdate), "proofread, status, flywire_xyz, cell_type, cell_class, ito_lee_hemilineage, top_nt", tab))
    fafbseg::flytable_update_rows(toupdate, table = tab)
    invisible(toupdate)
  }
}

# hidden
flytable_update_proofread <- function(DryRun=FALSE, update_roots=TRUE) {
  proofread_st=fafbseg::flytable_query("select _id, root_id, supervoxel_id, proofread from info")
  if(update_roots)
    proofread_st$root_id=fafbseg::flywire_updateids(proofread_st$root_id, svids = proofread_st$supervoxel_id)
  psp=fafbseg::flywire_cave_query("proofreading_status_public_v1")
  toupdate=proofread_st %>%
    dplyr::mutate(proofread.new=root_id %in% as.character(psp$pt_root_id) & !root_id==0) %>%
    dplyr::filter(proofread.new!=proofread) %>%
    dplyr::select(`_id`, proofread.new) %>%
    dplyr::rename(proofread=proofread.new)
  if(DryRun)
    toupdate
  else if(nrow(toupdate)){
    message(sprintf("updating %d rows for %s in %s", nrow(toupdate), "proofread", "info"))
    fafbseg::flytable_update_rows(toupdate, table = "info")
    invisible(toupdate)
  }
}

# hidden
flytable_update_status <- function(DryRun=FALSE, tab = "info") {
  ft=fafbseg::flytable_query(sprintf("select _id, root_id, status from %s", tab))
  ft$status.new = standard_statuses(ft$status )
  toupdate=ft %>%
    dplyr::filter(status.new!=status) %>%
    dplyr::select(`_id`, status.new) %>%
    dplyr::rename(status=status.new)
  if(DryRun)
    toupdate
  else if(nrow(toupdate)){
    message(sprintf("updating %d rows for %s in %s", nrow(toupdate), "status", tab))
    fafbseg::flytable_update_rows(toupdate, table = tab)
    invisible(toupdate)
  }
}

# hidden
flytable_update_top_nt <- function( df = NULL, DryRun=FALSE, update_roots=T, tab = "info") {
  if(is.null(df)){
    df=flywire_meta(flytable = FALSE)
  }
  nt_st=fafbseg::flytable_query(sprintf("select _id, root_id, supervoxel_id, top_nt from %s", tab))
  nt_st[is.na(nt_st)] = ""
  if(update_roots)
    nt_st$root_id=fafbseg::flywire_updateids(nt_st$root_id, svids = nt_st$supervoxel_id)
    df$root_id=fafbseg::flywire_updateids(df$root_id, svids = df$flywire_xyz, rawcoords = TRUE)
  toupdate=nt_st %>%
    dplyr::mutate(top_nt.new = df$top_nt[match(root_id,df$root_id)]) %>%
    dplyr::filter(top_nt.new!=top_nt) %>%
    dplyr::select(`_id`, top_nt.new) %>%
    dplyr::rename(top_nt=top_nt.new)
  if(DryRun)
    toupdate
  else if(nrow(toupdate)){
    message(sprintf("updating %d rows for %s in %s", nrow(toupdate), "top_nt", tab))
    fafbseg::flytable_update_rows(toupdate, table = tab)
    invisible(toupdate)
  }
}

# hidden
flytable_update_side <- function(df = NULL, DryRun=FALSE, update_roots=T, tab = "info") {
  if(is.null(df)){
    df=flywire_meta(flytable = FALSE)
  }
  nt_st=fafbseg::flytable_query(sprintf("select _id, root_id, supervoxel_id, side from %s", tab))
  if(update_roots)
    nt_st$root_id=fafbseg::flywire_updateids(nt_st$root_id, svids = nt_st$supervoxel_id)
  toupdate=nt_st %>%
    dplyr::mutate(side.new = df$side[match(root_id,df$root_id)]) %>%
    dplyr::filter(side.new!=side) %>%
    dplyr::select(`_id`, side.new) %>%
    dplyr::rename(side=side.new)
  if(DryRun)
    toupdate
  else if(nrow(toupdate)){
    message(sprintf("updating %d rows for %s in %s", nrow(toupdate), "side", tab))
    fafbseg::flytable_update_rows(toupdate, table = tab)
    invisible(toupdate)
  }
}

# hidden
flytable_update_hl_predictions <- function(df = NULL, DryRun=FALSE, update_roots=TRUE, tab = "info") {
  if(is.null(df)){
    df=flywire_meta(flytable = FALSE)
  }
  nt_st=fafbseg::flytable_query(sprintf("select _id, root_id, supervoxel_id, hemilineage_nblast_1, hemilineage_nblast_2, hemilineage_nblast_3 from %s", tab))
  if(update_roots)
    nt_st$root_id=fafbseg::flywire_updateids(nt_st$root_id, svids = nt_st$supervoxel_id)
  toupdate=nt_st %>%
    dplyr::mutate(hemilineage_nblast_1.new = df$hemilineage_nblast_1[match(root_id,df$root_id)]) %>%
    dplyr::mutate(hemilineage_nblast_2.new = df$hemilineage_nblast_2[match(root_id,df$root_id)]) %>%
    dplyr::mutate(hemilineage_nblast_3.new = df$hemilineage_nblast_3[match(root_id,df$root_id)]) %>%
    dplyr::filter(hemilineage_nblast_1.new!=hemilineage_nblast_1) %>%
    dplyr::filter(hemilineage_nblast_2.new!=hemilineage_nblast_2) %>%
    dplyr::filter(hemilineage_nblast_3.new!=hemilineage_nblast_3) %>%
    dplyr::select(`_id`, hemilineage_nblast_1.new, hemilineage_nblast_2.new, hemilineage_nblast_3.new) %>%
    dplyr::rename(hemilineage_nblast_1=hemilineage_nblast_1.new) %>%
    dplyr::rename(hemilineage_nblast_2=hemilineage_nblast_2.new) %>%
    dplyr::rename(hemilineage_nblast_3=hemilineage_nblast_3.new)
  if(DryRun)
    toupdate
  else if(nrow(toupdate)){
    message(sprintf("updating %d rows for %s in %s", nrow(toupdate), "hemilineage_nblast", tab))
    fafbseg::flytable_update_rows(toupdate, table = tab)
    invisible(toupdate)
  }
}

