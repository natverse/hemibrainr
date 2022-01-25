##### Flywire-Hemibrain matches

flytable_meta <- function(){
  ft=fafbseg::flytable_query("select _id, root_id, supervoxel_id,
                             proofread, status, pos_x, pos_y, pos_z, nucleus_id,
                             side, top_nt, top_p, ito_lee_hemilineage, hartenstein_hemilineage,
                             cell_type, hemibrain_match, root_duplicated from info")
  ft = as.data.frame(ft)
  #rownames(ft) = ft$root_id
  ft
}

get_matching_flytable <- function(dataset=c("flywire","hemibrain")){
  tab = switch(dataset, flywire = {"matching"}, hemibrain = {"hemibrain_matches"})
  ft=fafbseg::flytable_query(sprintf("select _id, id, side, cell, cell_type,cell_body_fiber,
                             ito_lee_hemilineage, match, quality, fafb_hemisphere_match,
                             fafb_hemisphere_match_quality, d
                             ataset, match_dataset, user from %s", "matches"), base = tab)
  ft = subset(ft, ft$dataset == dataset)
  if(dataset=="flywire"){
    ft$root_id = ft$id
  }else if(dataset=="hemibrain"){
    ft$bodyid = ft$id
  }
}

flytable_matches <- function(dataset=c("flywire","hemibrain")){
  ft=get_matching_flytable(dataset=dataset)
  ft = as.data.frame(ft)
  rownames(ft) = ft$id
  ft
}

flytable_matches_update <- function(df, dataset=c("flywire","hemibrain"), DryRun = FALSE){
  ft=flytable_matches(dataset=dataset)
  tab = switch(dataset, flywire = {"matching"}, hemibrain = {"hemibrain_matches"})
  toupdate = setdiff(df, ft)
  if(DryRun)
    toupdate
  else {
    fafbseg::flytable_update_rows(toupdate, table = tab)
    invisible(toupdate)
  }
}

flytable_matches_add <- function(df = NULL, DryRun = FALSE){
  if(is.null(df)){
    df = flytable_meta()
  }
  ft=fafbseg::flytable_query("select _id, id, side, cell, cell_type, cell_body_fiber, ito_lee_hemilineage, match, quality, fafb_hemisphere_match, fafb_hemisphere_match_quality, dataset, match_dataset, user from matching")
  toupdate = setdiff(df, ft)
  if(DryRun)
    toupdate
  else {
    fafbseg::flytable_update_rows(toupdate, table = "matching")
    invisible(toupdate)
  }
}

flytable_update_proofread <- function(DryRun=TRUE, update_roots=TRUE) {
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
  else {
    fafbseg::flytable_update_rows(toupdate, table = "info")
    invisible(toupdate)
  }
}

flytable_update_status <- function(DryRun=TRUE) {
  ft=ft.orig=fafbseg::flytable_query("select _id, root_id, supervoxel_id, status from info")
  ft$status = standard_statuses(status)
  toupdate=ft %>%
    dplyr::mutate(ft.new=root_id %in% as.character(psp$pt_root_id) & !root_id==0) %>%
    dplyr::filter(proofread.new!=proofread) %>%
    dplyr::select(`_id`, proofread.new) %>%
    dplyr::rename(proofread=proofread.new)
  if(DryRun)
    toupdate
  else {
    fafbseg::flytable_update_rows(toupdate, table = "info")
    invisible(toupdate)
  }
}

flytable_update_transmitter <- function(DryRun=T, update_roots=T) {
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
  else {
    fafbseg::flytable_update_rows(toupdate, table = "info")
    invisible(toupdate)
  }
}


