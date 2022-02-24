#' Create a flywire scene that shows your chosen flywire and hemibrain neuron meshes
#'
#'@description Return a sample Neuroglancer scene URL for FlyWire + hemibrain  dataset
#'
#' @param bodyid character,  a vector of one or more valid hemibrain body IDs
#' @param flywire.id character, a vector of one or more valid flywire root IDs
#' @param ... Additional arguments passed to \code{fafbseg::ngl_encode_url}.
#'
#' @examples
#' \donttest{
#' \dontrun{
#' flywire.id = c("720575940611249010","720575940626212584")
#' bodyid = c("1540129092","1513363614")
#' scene = flywire_hemibrain_scene(bodyid, flywire.id)
#' scene
#'}}
#'@return A URL
#'
#'@export
#'@seealso \code{\link{hemibrain_read_neurons}}
flywire_hemibrain_scene <- function(bodyid,
                                    flywire.id,
                                    ...
                                         ){
  if(missing(flywire.id)) {
    stopifnot(is.data.frame(bodyid) && ncol(bodyid)==2)
    flywire.id=bodyid[[2]]
    bodyid=bodyid[[1]]
  }
  # Colours, hemibrain blues
  blues = grDevices::colorRampPalette(c(hemibrain_bright_colors[["cyan"]], hemibrain_bright_colors[["blue"]],hemibrain_bright_colors[["navy"]]))
  hb.blues = blues(length(bodyid))
  names(hb.blues) = bodyid
  reds = grDevices::colorRampPalette(c(hemibrain_bright_colors[["pink"]],hemibrain_bright_colors[["cerise"]],hemibrain_bright_colors[["red"]], hemibrain_bright_colors[["darkred"]]))
  fw.reds = reds(length(flywire.id))
  names(fw.reds) = flywire.id
  # Colours, flywire reds
  sriurl2 <- "https://ngl.flywire.ai/?json_url=https://globalv1.flywire-daf.com/nglstate/5891269037391872"
  sriscene <- fafbseg::ngl_decode_scene(fafbseg::flywire_expandurl(sriurl2))
  sriscene$layers$`Production-segmentation_with_graph`$segments=as.list(as.character(flywire.id))
  sriscene$layers$hemibrain_meshes$segments=as.list(as.character(bodyid))
  sriscene <- ngl_add_colours2(sriscene, hb.blues,
                               layer = "hemibrain_meshes")
  sriscene <- ngl_add_colours2(sriscene, fw.reds,
                              layer = "Production-segmentation_with_graph")
  fafbseg::ngl_encode_url(sriscene, ...)
}

#'@export
flywire_gradient_scene <- function(flywire.id, ...){
  f <- fafbseg::flywire_scene(flywire.id, open = FALSE)
  d <- fafbseg::ngl_decode_scene(f)
  cols <- rainbow(length(flywire.id))
  names(cols) <- flywire.id
  sriscene <- ngl_add_colours2(d, cols,
                               layer = "Production-segmentation_with_graph")
  fafbseg::ngl_encode_url(sriscene, ...)
}

# hidden
make_partner_scene <- function(target.bodyid, partner.bodyid, partner.flywire.id, relationship=c("us", "ds"), ...) {
  relationship=match.arg(relationship)
  sriurl2 <- "https://ngl.flywire.ai/?json_url=https://globalv1.flywire-daf.com/nglstate/5064012857344000"
  sriscene <- fafbseg::ngl_decode_scene(flywire_expandurl(sriurl2))
  sriscene$layers$flywire_partner$segments=as.list(as.character(partner.flywire.id))
  sriscene$layers$hemibrain_target$segments=as.list(as.character(target.bodyid))
  sriscene$layers$hemibrain_partner$segments=as.list(as.character(partner.bodyid))
  fwcands=flywire_partner_summary(partner.flywire.id, Verbose = FALSE,
                                  partners = ifelse(relationship=="us", "inputs", "outputs"))
  if(nrow(fwcands)>20) fwcands=fwcands[1:20,]
  idfield=ifelse(relationship=="us", "pre_id", "post_id")

  sriscene$layers$candidates$segments=as.list(as.character(fwcands[[idfield]]))
  fafbseg::ngl_encode_url(sriscene, ...)
}

tmake_partner_scene <- function(...) {
  tryCatch(make_partner_scene(...), error=function(e) {warning(e); NA_character_})
}

#' Set colour for a specified layer properly
ngl_add_colours2 <- function(sc, colours, layer){
  stopifnot(fafbseg:::is.ngscene(sc))
  layers=ngl_layers(sc)
  stopifnot(layer %in% names(layers))
  sc.copy=sc
  fafbseg::ngl_layers(sc.copy)=layers[layer]
  # this is to work around an issue where there must be at least one segmen
  # for ngl_segments (IIRC) to wrk
  if(length(ngl_segments(sc.copy))==0)
    ngl_segments(sc.copy)=0
  sc.copy <- fafbseg::ngl_add_colours(sc.copy, colours)
  layers[[layer]]=sc.copy$layers[[1]]
  ngl_layers(sc)=layers
  sc
}

# hidden
make_top_partners_scene <- function(bodyid, flywire.id,
                                    partners=c("outputs", "inputs"),
                                    cleft.threshold = 30, n=20,
                                    col=rainbow, ...) {
  if(missing(flywire.id)) {
    stopifnot(is.data.frame(bodyid) && ncol(bodyid)==2)
    flywire.id=bodyid[[2]]
    bodyid=bodyid[[1]]
  }
  partners=match.arg(partners)
  cols=if(is.function(col)) col(n) else col

  sriurl2 <- "https://ngl.flywire.ai/?json_url=https://globalv1.flywire-daf.com/nglstate/5217719817338880"
  sriscene <- fafbseg::ngl_decode_scene(flywire_expandurl(sriurl2))

  sriscene$layers$fw_query$segments=as.list(as.character(flywire.id))
  sriscene$layers$hb_query$segments=as.list(as.character(bodyid))
  hbp=neuprintr::neuprint_connection_table(bodyid,
                                           prepost = ifelse(partners=='inputs', 'PRE', "POST"))
  if(nrow(hbp)>n) hbp=hbp[seq_len(n),]
  hbp$col=cols[seq_len(nrow(hbp))]
  sriscene <- ngl_add_colours2(sriscene, hbp[c("partner", "col")],
                               layer = "hb_partners")
  fwcands=fafbseg::flywire_partner_summary(flywire.id, Verbose = FALSE,
                                  partners = partners,
                                  cleft.threshold = cleft.threshold)
  if(nrow(fwcands)>n) fwcands=fwcands[seq_len(n),]
  idfield=ifelse(partners=="inputs", "pre_id", "post_id")
  fwcands$col=cols[seq_len(nrow(fwcands))]
  sriscene <- ngl_add_colours2(sriscene, fwcands[c(idfield, "col")],
                               layer = "fw_partners")
  fafbseg::ngl_encode_url(sriscene, ...)
}

# hidden
tmake_toppartners_scene=function(...) {
  tryCatch(make_top_partners_scene(...), error=function(e) {warning(e); NA_character_})
}
