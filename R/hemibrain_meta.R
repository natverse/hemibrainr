#######################################################################################
################################ hemibrain meta data ##################################
#######################################################################################

#' Get meta information for hemibrain neurons
#'
#' @description Read meta data for hemibrain neurons from NeuPrint and supplement this with other data available in this package.
#' Specifically, neuron lineage, class, olfactory system layer and FAFB match information,
#' as well as numbers related to neurons' putative axon-dendrite split, e.g. for synapses or cable_length.
#'
#' @param x a vector of bodyids that can be read from \url{'https://neuprint.janelia.org/'}.
#' @param ... arguments passed to \code{neuprintr::neuprint_get_meta}.
#'
#' @return a \code{data.frame} with columns that can give a user neuron lineage and class information, as well as numbers
#' related to neurons' putative axon-dendrite split.
#'
#' #' \itemize{
#'
#'   \item{"bodyid"} { - unique identifier for a single hemibrain neuron.}
#'
#'   \item{"pre"} { - the number of presynapses (outputs) a neuron contains. Each of this is polyadic.}
#'
#'   \item{"post"} { - the number of postsynapses (inputs) to the neuron.}
#'
#'   \item{"upstream"} { - the number of incoming connections to a neuron.}
#'
#'   \item{"downstream"} { - the number of outgoing connections from a neuron.}
#'
#'   \item{"soma"} { - neuron size in voxels.}
#'
#'   \item{"sid"} { - which brain hemisphere contains the neuron's soma.}
#'
#'   \item{"voxels"} { - whether the neuron has a soma in the hemibrain volume.}
#'
#'   \item{"name"} { - the name of this neuron.}
#'
#'   \item{"connectivity.type"}{ - the neuPrint designated 'type' for the neuron. A subset of neurons within a cell_type that share similar connectivity. A connectivity
#'   type is distinguished from a cell_type by an ending _letter, defined using CBLAST (Scheffer et al. 2020).}
#'
#'   \item{"cell_type"}{ - neurons of a shared morphology that take the same cell body fibre tract and come from the same hemilineage (Bates 2019).}
#'
#'   \item{"cellBodyFiber"}{ - The cell body fiber for a neuron, as read from neuPrint and annotated by a team under \href{https://www.janelia.org/people/kei-ito}{Kei Ito}.}
#'
#'   \item{"putative.classic.transmitter"}{ - putative neurotransmitter based on what neurons in the hemilineage in question has been shown to express, out of acetylcholine, GABA, glutamate.}
#'
#'   \item{"ito_lee_hemilineage"}{ - the hemilineage that we reckon this cell_type belongs to, based on expert review of light level data from the K. Ito and T. Lee groups (Yu 2013, Ito 2013).
#'    See \code{\link{hemibrain_hemilineages}}.}
#'
#'   \item{"hartenstein_hemilineage"}{ - the hemilineage that we reckon this cell_type belongs to, based on expert review of light level data from the V. Hartenstein Lee group (Wong 2013, Lovick 2013).
#'    See \code{\link{hemibrain_hemilineages}}.}
#'
#'   \item{"FAFB.match"}{ - the ID of the manual match from the FAFB data set. ID indicates a neuron reconstructed in FAFBv14 CATMAID. Many of these neurons will be available through Virtual Fly Brain.}
#'
#'   \item{"FAFB.match.quality"}{ - the matcher makers qualitative assessment of how good this match is: a poor match could be a neuron from a very similar same cell_type or a highly untraced neuron that may be the correct cell_type. An okay match should be a neuron that looks to be from the same morphological cell_type but there may be some discrepancies in its arbour. A good match is a
#'   neuron that corresponds well between FAFB and the hemibrain data.}
#'
#'   \item{"layer"}{ - probabilistic mean path length from neuron from ALRNs, depends on connection strengths.}
#'
#'   \item{"layer.ct"}{- the mean layer for cell_type, rounded to the nearest whole number.}
#'
#'   \item{"axon_outputs"}{- number of outgoing connections from the neuron's predicted axon.}
#'
#'   \item{"dend_outputs"}{- number of outgoing connections from the neuron's predicted dendrite.}
#'
#'   \item{"axon_inputs"}{- number of incoming connections from the neuron's predicted axon.}
#'
#'   \item{"dend_inputs"}{- number of incoming connections from the neuron's predicted dendrite.}
#'
#'   \item{"total_length"}{- total cable_length of the neuron in micrometers.}
#'
#'   \item{"axon_length"}{- total axon cable_length of the neuron in micrometers.}
#'
#'   \item{"dend_length"}{- total dendrite cable_length of the neuron in micrometers.}
#'
#'   \item{"pd_length"}{- total  cable_length of the primary dendrite 'linker' between axon and dendrite.}
#'
#'   \item{"segregation_index"}{- a quantification of how polarised a neuron is, in terms of its segregation of inputs onto its predicted dendrite and outputs onto its axon,
#'    where 0 is no-polarisation and 1 is totally polarised (Schneider-Mizell 2016).}
#'
#'}
#'
#' @examples
#' \donttest{
#' \dontrun{
#'
#' # Get neurons
#' meta = hemibrain_get_meta("763730414")
#' View(meta)
#'
#' }}
#' @rdname hemibrain_get_meta
#' @export
#' @seealso \code{\link{hemibrain_splitpoints}}, \code{\link{hemibrain_flow_centrality}}, \code{\link{hemibrain_somas}},
#' \code{\link{hemibrain_precomputed_splitpoints}}, \code{\link{hemibrain_metrics}},\code{\link{hemibrain_olfactory_layers}}
#' \code{\link{hemibrain_hemilineages}}, \code{\link{classed.ids}}, \code{\link{hemibrain_read_neurons}}
hemibrain_get_meta <- function(x, ...){
  # Get information from neuprint
  nmeta = neuprintr::neuprint_get_meta(x, possibleFields = neuprintr::neuprint_get_fields(), ...)
  if(!nrow(nmeta)){
    warning("No meta data could be retieved, returning NULL")
    NULL
  }else{
    # Add cell. type
    nmeta$connectivity.type = nmeta$type
    nmeta$cell_type = gsub("_[a-z]{1}$","",nmeta$type)
    nmeta$cell_type[is.na(nmeta$cell_type)] = "uncertain"
    nmeta$connectivity.type[is.na(nmeta$connectivity.type)] = "uncertain"

    # Add lineage information
    hemibrain_hemilineages <- dplyr::rename(.data = hemibrainr::hemibrain_hemilineages, lineage_type = type)
    hemibrain_hemilineages <- hemibrain_hemilineages[!duplicated(hemibrain_hemilineages$cell_body_fiber),]
    nmeta2 <- dplyr::left_join(nmeta,
                              hemibrain_hemilineages,
                              by=c("cellBodyFiber"="cell_body_fiber"))
    nmeta2$FAFB <- NULL

    # Neuron class
    nmeta$class = NA
    nmeta$class[nmeta$bodyid%in%hemibrainr::dn.ids] = "DN"
    nmeta$class[nmeta$bodyid%in%hemibrainr::ton.ids] = "TOON"
    nmeta$class[nmeta$bodyid%in%hemibrainr::lhn.ids] = "LHN"
    nmeta$class[nmeta$bodyid%in%hemibrainr::rn.ids] = "ALRN"
    nmeta$class[nmeta$bodyid%in%hemibrainr::orn.ids] = "ORN"
    nmeta$class[nmeta$bodyid%in%hemibrainr::hrn.ids] = "HRN"
    nmeta$class[nmeta$bodyid%in%hemibrainr::pn.ids] = "ALPN"
    nmeta$class[nmeta$bodyid%in%hemibrainr::upn.ids] = "uPN"
    nmeta$class[nmeta$bodyid%in%hemibrainr::mpn.ids] = "mPN"
    nmeta$class[nmeta$bodyid%in%hemibrainr::vppn.ids] = "VPPN"
    nmeta$class[nmeta$bodyid%in%hemibrainr::alln.ids] = "ALLN"
    nmeta$class[nmeta$bodyid%in%hemibrainr::dan.ids] = "DAN"
    nmeta$class[nmeta$bodyid%in%hemibrainr::mbon.ids] = "MBON"
    nmeta$class[nmeta$bodyid%in%hemibrainr::kc.ids] = "KC"
    nmeta$class[nmeta$bodyid%in%hemibrainr::apl.ids] = "APL"
    nmeta$class[nmeta$bodyid%in%hemibrainr::cent.ids] = "LHCENT"
    nmeta$class[nmeta$bodyid%in%hemibrainr::lc.ids] = "LCPN"

    # Correct soma location
    if(!is.null(nmeta2$somaLocation)){
      nmeta2$somaLocation <- sapply(nmeta2$somaLocation, function(entry){
        if(grepl('list',entry)){
          gsub("^list\\(|\\).*","",entry)
        }else{
          entry
        }
      })
      nmeta2$somaLocation[nmeta2$somaLocation=="NA"]=NA
    }

    # # Add match information
    # nmeta2$FAFB.match = hemibrainr::hemibrain_matched[as.character(nmeta2$bodyid),"match"]
    # nmeta2$FAFB.match.quality = hemibrainr::hemibrain_matched[as.character(nmeta2$bodyid),"quality"]

    # Add olfactory layer information
    nmeta2$layer = hemibrainr::hemibrain_olfactory_layers[match(nmeta2$bodyid,hemibrainr::hemibrain_olfactory_layers$node),"layer_mean"]
    nmeta2$ct.layer = NA
    for(ct in unique(nmeta2$type)){
      layer = round(mean(subset(nmeta2,nmeta2$type==ct)$layer))
      nmeta2$ct.layer[nmeta2$type==ct] = layer
    }

    # Add split information
    selcols = c("soma_edit", "skeletonization", "edited_cable",
                "axon_outputs", "dend_outputs", "axon_inputs",
                "dend_inputs", "total_outputs_density", "total_inputs_density",
                "axon_outputs_density", "dend_outputs_density", "axon_inputs_density",
                "dend_inputs_density", "total_length", "axon_length", "dend_length",
                "pd_length", "segregation_index")
    selcols = intersect(selcols, colnames(hemibrainr::hemibrain_metrics))
    selcols = setdiff(selcols, colnames(nmeta2))
    hemibrain_metrics_sel = hemibrainr::hemibrain_metrics[as.character(nmeta2$bodyid), selcols]
    nmeta2 = cbind(nmeta2, hemibrain_metrics_sel)
    rownames(nmeta2) = nmeta2$bodyid

    # Return
    nmeta2
  }
}

#' Set annotations for FAFB neurons in CATMAID based on hemibrain results
#'
#' @description Set annotations for FAFB neurons in CATMAID based on matches made to hemibrain neurons.
#' E.g. transfer information on matches and cell body fibers, and also update lineage related information.
#' Not that \code{catmaid::flywire_matching_rewrite()} will write annotations related to flywire.
#'
#' @param x a vector of FAFB skeleton IDs (skids) that can be read with \code{catmaid::catmaid_skids()}
#' @param find an annotation/search term/vector of skids passed to \code{catmaid::catmaid_skids}.
#' @param dataset either FAFB or hemibrain. If hemibrain, then \code{find} must be a vector of valid bodyids. Matches to FAFB neurons will be annotated for these bodyids, using \code{hemibrain_matches}.
#' @param ito_lee_hemilineage character, the correct K. Ito / T. Lee hemilineage. Must be an entry in \code{\link{hemibrain_hemilineages}}. If set to 'wipe' all lineage related annotations will be removed.
#' @param transmitter character, the known or putative transmitter identity of the neurons that will be fetched using \code{find}.
#' @param delete.find logical, is \code{TRUE} then if \code{find} is an annotation, it will be wiped from the neuron after this function is used (if you have permission to remove it)
#' @param putative if \code{TRUE} the word 'putative' is added to all lineage annotations to indicate to users that the labelled neurons may not have been given the correct lineages. More work may need to be done to solidify these assignments.
#' @param flywire logical. Whether or not to add flywire-related annotation information, e.g. a corresponding position in flywire voxel space and a matching, cognate flywire ID.
#' @param ... arguments passed to \code{neuprintr::neuprint_get_meta} and \code{catmaid::catmaid_login}.
#'
#' @return annotations set on CATMAID neurons in specified CATMAID instance, see your \code{catmaid::catmaid_connection()}
#'
#' @examples
#' \donttest{
#' \dontrun{
#'
#' # Transfer information on the VA6 uPN
#' fafb_hemibrain_annotate("16")
#'
#' # Annotate all GABAergic MBONs in FAFB as GABAergic
#' bodyids = c(1016835041,
#'          1048215779,
#'         1078693835,
#'          672352543,
#'          613719036,
#'          581678043,
#'          673085197,
#'          706840303,
#'          768555687,
#'          5812982924,
#'          424767514,
#'          517518166);
#'fafb_set_transmitter(find=bodyids,dataset="hemibrain",transmitter="GABA", putative = FALSE)
#'
#' }}
#' @rdname fafb_hemibrain_annotate
#' @export
#' @seealso \code{\link{hemibrain_get_meta}}
fafb_hemibrain_annotate <- function(x,
                                    flywire = TRUE,
                                    ...){
  if(flywire){
    if(!requireNamespace("fafbseg", quietly = TRUE)) {
      stop("Please install fafbseg using:\n", call. = FALSE,
           "remotes::install_github('natverse/fafbseg')")
    }
  }

  # Get matches
  matches = hemibrain_matches()
  matches = subset(matches,
                   matches$quality %in% c("good","medium")
                   & matches$dataset =="FAFB"
                   & ! matches$match %in% c("none","","NA"," ")) # Exclude poor for now
  hemi.meta = hemibrain_get_meta(unique(matches$match), OmitFailures = TRUE, ...)
  matches = subset(matches, matches$match %in% hemi.meta$bodyid)
  x = intersect(x,rownames(matches))

  for(i in x){

    # Print neuron under consideration
    print(knitr::kable(matches[i,c("id", "cell_type", "ito_lee_hemilineage", "match", "quality")]))

    # Get old annotations
    a = catmaid::catmaid_get_annotations_for_skeletons(i, ...)
    mat = a$annotation[grepl("hemibrain_match: |Hemibrain_match: ",a$annotation)]
    ct = a$annotation[grepl("Cell_type: |cell_type: ",a$annotation)]
    cbf = a$annotation[grepl("cellBodyFiber: |CellBodyFiber: ",a$annotation)]

    # Get flywire information
    if(flywire){
      fid = a$annotation[grepl("root_id: |root_id: ",a$annotation)]
      cat = catmaid::read.neurons.catmaid(i, ...)
      # Get xyz for root points
      roots = sapply(cat, function(y) nat::xyzmatrix(y)[nat::rootpoints(y),])
      roots = t(roots)
      fafb_xyz = apply(roots, 1, paste, collapse = ",")

      # Get FlyWire voxel coordinates
      roots.flywire = nat.templatebrains::xform_brain(roots, reference = "FlyWire", sample = "FAFB14", .parallel = TRUE, verbose = TRUE)
      rownames(roots.flywire) = rownames(roots)
      roots.flywire.raw = scale(roots.flywire, scale = c(4, 4, 40), center = FALSE)
      fw.ids = fafbseg::flywire_xyz2id(roots.flywire.raw, rawcoords = TRUE)
      fw.ids[fw.ids=="0"] = NA
      flywire_xyz = apply(roots.flywire.raw, 1, paste, collapse = ",")

      # Make annotation
      afid = paste0('root_id: ', fw.ids[1])
      message(afid)
      afw = paste0('flywire_xyz: ', flywire_xyz[1])
      message(afw)
      afafb = paste0('fafb_xyz: ', fafb_xyz[1])
      message(afafb)

    }

    # Get meta for matches
    y = matches[i,"match"]
    hmeta = subset(hemi.meta, hemi.meta$bodyid == y)

    # Assemble annotations
    am = paste0('hemibrain_match: ', y)
    message(am)
    act = paste0('cell_type: ', hmeta[,"type"])
    message(act)
    acbf = paste0('cellBodyFiber: ', hmeta[,"cellBodyFiber"])
    message(acbf)

    # next
    if(am%in%mat & act%in%ct & acbf%in%cbf){
      message("already annotated, skipping ...")
      next
    }

    # Set annotations
    annotations = c(am,
                    act,
                    acbf)
    if(flywire){
      annotations = c(am,
                      act,
                      acbf,
                      afid,
                      afw,
                      afafb)
    }
    catmaid::catmaid_set_annotations_for_skeletons(skids = i,
                                                   annotations = c(am,
                                                                   act,
                                                                   acbf),
                                                   force = TRUE, ...)

    # Meta-annotate
    catmaid::catmaid_set_meta_annotations(meta_annotations = "cell_type", annotations = act, ...)
    catmaid::catmaid_set_meta_annotations(meta_annotations = "cell body fiber", annotations = acbf, ...)
    catmaid::catmaid_set_meta_annotations(meta_annotations = "hemibrain_match", annotations = am, ...)
    if(flywire){
      catmaid::catmaid_set_meta_annotations(meta_annotations = "root_id", annotations = afid, ...)
      catmaid::catmaid_set_meta_annotations(meta_annotations = "flywire_xyz", annotations = afw, ...)
      catmaid::catmaid_set_meta_annotations(meta_annotations = "fafb_xyz", annotations = afafb, ...)
    }

    # remove annotions
    if(length(ct)){
      tryCatch(catmaid::catmaid_remove_annotations_for_skeletons(skids = i, annotations = paste0("annotation:^",ct,"$"), force = TRUE, ...),
               error = function(e) NULL)
    }
    if(length(cbf)){
      tryCatch(catmaid::catmaid_remove_annotations_for_skeletons(skids = i, annotations = paste0("annotation:^",cbf,"$"), force = TRUE, ...),
               error = function(e) NULL)
    }
  }

}

#' @rdname fafb_hemibrain_annotate
#' @export
fafb_set_hemilineage <- function(find,
                                 ito_lee_hemilineage,
                                 dataset = c("FAFB","hemibrain"),
                                 putative = TRUE,
                                 delete.find = FALSE,
                                 ...){
  dataset = match.arg(dataset)
  if(dataset=="hemibrain"){
    matches = hemibrain_matches()
    matches = subset(matches, matches$quality %in% c("good","medium") & matches$dataset == dataset)
    find = as.character(find)
    find = unique(matches[find,"match"])
    find = find[!is.na(find)]
    skds = find[!find%in%c("none",""," ")]
  }else{
    # Find neurons based on an annotation or as skeleton IDs
    skds = catmaid::catmaid_skids(x = find, ...)
    if(grepl("ItoLee|Hartenstein",find)){
      a = catmaid::catmaid_get_annotations_for_skeletons(find)
      b = subset(a, a$uid==109 & a$annotation == find) # Just trust annotations from asbates
      skds = unique(b$skid)
    }
  }
  if(!length(skds)){
    stop("FAFB skeleton IDs could not be found")
  }
  if(ito_lee_hemilineage!="wipe"){
    # Make sure we have a valid hemilineage
    if(length(ito_lee_hemilineage)>1|length(ito_lee_hemilineage)==0){
      stop("Please provide a single ito_lee_hemilineage")
    }
    if(!ito_lee_hemilineage%in%hemibrainr::hemibrain_hemilineages$ito_lee_hemilineage){
      stop("Please provide a valid ito_lee_hemilineage. See 'hemibrain_hemilineages'")
    }
    # Get meta information
    i = ito_lee_hemilineage
    ito_lee_lineage = subset(hemibrainr::hemibrain_hemilineages, hemibrainr::hemibrain_hemilineages$ito_lee_hemilineage == i)$ito_lee_lineage[1]
    hartenstein_lineage = subset(hemibrainr::hemibrain_hemilineages, hemibrainr::hemibrain_hemilineages$ito_lee_hemilineage == i)$hartenstein_lineage[1]
    hartenstein_hemilineage = subset(hemibrainr::hemibrain_hemilineages, hemibrainr::hemibrain_hemilineages$ito_lee_hemilineage == i)$hartenstein_hemilineage[1]
    # Make new annotations
    ito_lee_lineage = paste0("ito_lee_lineage: " , ito_lee_lineage)
    ito_lee_hemilineage = paste0("ito_lee_hemilineage: " , ito_lee_hemilineage)
    hartenstein_lineage = paste0("hartenstein_lineage: " , hartenstein_lineage)
    hartenstein_hemilineage = paste0("hartenstein_hemilineage: " , hartenstein_hemilineage)
    # Set these as putative?
    if(putative){
      ito_lee_lineage = paste0(ito_lee_lineage, "_putative")
      ito_lee_hemilineage = paste0(ito_lee_hemilineage, "_putative")
      hartenstein_lineage = paste0(hartenstein_lineage, "_putative")
      hartenstein_hemilineage = paste0(hartenstein_hemilineage, "_putative")
      lins = c(ito_lee_lineage,
               ito_lee_hemilineage,
               hartenstein_lineage,
               hartenstein_hemilineage,
               "Lineage_putative")
      message("Adding annotations:   ", paste(lins, collapse = "  ", sep = "  "))
      # Set annotations
      catmaid::catmaid_set_annotations_for_skeletons(skids = skds,
                                                     annotations = lins, force = TRUE, ...)
    }else{
      # Set annotations
      lins = c(ito_lee_lineage,
               ito_lee_hemilineage,
               hartenstein_lineage,
               hartenstein_hemilineage,
               "Lineage_annotated")
      message("Adding annotations:   ", paste(lins, collapse = "  ", sep = "  "))
      catmaid::catmaid_set_annotations_for_skeletons(skids = skds,
                                                     annotations = lins,
                                                     force = TRUE, ...)
      catmaid::catmaid_set_meta_annotations(meta_annotations = "ito_lee_lineage", annotations = ito_lee_lineage,...)
      catmaid::catmaid_set_meta_annotations(meta_annotations = "ito_lee_hemilineage", annotations = ito_lee_hemilineage,...)
      catmaid::catmaid_set_meta_annotations(meta_annotations = "hartenstein_lineage", annotations = hartenstein_lineage,...)
      catmaid::catmaid_set_meta_annotations(meta_annotations = "hartenstein_hemilineage", annotations = hartenstein_hemilineage,...)
    }
  }

  # Delete old annotations
  a = catmaid::catmaid_get_annotations_for_skeletons(skds)
  old = unique(a$annotation[grepl("ItoLee_|Hartenstein_|Volker_|volker_|Itolee_|itolee_|itoLee_",a$annotation)])
  old = setdiff(old,lins)
  # Remove
  if(length(old)){
    message("Removing annotations:   ", paste(old, collapse = "  ", sep = "  "))
    old = paste0("^",old,"$")
    catmaid::catmaid_remove_annotations_for_skeletons(skids = skds,
                                                      annotations = old,
                                                      force = TRUE, ...)
  }
  if(delete.find){
    catmaid::catmaid_remove_annotations_for_skeletons(skids = skds,
                                                      annotations = paste0("^",find,"$"),
                                                      force = TRUE, ...)
  }

}

#' @rdname fafb_hemibrain_annotate
#' @export
fafb_set_transmitter <-function(find,
                                dataset = c("FAFB","hemibrain"),
                                putative = TRUE,
                                transmitter = c("acetylcholine","GABA","glutamate"),
                                delete.find = FALSE,
                                ...){
  dataset = match.arg(dataset)
  transmitter = match.arg(transmitter)
  if(dataset=="hemibrain"){
    matches = hemibrain_matches()
    matches = subset(matches, matches$quality %in% c("good","medium") & matches$dataset == dataset)
    find = as.character(find)
    find = unique(matches[find,"match"])
    find = find[!is.na(find)]
    skds = find[!find%in%c("none",""," ","NA","unknown")]
  }else{
    skds = catmaid::catmaid_skids(find, ...)
  }
  if(!length(skds)){
    stop("FAFB skeleton IDs could not be found")
  }
  possible.transmitters = c("Vglut","ACh","Ach","acetylcholine","ChAT", "ChA", "GABA","glutamate","Glutamate",
                            "serotonine","dopamine","octopammine","SIFamide","unknown","NA")
  if(transmitter!="wipe"){
    # Set these as putative?
    if(putative){
      ann = paste0("putative transmitter: " , transmitter)
      anns = c(ann,
               "transmitter_putative")
      message("Adding annotations:   ", paste(anns, collapse = "  ", sep = "  "))
      # Set annotations
      catmaid::catmaid_set_annotations_for_skeletons(skids = skds,
                                                     annotations = anns, force = TRUE, ...)
      catmaid::catmaid_set_meta_annotations(meta_annotations = "transmitter", annotations = ann,...)
    }else{
      # Set annotations
      ann = paste0("known transmitter: " , transmitter)
      anns = c(ann,
               "transmitter_annotated")
      message("Adding annotations:   ", paste(anns, collapse = "  ", sep = "  "))
      catmaid::catmaid_set_annotations_for_skeletons(skids = skds,
                                                     annotations = anns,
                                                     force = TRUE, ...)
      catmaid::catmaid_set_meta_annotations(meta_annotations = "transmitter", annotations = ann,...)
    }
  }
  # Delete old annotations
  a = catmaid::catmaid_get_annotations_for_skeletons(skds)
  b = unique(a[grepl("Putative_Neurotransmitter|putative neurotransmitter|known transmitter|putative transmitter|known neurotransmitter",a$annotation),])
  old = unique(b$annotation[grepl(paste(possible.transmitters,collapse="|"),b$annotation)])
  old = setdiff(old,anns)
  if(length(old)){
    message("Removing annotations:   ", paste(old, collapse = "  ", sep = "  "))
    old = paste0("^",old,"$")
    # Remove
    if(length(old)){
      catmaid::catmaid_remove_annotations_for_skeletons(skids = skds,
                                                        annotations = old,
                                                        force = TRUE, ...)
    }
    if(delete.find){
      catmaid::catmaid_remove_annotations_for_skeletons(skids = skds,
                                                        annotations = paste0("^",find,"$"),
                                                        force = TRUE, ...)
    }
  }
}

# Hidden
## Don't use unless very necessary
fafb_lineage_complete_wipe <- function(server = "v14seg-Li-190411.0", ...){
  # Choose CATMAID server
  conn = catmaid::catmaid_login(...)
  if(!is.null(server)){
    conn$server = paste0("https://neuropil.janelia.org/tracing/fafb/",server,"/")
  }
  a = catmaid::catmaid_get_annotationlist(..., conn = conn)
  instance = catmaid::catmaid_login(conn=conn, ...)$server
  message("Looking at ", instance)
  old = a$annotations$name[grepl("ItoLee_|Hartenstein_|Volker_|volker_|Itolee_|itolee_|itoLee_|Lineage_",unlist(a$annotations$name))]
  for(ann in old){
    message("Removing all instances of annotation: ", ann)
    ann2 = paste0("annotation:^",ann,"$")
    skds = catmaid::catmaid_skids(ann2, conn = conn, ...)
    catmaid::catmaid_remove_annotations_for_skeletons(skids = skds,
                                                      annotations = ann2,
                                                      force = TRUE,
                                                      conn = conn,
                                                      ...)
  }
}
