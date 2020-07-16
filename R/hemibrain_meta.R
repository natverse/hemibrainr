#######################################################################################
################################ hemibrain meat data ##################################
#######################################################################################

#' Get meta information for hemibrain neurons
#'
#' @description Read meta data for hemibrain neurons from NeuPrint and supplement this with other data available in this package.
#' Specifically, neuron lineage, class, olfactory system layer and FAFB match information,
#' as well as numbers related to neurons' putative axon-dendrite split, e.g. for synapses or cable length.
#'
#' @param x a vector of bodyids that can be read from \url{'https://neuprint.janelia.org/'}.
#' @param ... arguments passed to \code{neuprintr::neuprint_get_meta}.
#'
#' @return a \code{data.frame} with columns that can give a user neuron lineage and class information, as well as numbers
#' related to neurons' putative axon-dendrite split.
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
  nmeta = neuprintr::neuprint_get_meta(x, ...)

  # Add lineage information
  nmeta2 = merge(nmeta,hemibrain_hemilineages, all.x = TRUE, all.y = FALSE)
  nmeta2$FAFB = NULL

  # Neuron class
  nmeta$class = NA
  nmeta$class[nmeta$bodyid%in%dn.ids] = "DN"
  nmeta$class[nmeta$bodyid%in%ton.ids] = "TON"
  nmeta$class[nmeta$bodyid%in%lhn.ids] = "LHN"
  nmeta$class[nmeta$bodyid%in%rn.ids] = "RN"
  nmeta$class[nmeta$bodyid%in%orn.ids] = "ORN"
  nmeta$class[nmeta$bodyid%in%hrn.ids] = "HRN"
  nmeta$class[nmeta$bodyid%in%pn.ids] = "PN"
  nmeta$class[nmeta$bodyid%in%upn.ids] = "uPN"
  nmeta$class[nmeta$bodyid%in%mpn.ids] = "mPN"
  nmeta$class[nmeta$bodyid%in%vppn.ids] = "VPPN"
  nmeta$class[nmeta$bodyid%in%alln.ids] = "ALLN"
  nmeta$class[nmeta$bodyid%in%dan.ids] = "DAN"
  nmeta$class[nmeta$bodyid%in%mbon.ids] = "MBON"

  # Add match information
  nmeta2$FAFB.match = hemibrain_matched[as.character(nmeta2$bodyid),"match"]
  nmeta2$FAFB.match.quality = hemibrain_matched[as.character(nmeta2$bodyid),"quality"]

  # Add olfactory layer information
  nmeta2$layer = hemibrain_olfactory_layers[match(nmeta2$bodyid,hemibrain_olfactory_layers$node),"layer_mean"]
  nmeta2$ct.layer = NA
  for(ct in unique(nmeta2$type)){
    layer = round(mean(subset(nmeta2,type==ct)$layer))
    nmeta2$ct.layer[nmeta2$type==ct] = layer
  }

  # Add split information
  selcols=setdiff(colnames(hemibrainr::hemibrain_metrics), colnames(nmeta2))
  hemibrain_metrics_sel = hemibrainr::hemibrain_metrics[as.character(nmeta2$bodyid), selcols]
  nmeta2 = cbind(nmeta2, hemibrain_metrics_sel)
  rownames(nmeta2) = nmeta2$bodyid

  # Return
  nmeta2
}


#' Set annotations for FAFB neurons in CATMAID based on hemibrain results
#'
#' @description Set annotations for FAFB neurons in CATMAID based on matches made to hemibrain neurons.
#' E.g. transfer information on matches and cell body fibers, and also update lineage related information.
#'
#' @param x a vector of FAFB skeleton IDs (skids) that can be read with \code{catmaid::catmaid_skids()}
#' @param find an annotation/search term/vector of skids passed to \code{catmaid::catmaid_skids}.
#' @param ItoLee_Hemilineage character, the correct K. Ito / T. Lee hemilineage. Must be an entry in \code{\link{hemibrain_hemilineages}}. If set to 'wipe' all lineage related annotations will be removed.
#' @param delete.find logical, is \code{TRUE} then if \code{find} is an annotation, it will be wiped from the neuron after this function is used (if you have permission to remove it)
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
#' }}
#' @rdname fafb_hemibrain_annotate
#' @export
#' @seealso \code{\link{hemibrain_get_meta}}
fafb_hemibrain_annotate <- function(x, ...){

  # Get matches
  matches = hemibrain_matches()
  matches = subset(matches, match.quality %in% c("good","medium","poor"))

  # Go by neuron and relay the results
  amatches = acts = acbfs = c()
  for(i in x){
    # Get old annotations
    a = catmaid_get_annotations_for_skeletons(sk, ...)
    ct = a$annotation[grepl("Cell_type: |cell_type: ",a$annotation)]
    cbf = a$annotation[grepl("cellBodyFiber: |CellBodyFiber: ",a$annotation)]

    # remove annotions
    if(length(ct)){
      tryCatch(catmaid::catmaid_remove_annotations_for_skeletons(skids = i, annotations = paste0("annotation:^",ct,"$"), force = TRUE, ...),
               error = function(e) NULL)
    }
    if(length(cbf)){
      tryCatch(catmaid::catmaid_remove_annotations_for_skeletons(skids = i, annotations = paste0("annotation:^",cbf,"$"), force = TRUE, ...),
               error = function(e) NULL)
    }

    # Get meta for matches
    print(x)
    y = matches[x,"match"]
    hmeta = hemibrain_get_meta(y, ...)

    # Assemble annotations
    am = paste0('hemibrain_match: ', y)
    message(am)
    act = paste0('cell_type: ', hmeta[,"type"])
    message(act)
    acbf = paste0('cellBodyFiber: ', hmeta[,"cellBodyFiber"])
    message(acbf)

    # Set annotations
    catmaid::catmaid_set_annotations_for_skeletons(skids = i,
                                                   annotations = c(am,
                                                                   act,
                                                                   acbf),
                                                   force = TRUE, ...)

    # Collect for meta-annotation
    amatches = c(amatches, am)
    acts = c(act, act)
    acbfs = c(acbfs, acbf)
  }

  # Meta-annotate
  catmaid::catmaid_set_meta_annotations(meta_annotations = "cell type", annotations = acts, ...)
  catmaid::catmaid_set_meta_annotations(meta_annotations = "cell body fiber", annotations = acbfs, ...)
  catmaid::catmaid_set_meta_annotations(meta_annotations = "hemibrain match", annotations = amatches, ...)

}

#' @rdname fafb_hemibrain_annotate
#' @export
fafb_set_hemilineage <- function(find,
                                 putative = TRUE,
                                 ItoLee_Hemilineage,
                                 delete.find = FALSE,
                                 ...){
  if(ItoLee_Hemilineage!="wipe"){
    # Make sure we have a valid hemilineage
    if(length(ItoLee_Hemilineage)>1|length(ItoLee_Hemilineage)==0){
      stop("Please provide a single ItoLee_Hemilineage")
    }
    if(!ItoLee_Hemilineage%in%hemibrain_hemilineages$ItoLee_Hemilineage){
      stop("Please provide a valid ItoLee_Hemilineage. See 'hemibrain_hemilineages'")
    }
    # Find neurins based on an annotation or as skeleton IDs
    skds = catmaid::catmaid_skids(x = find, ...)
    # Get meta information
    i = ItoLee_Hemilineage
    ItoLee_Lineage = subset(hemibrain_hemilineages, hemibrain_hemilineages$ItoLee_Hemilineage == i)$ItoLee_Lineage[1]
    Hartenstein_Lineage = subset(hemibrain_hemilineages, hemibrain_hemilineages$ItoLee_Hemilineage == i)$Hartenstein_Lineage[1]
    Hartenstein_Hemilineage = subset(hemibrain_hemilineages, hemibrain_hemilineages$ItoLee_Hemilineage == i)$Hartenstein_Hemilineage[1]
    # Make new annotations
    ItoLee_Lineage = paste0("ItoLee_Lineage: " , ItoLee_Lineage)
    ItoLee_Hemilineage = paste0("ItoLee_Hemilineage: " , ItoLee_Hemilineage)
    Hartenstein_Lineage = paste0("Hartenstein_Lineage: " , Hartenstein_Lineage)
    Hartenstein_Hemilineage = paste0("Hartenstein_Hemilineage: " , Hartenstein_Hemilineage)
    # Set these as putative?
    if(putative){
      ItoLee_Lineage = paste0(ItoLee_Lineage, "_putative")
      ItoLee_Hemilineage = paste0(ItoLee_Hemilineage, "_putative")
      Hartenstein_Lineage = paste0(Hartenstein_Lineage, "_putative")
      Hartenstein_Hemilineage = paste0(Hartenstein_Hemilineage, "_putative")
      lins = c(ItoLee_Lineage,
               ItoLee_Hemilineage,
               Hartenstein_Lineage,
               Hartenstein_Hemilineage,
               "Lineage_putative")
      message("Adding annotations:   ", paste(lins, collapse = "  ", sep = "  "))
      # Set annotations
      catmaid::catmaid_set_annotations_for_skeletons(skids = skds,
                                                     annotations = lins, force = TRUE, ...)
    }else{
      # Set annotations
      lins = c(ItoLee_Lineage,
               ItoLee_Hemilineage,
               Hartenstein_Lineage,
               Hartenstein_Hemilineage,
               "Lineage_annotated")
      message("Adding annotations:   ", paste(lins, collapse = "  ", sep = "  "))
      catmaid::catmaid_set_annotations_for_skeletons(skids = skds,
                                                     annotations = lins,
                                                     force = TRUE, ...)
      catmaid::catmaid_set_meta_annotations(meta_annotations = "ItoLee_Lineage", annotations = ItoLee_Lineage,...)
      catmaid::catmaid_set_meta_annotations(meta_annotations = "ItoLee_Hemilineage", annotations = ItoLee_Hemilineage,...)
      catmaid::catmaid_set_meta_annotations(meta_annotations = "Hartenstein_Lineage", annotations = Hartenstein_Lineage,...)
      catmaid::catmaid_set_meta_annotations(meta_annotations = "Hartenstein_Hemilineage", annotations = Hartenstein_Hemilineage,...)
    }
  }

  # Delete old annotations
  a = catmaid::catmaid_get_annotations_for_skeletons(find)
  old = unique(a$annotation[grepl("ItoLee_|Hartenstein_|Volker_|volker_|Itolee_|itolee_|itoLee_",a$annotation)])
  old = setdiff(old,lins)
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


# Hidden
## Don't use unless very necessary
fafb_lineage_complete_wipe <- function(server = "v14seg-Li-190411.0", ...){
  # Choose CATMAID server
  conn = catmaid::catmaid_login(...)
  if(!is.null(server)){
    conn$server = paste0("https://neuropil.janelia.org/tracing/fafb/",server,"/")
  }
  a = catmaid::catmaid_get_annotationlist(..., conn = conn)
  instance = catmaid_login(conn=conn, ...)$server
  message("Looking at ", instance)
  old = a$annotations$name[grepl("ItoLee_|Hartenstein_|Volker_|volker_|Itolee_|itolee_|itoLee_|Lineage_",a$annotations$name)]
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
