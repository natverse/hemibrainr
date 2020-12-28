###############################################################################
################################ Standardise ##################################
###############################################################################

#' Use standard names and spellings
#'
#' @description Standardise the names of lineage groups, neuron compartments and transmitters.
#'
#' @param x a character vector to be standardised.
#' @param invert return compartment numbers rather than names.
#'
#' @return a character vector
#' @export
#' @rdname standardise
standard_transmitters <- function(x){
  x[grepl("^Neurotrans|^neurotrans",x)] = "transmitter"
  x[grepl("^ACh|^Ach|^ach|^acet|^ChA|^CHa|^cholin|^ACH",x)] = "acetylcholine"
  x[grepl("^gaba|^GABA|^GAD",x)] = "GABA"
  x[grepl("^glut|^vGlut|^Vglut|^glutamate|^Glutamate|^GLUT",x)] = "glutamate"
  x[grepl("^5-HT|^5HT|^Dope|^dope|^Dopa|^dopa|^DOP",x)] = "dopamine"
  x[grepl("^Sero|^sero|^TH-|^SER",x)] = "serotonin"
  x[grepl("^Oct|^oct|^OCT",x)] = "octopamine"
  x[grepl("^Unknow|NA|unknow|^None|^none",x)] = "unknown"
  x[is.na(x)] = "unknown"
  x
}

#' @export
#' @rdname standardise
standard_statuses <- function(x, invert= FALSE){
  standard_status <-function(z, invert = FALSE){
    if(invert){
      z[z=="incomplete"] = "i"
      z[z=="complete"] = "c"
      z[z=="adequate"] = "a"
      z[z=="merge_error"] = "m"
      z[z=="needs_extending"] = "e "
      z[z=="wrong_hemilineage"] = "w"
      z[z=="wrong_side"] = "s"
      z[z=="not_neuron"] = "n"
    }else{
      z[z=="i"] = "incomplete"
      z[z=="c"] = "complete"
      z[z=="a"] = "adequate"
      z[z=="m"] = "merge_error"
      z[z=="e"] = "needs_extending "
      z[z=="w"] = "wrong_hemilineage"
      z[z=="s"] = "wrong_side"
      z[z=="n"] = "not_neuron"
    }
    paste(sort(z),collapse="/",sep="/")
  }
  y = strsplit(x=x,split="/")
  z = sapply(y,standard_status)
  x
}
# i = incomplete [very small fragment]
# c = complete [well fleshed out neuron, may even have most medium/small branches]
# a = adequate [there is a cell body fibre, axon and dendrite]
# m = noticable merge error [this neuron is merged to another]
# e = needs extending [not quite adequate, but more than a tiny fragment]
# w = wrong hemilineage [based on its soma position and cell body fibre, this neuron looks like it is not in the same hemilineage as others of this tab]
# s = wrong side [soma is on the wrong hemisphere, given the name of this tab]
# n = not a neuron [this segmentation is not a neuron, i.e. glia, erroneous]

#' @export
#' @rdname standardise
standard_lineages <- function(x){
  x[grepl("^ItoLee_l|^itolee_l|^ItoLee_L|^itolee_L",x)] = "ItoLee_Lineage"
  x[grepl("^hartenstein_l|^Hartenstein_l|^Volker_l|^volker_l|
          ^hartenstein_L|^Hartenstein_L|^Volker_L|^volker_L",x)] = "Hartenstein_Lineage"
  x[grepl("^ItoLee_h|^itolee_h",x)] = "ItoLee_Hemilineage"
  x[grepl("^hartenstein_h|^Hartenstein_h|^Volker_h|^volker_h|
          ^hartenstein_H|^Hartenstein_H|^Volker_h|^volker_H",x)] = "Hartenstein_Hemilineage"
  x[is.na(x)] = "unknown"
  x
}

#' @export
#' @rdname standardise
standard_compartments <- function(x, invert = FALSE){
  if(invert){
    x[x=="dendrite"] = 3
    x[x=="axon"] = 2
    x[x=="soma"] = 1
    x[x=="primary.dendrite"] = 4
    x[x=="primary.neurite"] = 7
  }else{
    x[x==0] = "unknown"
    x[x==3] = "dendrite"
    x[x==2] = "axon"
    x[x==1] = "soma"
    x[x==4] = "primary.dendrite"
    x[x==7] = "primary.neurite"
  }
  x
}

#' @export
#' @rdname standardise
standardise <- standardize <- function(x){
  x <- standard_transmitters(x)
  x <- standard_lineages(x)
  x <- standard_compartments(x)
  x <- standard_statuses(x)
  x
}

# hidden
standardise_quality <- function(x){
  x[x=="e"] = "good"
  x[x=="o"] = "medium"
  x[x=="p"] = "poor"
  x[x=="t"] = "tract"
  x[x=="n"] = "none"
  x
}
