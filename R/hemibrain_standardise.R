###############################################################################
################################ Standardise ##################################
###############################################################################

#' Use standard names and spellings
#'
#' @description Standardise the names of lineage groups and transmitters.
#'
#' @param x a character vector to be standardised.
#'
#' @return a character vector
#' @export
#' @rdname standardise
standard_transmitters <- function(x){
  x[grepl("^Neurotrans|^neurotrans",x)] = "transmitter"
  x[grepl("^ACh|^Ach|^acet|^ChA|^CHa|^cholin",x)] = "acetylcholine"
  x[grepl("^gaba|^GABA|^GAD",x)] = "GABA"
  x[grepl("^glut|^vGlut|^Vglut|^glutamate|^Glutamate",x)] = "glutamate"
  x[grepl("^5-HT|^5HT|^Dope|^dope|^Dopa|^dopa",x)] = "dopamine"
  x[grepl("^Sero|^sero|^TH-",x)] = "serotonin"
  x[grepl("^Oct|^oct",x)] = "octopamine"
  x[grepl("^Unknow|NA|unknow",x)] = "unknown"
  x[is.na(x)] = "unknown"
  x
}

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
standard_compartments <- function(x){
  x[x==3] = "dendrite"
  x[x==2] = "axon"
  x[x==1] = "soma"
  x[x==4] = "primary.dendrite"
  x[x==7] = "primary.neurite"
  x
}

#' @export
#' @rdname standardise
standardise <- standardize <- function(x){
  x <- standard_transmitters(x)
  x <- standard_lineages(x)
  x <- standard_compartments(x)
  x
}
