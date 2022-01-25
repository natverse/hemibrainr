# Lineage-related light-level data

#' Download light-level lineage clone data from Ito et al. 2013 or Yu et al. 2013
#'
#'@description Get registered light-level clonal data of central brain, secondary developmental lineages. Most neurons in the fly brain hail from one of these lineages.
#'You can find their meta information in \code{\link{hemibrain_hemilineages}}.
#'
#' @param x a specific lineage, e.g. ALlv1. If \code{NULL}, the default, all lineages are fetched
#' @param brain the brainspace in which hemibrain neurons have been registered. Defaults to raw voxel space for the FlyWire project.
#' @param dataset whether to get clonal data from Ito et al. 2013 \code{"Ito"} or Yu et al. 2013 \code{"Lee"}.
#' @param local \code{FALSE} or path. By default (\code{FALSE}) data is read from \code{options()$Drive_hemibrain_data}), but the user can specify an alternative path.
#' @param from the name system in which you can find \code{x}. One of: cellBodyFiber, ito_lee_lineage, ito_lee_hemilineage, hartenstein_lineage, hartenstein_hemilineage
#' @param to the lineage name system you want to move to, can be multiple of: cellBodyFiber, ito_lee_lineage, ito_lee_hemilineage, hartenstein_lineage, hartenstein_hemilineage
#'
#' @references Ito M, Masuda N, Shinomiya K, Endo K, Ito K.
#' Systematic analysis of neural projections reveals clonal composition of the Drosophila brain.
#' Curr Biol. 2013 Apr 22;23(8):644-55. doi: 10.1016/j.cub.2013.03.015. Epub 2013 Mar 28. PMID: 23541729.
#'
#' @references Yu HH, Awasaki T, Schroeder MD, Long F, Yang JS, He Y, Ding P, Kao JC, Wu GY, Peng H, Myers G, Lee T.
#' Clonal development and organization of the adult Drosophila central brain. Curr Biol. 2013 Apr 22;23(8):633-43.
#' doi: 10.1016/j.cub.2013.02.057. Epub 2013 Mar 28. PMID: 23541733; PMCID: PMC3637848.
#'
#' @examples
#' \donttest{
#' \dontrun{
#' # Loads all processed flywire neurons as a neuronlistfh object
#' lee.lins = lineage_clones(brain = "FAFB14", dataset = "Lee")
#' ito.lins = lineage_clones(brain = "FAFB14", dataset = "Ito")
#'
#' # Set the scene
#' library(elmr)
#' nopen3d()
#' plot3d(FAFB)
#'
#' # Plot WEDd1
#' lee.lins.wedd1 = subset(lee.lins, ito_lee_lineage == "WEDd1")
#' plot3d(xyzmatrix(lee.lins.wedd1),
#' add = TRUE, col = hemibrain_bright_colors["green"])
#'
#' # Plot WEDd1 from the other clonal data set
#' ito.lins.wedd1 = subset(ito.lins, ito_lee_lineage == "WEDd1")
#' plot3d(xyzmatrix(ito.lins.wedd1),
#' add = TRUE, col = hemibrain_bright_colors["darkgreen"])
#'
#' # Get flywire neurons
#' fw.neurons = flywire_neurons()
#' fw.neurons.wedd1 = subset(fw.neurons, ito_lee_hemilineage == "WEDd1")
#' plot3d(fw.neurons.wedd1,
#' col = hemibrain_bright_colors["marine"], soma = 4000)
#'
#' # What are the equivalent hemibrain neurons?
#' cbf.hl = subset(hemibrain_hemilineages, ito_lee_lineage == "WEDd1")
#' hb.neurons = hemibrain_neurons(brain="FAFB14")
#' hb.neurons.wedd1 = subset(hb.neurons,cellBodyFiber%in%cbf.hl$cellBodyFiber)
#' plot3d(hb.neurons.wedd1,
#' col = hemibrain_bright_colors["orange"], soma = 4000)
#'
#'}}
#'@return A \code{neuronlist} object containing point clouds for neuron lineage clones. In the meta-data, it might be useful for some users to note that
#' you will get:
#'
##' \itemize{
##'  \item{"hartenstein_lineage"}{ The Hartenstein lineage name based on expert review of light level data from the V. Hartenstein group (Wong 2013, Lovick 2013).
#'    See \code{\link{hemibrain_hemilineages}}.}
##'  \item{"ito_lee_lineage"}{ Lineage name, data from the K. Ito and T. Lee groups (Yu 2013, Ito 2013).
#'    See \code{\link{hemibrain_hemilineages}}}
##'  \item{"side"}{ The hemisphere of the brain onto which this clone is registered}
##' }
#'@name lineage_clones
#'@export
#'@seealso \code{\link{hemibrain_neurons}}, \code{\link{flywire_neurons}}, \code{\link{hemibrain_hemilineages}}
lineage_clones <- function(x = NULL,
                           local = FALSE,
                           brain = c("FAFB14", "FlyWire", "JRCFIB2018Fraw","JRCFIB2018F","FAFB","JFRC2", "JFRC2013","JRC2018F","FCWB"),
                           dataset = c("Ito","Lee")
                           ){
  brain = match.arg(brain)
  dataset = match.arg(dataset)
  if(brain == "JRCFIB2018Fraw"){
    brain = "JRCFIB2018F"
    scale = TRUE
  }else if (brain %in% c("FAFB","FlyWire")){
    brain = "FAFB14"
    scale = FALSE
  }else{
    scale = FALSE
  }

  # Get Google drive folder
  savedir = good_savedir(local = local)
  neuronsdir = file.path(savedir,"light_level/lineages/")
  fhdir = file.path(neuronsdir,dataset,brain)

  # Exists?
  if(!file.exists(fhdir)){
    stop("Cannot find file: ", fhdir)
  }

  # Read
  message("Loading ", fhdir)
  filelist = list.files(path = fhdir, pattern = ".rds", full.names = TRUE)
  filelist = sort(filelist,decreasing = TRUE)
  if(length(filelist)){
    fh.file = filelist[1]
    neurons.fh = nat::read.neuronlistfh(fh.file)
    test = tryCatch(neurons.fh[[1]], error = function(e){
      warning(e)
      try(file.remove(paste0(attributes(neurons.fh)$db@datafile,"___LOCK")), silent = FALSE)
    })
    attr(neurons.fh,"df") = neurons.fh[,]
  }else{
    warning("neuronlistfh (.rds) file not found at: ", fhdir)
    return(NULL)
  }

  # Scale neurons if needs be
  if(scale){
    neurons.fh = scale_neurons(neurons.fh,scaling=8/1000)
  }

  # Return
  if(!is.null(x)){
    neurons.fh = neurons.fh[names(neurons.fh)%in%x,]
  }
  neurons.fh
}

#'@export
#'@rdname lineage_clones
lineage_converter <- function(x, from = c("ito_lee_hemilineage", "cellBodyFiber", "ito_lee_lineage", "hartenstein_lineage", "hartenstein_hemilineage"), to = 'hartenstein_hemilineage'){
  from = match.arg(from)
  df = hemibrainr::hemibrain_hemilineages[match(x, hemibrainr::hemibrain_hemilineages[[from]]),colnames(hemibrainr::hemibrain_hemilineages)%in%to]
  as.vector(df)
}

