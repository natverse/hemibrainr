#' Download a large set of well-traced skeletonised neurons from FlyWire
#'
#'@description Get a large number of skeletonised neurons from FlyWire. The neuronlist is saved as a
#'  \code{nat::neuronlistfh} object so certain neurons may be read from it
#'  without loading the entire, large neuronlist into memory. You will need access to the hemibrain Google Team Drive and
#'  have it mounted with Google Filestream.
#'
#' @param brain the brainspace in which hemibrain neurons have been registered. Defaults to raw voxel space for the FlyWire project.
#' @param local logical, whether to try to read locally saved neurons (by default at: \code{options()$hemibrain_data}) or neurons from Google Drive (\code{options()$Gdrive_hemibrain_data}).
#' @param mirror logical, whether or not to read neurons that have been mirrored (i.e. flipped to the 'other' brain hemisphere).
#' @param ... Additional arguments passed to \code{nat::\link{nlapply}}.
#'
#' @examples
#' \donttest{
#' \dontrun{
#' fw.neurons = flywire_neurons()
#' }}
#' @return A \code{neuronlist} object containing flywire skeletons. In the meta-data, it might be useful for some users to note that
#' you will get:
#'
##' \itemize{
##'  \item{"flywire.id"}{ The ID given to the corresponding volumetric body in flywire.
##'  These are used to do things like fetch volumes and are the input to the \code{skeletor} function. However, they are highly volatile and
##'  change a lot with active tracing.}
##'  \item{"flywire.xyz"}{ The voxel coordinate of a single point in this neuron, usually a cell body fiber tract position. This is a more accurate way
##'  of keeping tract of neuron as it will always correspond to the same 'neuron' even though its related flywire.id will change with merge/split events in flywire.}
##'  \item{"hemilineage"}{ An estimated hemilineage identity from both of two naming systems, Ito et al. 2013 and Wong et al. 2013}
##'  \item{"side"}{ An estimate as to the 'side', i.e. brain hemisphere, in which the neuron lies}
##'  \item{"skid"}{ The 'skeleton ID' of this neuron's match in CATMAID for FAFBv14}
##'  \item{"FAFB.xyz"}{ The coordinates in nanometres of a point in the neuron, in FAFBv14 space}
##'  \item{"hemibrain_match"}{ The bodyid of an estimated hemibrain match}
##' }
#'
#'@export
#'@seealso \code{\link{hemibrain_read_neurons}}
#'@importFrom utils download.file
#'@importFrom googledrive drive_ls as_id
flywire_neurons <- function(local = FALSE,
                            brain = c("FlyWire", "JRCFIB2018F","JRCFIB2018F","FAFB","JFRC2","JRC2018F","FCWB"),
                            mirror = FALSE,
                            ...){
  brain = match.arg(brain)
  if(brain == "JRCFIB2018Fraw"){
   scale = TRUE
  }else{
    scale = FALSE
  }

  # Get google drive folder
  savedir = good_savedir(local = local)
  neuronsdir = paste0(savedir,"flywire_neurons/")
  fhdir = paste0(neuronsdir,brain,"/")

  # Read
  filelist = list.files(path = fhdir, pattern = ".rds", full.names = TRUE)
  filelist = filelist[grepl("mirror",filelist)==mirror]
  filelist = sort(filelist,decreasing = TRUE)
  if(length(filelist)){
    fh.file = filelist[1]
    neurons.fh = nat::read.neuronlistfh(fh.file)
    attr(neurons.fh,"df") = neurons.fh[,]
  }else{
    warning("neuronlistfh (.rds) file not found at: ", fhdir)
    return(NULL)
  }

  # Scale neurons if needs be
  if(scale){
    neurons.fh = scale_neurons(neurons.fh,scaling=8/1000, ...)
  }

  # Return
  neurons.fh
}
