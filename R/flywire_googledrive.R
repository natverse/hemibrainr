#' Download a large set of well-traced skeletonised neurons from FlyWire
#'
#'@description Get a large number of skeletonised neurons from FlyWire. The neuronlist is saved as a
#'  \code{nat::neuronlistfh} object so certain neurons may be read from it
#'  without loading the entire, large neuronlist into memory. You will need access to the hemibrain Google Team Drive and
#'  have it mounted with Google filestream.
#'
#' @param brain the brainspace in which hemibrain neurons have been registered. Defaults to raw voxel space for the FlyWire project.
#'@param ... Additional arguments passed to \code{\link{nlapply}} when iterating
#'  over neurons to download.
#'
#' @examples
#' \donttest{
#' \dontrun{
#' fw.neurons = flywire_neurons()
#'}}
#'@export
#'@seealso \code{\link{hemibrain_read_neurons}}
#'@importFrom utils download.file
#'@importFrom googledrive drive_ls as_id
flywire_neurons <- function(local = FALSE,
                            brain = c("FlyWire", "JRCFIB2018F","JRCFIB2018F","FAFB","JFRC2","JRC2018F","FCWB"),
                            mirror = FALSE){
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
    neurons.fh = scale_neurons(neurons.fh,scaling=8/1000)
  }

  # Return
  neurons.fh
}
