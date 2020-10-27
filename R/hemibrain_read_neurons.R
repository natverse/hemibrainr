################################################################################
################################ Read neurons ##################################
################################################################################

#' Read neurons from the hemibrain connectome project
#'
#' @description Read neurons from the hemibrain connectome project. This function
#' uses the package \code{neuprintr} to read neurons. It then uses
#' \code{\link{hemibrain_flow_centrality}} to re-root and split neurons into putative
#' axons and dendrites. Optionally, it may also convert neurons from their raw voxel
#' space (as they are stored in neuPrint) to microns.
#'
#' @param x a vector of bodyids that can be read from \url{'https://neuprint.janelia.org/'}.
#' @param brain the brainspace in which hemibrain neurons have been registered. Defaults to raw voxel space for the hemibrain project.
#' @param mirror logical, whether or not to read neurons that have been mirrored (i.e. flipped to the 'other' brain hemisphere).
#' @param microns convert dimensions from raw voxels into microns (template brain: \code{JRCFIB2018F}, else \code{JRCFIB2018Fraw}).
#' @param remove.bad.synapses whether or not to run \code{\link{hemibrain_remove_bad_synapses}} on neurons pulled from neuPrint.
#' @param reroot logical, whether or not somas should be re-rooted.
#' Note that if FALSE, re-rooting occurs anyway via \code{hemibrain_flow_centrality}. However, setting this argument to \code{TRUE}
#' @param googlesheet logical, whether or not manually checked somas should be read from the \href{https://docs.google.com/spreadsheets/d/1YjkVjokXL4p4Q6BR-rGGGKWecXU370D1YMc1mgUYr8E/edit#gid=1524900531}{Google Sheet}
#' @param clean whether or not to set synapse-less branches to \code{Label = 0}.
#' @param local \code{FALSE} or path. By default (\code{FALSE}) data is read from \code{options()$Drive_hemibrain_data}), but the user can specify an alternative path.
#' @param scaling the factor by which neuron coordinates in raw voxel space should be multiplied. The default scales to microns.
#' @param dotprops logical. Whether or not to retrieve a \code{nat::dotprops} object, i.e. vector cloud representations of neurons for NBLASTing.
#' The dotprops object will be in JRC2018FIBF micron space.
#' @param folder the sub-folder in which to look for a \code{nat::neuronlistfh} .rds file and its data, representing hemibrain neurons. The function will look for this
#' folder in the location: \code{hemibrainr:::good_savedir(local=local)},
#' by default the mounted Google drive (\code{options()$Gdrive_hemibrain_data}) or locally ((\code{options()$hemibrain_data}))
#' @param ... arguments passed to \code{neuprintr::neuprint_read_neurons}, \code{\link{hemibrain_remove_bad_synapses}}
#'  and \code{\link{hemibrain_flow_centrality}}
#'
#' @inherit flow_centrality return
#'
#' @examples
#' \donttest{
#' \dontrun{
#'
#' # Choose neurons
#' ## In this case some antennal lobe local neurons
#' al.local.neurons = c("1702323386", "2068966051", "2069311379", "1702305987", "5812996027",
#' "1702336197", "1793744512", "1976565858", "2007578510", "2101339904",
#' "5813003258", "2069647778", "1947192569", "1883788812", "1916485259",
#' "1887177026", "2101348562", "2132375072", "2256863785", "5813002313",
#' "5813054716", "5813018847", "5813055448", "1763037543", "2101391269",
#' "1794037618", "5813018729", "2013333009")
#'
#' # Get neurons
#' neurons = hemibrain_read_neurons(al.local.neurons)
#'
#' # Plot the split to check it
#' nat::nopen3d()
#' nlscan_split(neurons, WithConnectors = TRUE)
#' }}
#' @rdname hemibrain_read_neurons
#' @export
#' @seealso \code{\link{hemibrain_splitpoints}}, \code{\link{hemibrain_flow_centrality}}, \code{\link{hemibrain_somas}},
#' \code{\link{hemibrain_precomputed_splitpoints}}, \code{\link{hemibrain_metrics}},\code{\link{hemibrain_remove_bad_synapses}}
#' ,\code{\link{hemibrain_get_meta}},\code{\link{flywire_neurons}}
hemibrain_read_neurons<-function(x = NULL,
                                 local = FALSE,
                                 microns = FALSE,
                                 reroot = TRUE,
                                 googlesheet = FALSE,
                                 remove.bad.synapses = FALSE,
                                 clean = FALSE,
                                 ...){
  if(is.null(x) && isFALSE(local)){
    stop("You must either supply bodyids to read from neuPrint using x, or
         specify a location from which to read saved a save neuronlistfh object using 'local'. See
         ?hemibrain_download_neurons for details on the latter option.")
  }
  neurons.flow.fh = tryCatch(hemibrain_neurons(local = local), function(e) NULL)
  if(!is.null(neurons.flow.fh)){
    y = intersect(x,names(neurons.flow.fh))
    z = setdiff(x,names(neurons.flow.fh))
    if(length(z)>0){
      warning("The following bodyids could not be read from your saved neuronlistfh object. If they exist, you may need to read them from neuPrint by changing the savedir argument of this function to FALSE: ",
              paste(z,collapse=", "))
    }
    neurons.flow = neurons.flow.fh[as.character(y)]
  }else{
    neurons = neuprintr::neuprint_read_neurons(x, ...)
    if(reroot){
      neurons = hemibrain_reroot(neurons, method = "manual", googlesheet = googlesheet, ...)
    }
    neurons.flow = hemibrain_flow_centrality(neurons, ...)
    if(remove.bad.synapses){
      hemibrain.rois = hemibrain_roi_meshes()
      neurons.flow = hemibrain_remove_bad_synapses(neurons.flow, meshes = hemibrain.rois, ...)
    }
  }
  if(clean){
    neurons.flow = hemibrain_clean_skeleton(neurons.flow, rval = "neuron", ...)
  }
  if(microns){
    neurons.flow = scale_neurons(neurons.flow, scaling = (8/1000))
    nat.templatebrains::regtemplate(neurons.flow) = "JRCFIB2018F"
  }else{
    nat.templatebrains::regtemplate(neurons.flow) = "JRCFIB2018Fraw"
  }
  neurons.flow[,] = hemibrain_get_meta(names(neurons.flow))
  neurons.flow = add_field_seq(neurons.flow,neurons.flow[,"bodyid"],field="bodyid")
  neurons.flow = metadata_add_tags(neurons.flow)
  neurons.flow
}

#' @rdname hemibrain_read_neurons
#' @export
scale_neurons <-function(x, scaling = (8/1000), ...) UseMethod("scale_neurons")
#' @export
scale_neurons.neuron <- function(x, scaling = (8/1000), ...){
  nat::xyzmatrix(x$d) = nat::xyzmatrix(x$d)*scaling
  if(!is.null(x$d$W)){
    gt0=x$d$W > 0   # special case negative diameter is sometimes used as a signalling value
    x$d$W[gt0] = x$d$W[gt0]*scaling
  }
  if(!is.null(x$connectors)){
    nat::xyzmatrix(x$connectors) = nat::xyzmatrix(x$connectors)*scaling
  }
  x
}
#' @export
scale_neurons.neuronlist = function(x, scaling = (8/1000), ...){
  nat::nlapply(x,scale_neurons.neuron, scaling = scaling, ...)
}

#' Remove erroneous branchlets to beautify neuron skeletons
#'
#' @description Remove synapse-less branches from axons and dendrites, and clean up
#' the primary neurite and primary dendrite (linker) of a neuron. Given neurons must have been split
#' by either \code{\link{flow_centrality}} or \code{\link{hemibrain_flow_centrality}}.
#' \code{\link{hemibrain_flow_centrality}} to re-root and split neurons into putative
#' axons and dendrites. Optionally, it may also convert neurons from their raw voxel
#' space (as they are stored in NeuPrint) to microns.
#'
#' @param x a \code{nat::neuronlist} or \code{nat::neuron} object
#' @param rval whether to return a pruned neuron, a neuron with the synapse-less branches
#' set to \code{Label = 0} or the \code{PointNo} for those branches.
#' @param ... methods passed to \code{nat::nlapply}.
#'
#' @inherit flow_centrality return
#'
#' @examples
#' \donttest{
#' # Choose neurons
#' ## In this case some antennal lobe local neurons
#' al.local.neurons = c("1702323386", "2068966051", "2069311379", "1702305987", "5812996027",
#' "1702336197", "1793744512", "1976565858", "2007578510", "2101339904",
#' "5813003258", "2069647778", "1947192569", "1883788812", "1916485259",
#' "1887177026", "2101348562", "2132375072", "2256863785", "5813002313",
#' "5813054716", "5813018847", "5813055448", "1763037543", "2101391269",
#' "1794037618", "5813018729", "2013333009")
#'
#' # Only use first three neurons for documentation purposes
#' al.local.neurons=al.local.neurons[1:3]
#' # Get neurons
#' neurons = neuprintr::neuprint_read_neurons(al.local.neurons)
#'
#' # Split  neurons
#' neurons = hemibrain_flow_centrality(neurons)
#'
#' # Clean neurons
#' neurons.cleaned = hemibrain_clean_skeleton(neurons)
#'
#' \dontrun{
#' # Plot the split to check it
#' nat::nopen3d()
#' nlscan_split(neurons.cleaned, WithConnectors = TRUE)
#' }}
#' @export
#' @seealso \code{\link{hemibrain_splitpoints}}, \code{\link{hemibrain_flow_centrality}},
#' \code{\link{hemibrain_precomputed_splitpoints}}, \code{\link{hemibrain_metrics}},\code{\link{hemibrain_remove_bad_synapses}}
hemibrain_clean_skeleton <-function(x,
                                    rval = c("pruned","neuron","points"),
                                    ...) UseMethod("hemibrain_clean_skeleton")

#' @export
hemibrain_clean_skeleton.neuron <- function(x, rval = c("pruned","neuron","points"), ...){
  rval = match.arg(rval)
  a = tryCatch( axonic_cable(x), error = function(e) NULL )
  a.clean = tryCatch( prune_synapseless_branches(a, neuron = FALSE), error = function(e) NULL )
  d = tryCatch( dendritic_cable(x), error = function(e) NULL )
  d.clean = tryCatch( prune_synapseless_branches(d, neuron = FALSE), error = function(e) NULL )
  remove = c(a.clean, d.clean)
  if(rval=="pruned"){
    if(length(remove)>0){
      y = nat::prune_vertices(x, verticestoprune = remove, invert = FALSE)
      y = carryover_tags(x, y)
      y = carryover_labels(x, y)
    }else{
      y = x
    }
    y = hemibrain_neuron_class(y)
    y
  }else if(rval == "neuron"){
    x = add_Label(x, PointNo = remove, Label = 0, erase = FALSE)
  }else{
    remove
  }
}

#' @export
hemibrain_clean_skeleton.neuronlist <- function(x, rval = c("pruned","neuron","points"), ...){
  nat::nlapply(x, FUN = hemibrain_clean_skeleton.neuron, rval = rval, ...)
}

# hidden
# al.local.neurons="1702323386"
# neurons = neuprintr::neuprint_read_neurons(al.local.neurons)
# neurons.split = hemibrain_flow_centrality(neurons)
# neurons.cleaned = hemibrain_clean_skeleton(neurons.split)
prune_synapseless_branches <- function(x, neuron = TRUE){
  s = x$SubTrees
  prune = c()
  con.pos = unique(x$connectors$treenode_id)
  con.pos = x$d$PointNo[match(con.pos,x$d$PointNo)]
  if(is.null(con.pos)){
    stop("No connectors in neuron")
  }
  for(t in 1:x$nTrees){
    ss = unlist(s[[t]])
    ss = x$d$PointNo[ss]
    with.syns = sum(ss%in%con.pos)>0
    if(!with.syns){
      prune = c(prune, ss)
    }
  }
  if(neuron){
    if(length(prune)>0){
      y = nat::prune_vertices(x, verticestoprune = prune, invert = FALSE)
      y$connectors = x$connectors[con.pos %in% y$d$PointNo, ]
      relevant.points = subset(x$d, x$d$PointNo %in% y$d$PointNo)
      y$d = relevant.points[match(y$d$PointNo, relevant.points$PointNo), ]
    }else{
      y=x
    }
    y = hemibrain_neuron_class(y)
    y
  }else{
    x$d$PointNo[x$d$PointNo%in%prune]
  }
}


#'Download all automatically split neurons from the hemibrain project
#'
#'@description Download all the automatically split neurons from the hemibrain
#'  project, from Google Drive. The neuronlist is saved as a
#'  \code{nat::neuronlistfh} object so certain neurons may be read from it
#'  without loading the entire, large neuronlist into memory. Note, you do not
#'  need to do this if you have your Google Drive mounted on your computer and you
#'  have access to the Flyconnectome Repository 2020lhfib. In which case, you should
#'  be able to read neurons directly from the drive whiled you have an Internet connection
#'  using \code{hemibrain_read_neurons(savedir=TRUE)}. Neurons downloaded into your working directory.
#'
#'  In order to do download all neurons locally, you need to add the following folder to your Google Drive:
#'  \url{https://drive.google.com/drive/folders/14UPg7CvHDtvzNjvNgAULYnxZ018Xgf5H?usp=sharing}.
#'  Contact us if you do not have, but would like, permission.
#'  If this takes to much time, you can also download the relevant Google Drive folder manually. To do so, open this link: %s and then download the folder to this location on your computer: ",
#'  \url{https://drive.google.com/open?id=1px6o2R_heFLRCtTF4Q2SvJNR9CWFhud0}. Remember to unzip all files.
#'
#'@param overwrite whether or not to overwrite the existing saved files, within
#'  \code{savedir}.
#'@param ... Additional arguments passed to \code{\link{nlapply}} when iterating
#'  over neurons to download.
#'
#' @examples
#' \donttest{
#' \dontrun{
#' # Parallelise
#' numCores <- parallel::detectCores()
#' doMC::registerDoMC(numCores/2)
#' message("Using ", numCores/2, " cores")
#'
#' # Download neurons
#' hemibrain_download_neurons(.parallel = TRUE, overwrite = FALSE)
#'
#' # Get specific neurons
#' neurons = hemibrain_read_neurons("1702323386", savedir = TRUE, local = TRUE)
#'}}
#'@export
#'@seealso \code{\link{hemibrain_read_neurons}}
#'@importFrom utils download.file
#'@importFrom googledrive drive_ls as_id
hemibrain_download_neurons <- function(overwrite = FALSE,
                                       ...
                                       ){
  message("NOTE: This function is experimental. The recommended method of reading all hemibrain neurons quickly is to set up Google Drive using Google Filestream and use the function: hemibrain_neurons. You must have access to the Google Team Drive 'hemibrain', administered by the flyconnectome group at the University of Cambidge, to do this. If in doubt over access, you can contact us.")
  savedir = good_savedir(local = TRUE)
  message(sprintf("If this takes too much time, you can also download the relevant Google Drive folder manually. To do so, open this link: %s and then download the folder to this location on your computer: ",
          "https://drive.google.com/open?id=1px6o2R_heFLRCtTF4Q2SvJNR9CWFhud0. Remember to unzip all files.",
          paste0(savedir,"hemibrain_neurons/")))
  ### Download file hash meta data
  message("Finding data ...")
  datals = tryCatch(googledrive::drive_ls(googledrive::as_id("https://drive.google.com/drive/folders/1AHgTOzoYfyhsd6XimVuKq1jgLeyVBpWd")),
                    error = function(e) NULL)
  if(is.issue(datals$id[1])){
    utils::browseURL("https://drive.google.com/drive/folders/14UPg7CvHDtvzNjvNgAULYnxZ018Xgf5H?usp=sharing",
                     browser = getOption("browser"),
                     encodeIfNeeded = FALSE)
    stop("Google Drive download failed. This might be because you have not added our data folder to your drive. Please follow this link:","
    https://drive.google.com/drive/folders/14UPg7CvHDtvzNjvNgAULYnxZ018Xgf5H?usp=sharing",
            "\nand add this folder to your Google Drive. If you do not have permission, you can contact us for access.")
  }
  ### Download neuronlist as a .rds file
  message("Downloading metadata ...")
  nlfh = paste0(savedir,"hemibrain_neurons/polypre_centrifugal_synapses/hemibrain_all_neurons_flow_polypre_centrifugal_synapses.rds")
  download.file(sprintf("https://drive.google.com/uc?authuser=0&id=%s&export=download","19mpATVptH9n42yhmZybEYl9k9PndN4_G"),
                destfile = nlfh)
  ### Download actual data
  message("Downloading data, ", nrow(datals)," files ...")
  metafh = paste0(savedir,"hemibrain_neurons/polypre_centrifugal_synapses/data/")
  download = googledrive_downloadmany(ls = datals, dir = metafh, overwrite = overwrite, OmitFailures = TRUE, ...)
  message(length(unlist(download)),"/",nrow(datals)," data files successfully downloaded. Re-run (overwrite = FALSE) to try remainder.")
  message("Download complete, files can be found here: ", savedir)
}

# hidden
## download many Google drive files
googledrive_downloadmany <- function(ls, dir, overwrite = FALSE, ...){
  nat::nlapply(1:nrow(ls),
               function(i) googledrive_simpledownload(id=ls$id[i],
                                                      file = paste0(dir,ls$name[i]),
                                                      overwrite = overwrite),
               ...)
}

# hidden
googledrive_simpledownload <- function(id, file, overwrite = FALSE){
  if(!overwrite){
    if(!file.exists(file)){
      download.file(sprintf("https://drive.google.com/uc?authuser=0&id=%s&export=download",id),
                    destfile = file, quiet = TRUE)
    }
  }else{
    download.file(sprintf("https://drive.google.com/uc?authuser=0&id=%s&export=download",id),
                  destfile = file, quiet = TRUE)
  }
  file
}

#' @rdname hemibrain_read_neurons
#' @export
hemibrain_neurons <- function(local = FALSE,
                              brain = c("JRCFIB2018Fraw","JRCFIB2018F","FAFB14","JFRC2","JRC2018F","FCWB"),
                              mirror = FALSE,
                              dotprops = FALSE,
                              folder = "hemibrain_neurons/"){
  brain = match.arg(brain)
  savedir = good_savedir(local = local)
  neuronsdir = paste0(savedir,folder)
  if(dotprops){
    message("Vector cloud object only available for JRCFIB2018F")
    brain = "JRCFIB2018F"
    fhdir = paste0(neuronsdir,brain,"/dotprops/")
  }else{
    fhdir = paste0(neuronsdir,brain,"/")
  }
  filelist = list.files(path = fhdir, pattern = ".rds", full.names = TRUE)
  filelist = filelist[grepl("mirror",filelist)==mirror]
  filelist = sort(filelist,decreasing = TRUE)
  if(length(filelist)){
    fh.file = filelist[1]
    neurons.flow.fh = nat::read.neuronlistfh(fh.file)
  }else{
    warning("neuronlistfh (.rds) file not found at: ", fhdir)
    return(NULL)
  }
  neurons.flow.fh
}

# hidden
good_savedir <- function(local=FALSE){
  if(isFALSE(local)){
    if(is.null(options()$Gdrive_hemibrain_data)){
      options(Gdrive_hemibrain_data = "/Volumes/GoogleDrive/Shared\ drives/hemibrain/")
    }
    savedir = options()$hemibrain_data
    if(dir.exists(options()$Gdrive_hemibrain_data) & !local){
      message("Using Google Team Drive: ", hemibrainr_team_drive())
      savedir = options()$Gdrive_hemibrain_data
    }
  }else{
    message("Looking for data in ", local)
    savedir = local
  }
  savedir
}
