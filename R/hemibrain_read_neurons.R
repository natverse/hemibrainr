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
#' @param savedir directory in which \code{nat::neuronlistfh} object as been saved using \code{hemibrain_download_neurons}. This acts as a database,
#' from which hemibrain neurons can be read.
#' @param microns convert dimensions from raw voxels into microns (template brain: \code{JRCFIB2018F}, else \code{JRCFIB2018Fraw}).
#' @param ... methods passed to \code{neuprintr::neuprint_read_neurons}, \code{\link{hemibrain_remove_bad_synapses}}
#'  and \code{\link{hemibrain_flow_centrality}}
#'
#' @inherit flow_centrality return
#'
#' @examples
#' \donttest{
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
#' \dontrun{
#' # Plot the split to check it
#' nat::nopen3d()
#' nlscan_split(neurons, WithConnectors = TRUE)
#' }}
#' @export
#' @seealso \code{\link{hemibrain_splitpoints}}, \code{\link{hemibrain_flow_centrality}},
#' \code{\link{hemibrain_precomputed_splitpoints}}, \code{\link{hemibrain_metrics}},\code{\link{hemibrain_remove_bad_synapses}}
hemibrain_read_neurons<-function(x = NULL,
                                 savedir = FALSE,
                                 microns = TRUE,
                                 ...){
  if(is.null(x)&!savedir){
    stop("You must either supply bodyids to read from neuPrint using x, or
         specify a location from which to read saved a save neuronlistfh object using savedir. See
         ?hemibrain_download_neurons for details on the latter option.")
  }
  if(isTRUE(savedir)){
    savedir = options()$hemibrain_data
    if(length(list.files(path = savedir, pattern = ".rds"))){
      neurons.flow = hemibrain_read_neurons_local(savedir = savedir)
      neurons.flow = neurons.flow.fh[as.character(x)]
    }
  }else{
    neurons = neuprintr::neuprint_read_neurons(x, ...)
    neurons.flow = hemibrain_flow_centrality(neurons, ...)
    neurons.flow = hemibrain_remove_bad_synapses(neurons.flow, ...)
  }
  neurons.flow = hemibrain_clean_skeleton(neurons.flow, rval = "neuron", ...)
  hemibrain_metrics = hemibrain_metrics[,!colnames(hemibrain_metrics)%in%colnames(neurons.flow[,])]
  df = cbind(neurons.flow[,], hemibrain_metrics[names(neurons.flow),])
  rownames(df) = names(neurons.flow)
  if(microns){
    neurons.flow = scale_neurons(neurons.flow, scaling = (8/1000))
    nat.templatebrains::regtemplate(neurons.flow) = "JRCFIB2018F"
  }else{
    nat.templatebrains::regtemplate(neurons.flow) = "JRCFIB2018Fraw"
  }
  neurons.flow[,] = df
  neurons.flow = add_field_seq(neurons.flow,neurons.flow[,"bodyid"],field="bodyid")
  neurons.flow
}

# hidden
scale_neurons <-function(x, scaling = (8/1000), ...) UseMethod("scale_neurons")
scale_neurons.neuron <- function(x, scaling, ...){
  nat::xyzmatrix(x$d) = nat::xyzmatrix(x$d)*scaling
  nat::xyzmatrix(x$connectors) = nat::xyzmatrix(x$connectors)*scaling
  x
}
scale_neurons.neuronlist = function(x, scaling = (8/1000), ...){
  nat::nlapply(x,scale_neurons.neuron, scaling = scaling, ...)
}

#' Remove erroneous branchlets to beautify neuron skeletons
#'
#' @description Remove synapse-less branches from axons and dendrites, and clean up
#' the primary neurite and primary dendrite (linker) of a neuron. Givne neurons must have been split
#' by either \code{\link{flow_centrality}} or \code{\link{hemibrain_flow_centrality}}.
#' \code{\link{hemibrain_flow_centrality}} to re-root and split neurons into putative
#' axons and dendrites. Optionally, it may also convert neurons from their raw voxel
#' space (as they are stored in NeuPrint) to microns.
#'
#' @param x a \code{nat::neuronlist} or \code{nat::neuron} object
#' @param rval whether to return a pruned neuron, a neuron with the synapse-less branches
#' set to \code{Label = 0} or the \code{PointNo} for those brnahces.
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
  a = axonic_cable(x)
  a.clean = prune_synapseless_branches(a, neuron = FALSE)
  d = dendritic_cable(x)
  d.clean = prune_synapseless_branches(d, neuron = FALSE)
  remove = c(a.clean, d.clean)
  if(rval=="pruned"){
    if(length(remove)>0){
      y = nat::prune_vertices(x, verticestoprune = remove, invert = FALSE)
      y$connectors = x$connectors[x$connectors$treenode_id %in% y$d$PointNo, ]
      relevant.points = subset(x$d, PointNo %in% y$d$PointNo)
      y$d = relevant.points[match(y$d$PointNo, relevant.points$PointNo), ]
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
prune_synapseless_branches <- function(x, neuron = TRUE){
  s = x$SubTrees
  prune = c()
  con.pos = unique(x$connectors$treenode_id)
  con.pos = x$d$PointNo[match(con.pos,x$d$PointNo)]
  if(is.null(con.pos)){
    stop("No connectors in neuron,")
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
      relevant.points = subset(x$d, PointNo %in% y$d$PointNo)
      y$d = relevant.points[match(y$d$PointNo, relevant.points$PointNo), ]
    }else{
      y=x
    }
    y = hemibrain_neuron_class(y)
    y
  }else{
    x$d$PointNo[prune]
  }
}


#'Download all automatically split neurons from the hemibrain project
#'
#'@description Download all the automatically split neurons from the hemibrain
#'  project, from Google Drive. The neuronlist is saved as a
#'  \code{nat::neuronlistfh} object so certain neurons may be read from it
#'  without loading the entire, large neuronlist into memory. In order to do
#'  this, you need to add the following folder to your Google Drive:
#'  https://drive.google.com/drive/folders/14UPg7CvHDtvzNjvNgAULYnxZ018Xgf5H?usp=sharing.
#'   Contact us if you do not have, but would like, permission.
#'
#'@param savedir where to save the \code{.rds} and meta data files for the
#'  \code{nat::neuronlistfh} object
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
#' neurons = hemibrain_read_neurons("1702323386", savedir = TRUE)
#'}}
#'@export
#'@seealso \code{\link{hemibrain_splitpoints}},
#'  \code{\link{hemibrain_flow_centrality}},
#'  \code{\link{hemibrain_precomputed_splitpoints}},
#'  \code{\link{hemibrain_metrics}},\code{\link{hemibrain_remove_bad_synapses}}
#'@importFrom utils download.file
#'@importFrom googledrive drive_ls as_id
hemibrain_download_neurons <- function(savedir = TRUE,
                                       overwrite = FALSE,
                                       ...
                                       ){
  savedir = good_savedir(savedir=savedir)
  message(sprintf("If this takes to much time, you can also download the relevant Google Drive folder manually. To do so, open this link: %s and then download the folder to this location on your computer: ",
          "https://drive.google.com/open?id=1px6o2R_heFLRCtTF4Q2SvJNR9CWFhud0. Remember to unzip all files.",
          savedir))
  ### Download file hash meta data
  message("Finding data ...")
  datals = tryCatch(googledrive::drive_ls(googledrive::as_id("https://drive.google.com/drive/folders/1AHgTOzoYfyhsd6XimVuKq1jgLeyVBpWd")),
                    error = function(e) NULL)
  if(is.issue(datals$id[1])){
    utils::browseURL("https://drive.google.com/drive/folders/14UPg7CvHDtvzNjvNgAULYnxZ018Xgf5H?usp=sharing",
                     browser = getOption("browser"),
                     encodeIfNeeded = FALSE)
    stop("Google Drive download failed. This might be because you have not added our data folder to your drive.
            Please follow this link: https://drive.google.com/drive/folders/14UPg7CvHDtvzNjvNgAULYnxZ018Xgf5H?usp=sharing,
            and add this folder to your Google Drive. If you do not have permission, you can contact us for access.")
  }
  ### Download neuronlist as a .rds file
  message("Downloading metadata ...")
  nlfh = paste0(savedir,"hemibrain_all_neurons_flow_polypre_centrifugal_synapses/hemibrain_all_neurons_flow_polypre_centrifugal_synapses.rds")
  download.file(sprintf("https://drive.google.com/uc?authuser=0&id=%s&export=download","19mpATVptH9n42yhmZybEYl9k9PndN4_G"),
                destfile = nlfh)
  ### Download actual data
  message("Downloading data, ", nrow(datals)," files ...")
  metafh = paste0(savedir,"hemibrain_all_neurons_flow_polypre_centrifugal_synapses/data/")
  download = googledrive_downloadmany(ls = datals, dir = metafh, overwrite = overwrite, OmitFailures = TRUE, ...)
  message(length(unlist(download)),"/",nrow(datals)," data files successfully downloaded. Re-run to try remainder.")
  message("Download complete, files can be found here: ", savedir)
}

# hidden
## download many google drive files
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

# hidden
hemibrain_read_neurons_local <- function(savedir = TRUE,
                                         neuron.split = "hemibrain_all_neurons_flow_polypre_centrifugal_synapses"){
  savedir = good_savedir(savedir = savedir,
                      neuron.split = neuron.split)
  if(length(list.files(path = savedir, pattern = ".rds"))){
    neurons.flow = nat::read.neuronlistfh(savedir,neuron.split,"/",neuron.split,".rds")
  }else{
    stop("neuronlistfh (.rds) file not found at: ", savedir)
  }
  neurons.flow
}

# hidden
good_savedir <- function(savedir = TRUE,
                         neuron.split = "hemibrain_all_neurons_flow_polypre_centrifugal_synapses"
                         ){
  if(isTRUE(savedir)){
    savedir = options()$hemibrain_data
  }
  if(is.issue(savedir)){
    options(hemibrain_data = paste0(getwd(),"/data-raw/hemibrain_data/"))
    warning("The following option has been set: options(hemibrain_data = paste0(getwd(),'/data-raw/hemibrain_data/')) ")
  }
  if(!dir.exists(savedir)){
    dir.create(savedir, recursive = TRUE)
    warning("Made new hemibrain save directory. Neurons destined to be saved in ", savedir)
  }
  suppressWarnings(dir.create(paste0(savedir,"/",neuron.split,"/data/"), recursive = TRUE))
  savedir
}


