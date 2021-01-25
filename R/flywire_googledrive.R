#' Download a large set of well-traced skeletonised neurons from FlyWire
#'
#'@description Get a large number of skeletonised neurons from FlyWire. The neuronlist is saved as a
#'  \code{nat::neuronlistfh} object so certain neurons may be read from it
#'  without loading the entire, large \code{neuronlist} into memory. You will need access to the hemibrain Google Team Drive and
#'  have it mounted with Google filestream.The function \code{flywire_neurons_update} can be used to update the available data.
#'  If you want to flag flywire neurons that should be added to the Google drive, without doing this yourself, you can use
#'  \code{flywire_request}. The \code{flywire_basics} function will calculate some useful meta-data, namely
#' a point in the primary neurite tract that can be used as a stable reference for this neuron.
#'
#' @param brain the brainspace in which hemibrain neurons have been registered. Defaults to raw voxel space for the FlyWire project.
#' @param x flywire IDs for desired. If left as \code{NULL}, all flywire neurons that can be summoned from the hemibrainr Google drive are summoned.
#' @param local \code{FALSE} or path. By default (\code{FALSE}) data is read from \code{options()$Drive_hemibrain_data}), but the user can specify an alternative path.
#' @param mirror logical, whether or not to read neurons that have been mirrored (i.e. flipped to the 'other' brain hemisphere).
#' @param flywire.neurons a \code{neuronlist} of flywire neurons in FlyWire space. The \code{flywire_basics} function will calculate some useful meta-data, namely
#' a point in the primary neurite tract that can be used as a stable reference for this neuron.
#' @param x flywire IDs to update, for the saved Google drive \code{neuronlistfh} objects called with \code{flywire_neurons}.
#' @param request a neuronlist, matrix of x,y,z position or flywire ID to add to a
#' \href{https://docs.google.com/spreadsheets/d/1rzG1MuZYacM-vbW7100aK8HeA-BY6dWAVXQ7TB6E2cQ/edit#gid=0}{Google sheet} that records flywire positions
#' flagged to be processed into neuron skeletons that can be called by \code{flywire_neurons}.
#' @param nblast which flywire NBLAST to update on Google drive.
#' @param selected_file the Google sheet onto which to add new flywire coordinate. I.e. \href{https://docs.google.com/spreadsheets/d/1rzG1MuZYacM-vbW7100aK8HeA-BY6dWAVXQ7TB6E2cQ/edit#gid=0}{Google sheet}.
#' @param selected_sheet the Google sheet onto which to add new flywire coordinate. I.e. \href{https://docs.google.com/spreadsheets/d/1rzG1MuZYacM-vbW7100aK8HeA-BY6dWAVXQ7TB6E2cQ/edit#gid=0}{Google sheet}.
#' @param sheet the tab onto which to add your requests.
#' @param gsheet logical, whether or not the request is are googlesheet keys. If they are, every item in their \code{flywire.xyz} columns is added to
#' the sheet specified by \code{selected_sheet}.
#' @param swc logical. When using neurons with \code{flywire_neurons} from the Google drive, whether to read \code{.swc} files (if \code{TRUE}), or a neuronlistfh object (default).
#' @param ... Additional arguments passed to \code{nat::nlapply}.and/or \code{fafbseg::skeletor}.
#'
#' @examples
#' \donttest{
#' \dontrun{
#' # Loads all processed flywire neurons as a neuronlistfh object
#' fw.neurons = flywire_neurons()
#'
#' # Get them already bridged to the JRC2018F brainsapce
#' ## (Bogovic et al. 2018, high performance brain space)
#' fw.neurons.jrc2018f = flywire_neurons(brain = "JRC2018F")
#'
#' # Get them already mirrored to the other hemispehre, i.e. flipped
#' fw.neurons.jrc2018f.m = flywire_neurons(brain = "JRC2018F", mirror = TRUE)
#'
#' # Now say you have some flywire ID~s you wants to add
#' ## to the nightly processing so they are available from Google drive.
#' ## And they are stored in a .csv:
#' library(readr)
#' csv = read_csv("/Users/abates/Downloads/FlyWire_list.txt",
#' col_type = cols(.default = "c"))
#' ids = csv[,1][[1]]
#' neurons = skeletor(ids)
#' neurons.with.info = hemibrainr:::flywire_basics(neurons)
#' new.points = nat::xyzmatrix(neurons.with.info[,"flywire.xyz"])
#' flywire_request(new.points)
#' ## Now they are added to a google sheet, and will be read and
#' ##  processed as part of this nightly pipeline:
#' ### https://github.com/flyconnectome/fafbpipeline
#'
#'}}
#'@return A \code{neuronlist} object containing flywire skeletons. In the meta-data, it might be useful for some users to note that
#'you will get:
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
flywire_neurons <- function(x = NULL,
                            local = FALSE,
                            brain = c("FlyWire", "JRCFIB2018Fraw","JRCFIB2018F","FAFB","FAFB14","JFRC2", "JFRC2013","JRC2018F","FCWB"),
                            mirror = FALSE,
                            swc = FALSE,
                            ...){
  brain = match.arg(brain)
  if(brain == "JRCFIB2018Fraw"){
    brain = "JRCFIB2018F"
    scale = TRUE
  }else if (brain %in% c("FAFB")){
    brain = "FAFB14"
    scale = FALSE
  }else{
    scale = FALSE
  }

  # Get Google drive folder
  savedir = good_savedir(local = local)
  neuronsdir = file.path(savedir,"flywire_neurons/")
  fhdir = file.path(neuronsdir,brain,"/")

  # Read
  if(swc){
    warning("Only native flywire neurons in FlyWire space supported when swc = TRUE")
    swc.folder = file.path(neuronsdir,"FlyWire","swc")
    swc.list = list.files(swc.folder, full.names = TRUE, pattern = "swc")
    if(!is.null(x)){
      swc.list = swc.list[grepl(paste(x,collapse="|"),swc.list)]
    }
    neurons.fh = nat::read.neurons(paths = swc.list)
    fw.meta = flywire_meta(sql=FALSE)
    neurons.fh[,] = fw.meta[match(names(neurons.fh),fw.meta$flywire.id),]
  }else{
    message("Loading ", fhdir)
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
  }

  # Scale neurons if needs be
  if(scale){
    neurons.fh = scale_neurons(neurons.fh,scaling=8/1000, ...)
  }

  # Return
  if(!is.null(x)){
    neurons.fh = neurons.fh[names(neurons.fh)%in%x,]
  }
  neurons.fh
}

#' @rdname flywire_neurons
#' @export
flywire_neurons_update <- function(x,
                                brain = c("FlyWire", "JRCFIB2018F","FAFB14","JFRC2", "JFRC2013", "JRC2018F","FCWB"),
                                mirror = FALSE,
                                ...){
  if(brain=="FlyWire"){
    brain = "FAFB14"
  }
  if(!requireNamespace("fafbseg", quietly = TRUE)) {
    stop("Please install fafbseg using:\n", call. = FALSE,
         "remotes::install_github('natverse/fafbseg')")
  }
  brain = match.arg(brain)
  if(!requireNamespace("fafbseg", quietly = TRUE)) {
    stop("Please install fafbseg using:\n", call. = FALSE,
         "remotes::install_github('natverse/fafbseg')")
  }
  if(mirror){
    flip = "_mirrored"
    flip = NULL
  }
  message("Skeletonising skelected neurons using skeletor ...")
  y = fafbseg::skeletor(x, brain = elmr::FAFB14.surf, ...)
  y = flywire_basics(y)
  f = flywire_neurons(brain = brain, mirror = mirror)

  # Transform
  if(mirror){
    message("Mirroring neurons version ...")
    m = nat.templatebrains::xform_brain(y, reference = "JRC2", sample = "FAFB14", ...)
    t = nat.templatebrains::mirror_brain(m, brain = "JRC2", ...)
    if(brain!="FAFB14"){
      t = nat.templatebrains::xform_brain(t, sample = "JRC2", reference = brain, ...)
    }
  }else if (brain!="FAFB14"){
    message("Transforming skeletons ...")
    t = nat.templatebrains::xform_brain(y, reference = brain, sample = "FAFB14", ...)
  }else{
    t  = y
  }

  # Add to data
  message("Making neuronlist object ...")
  z = nat::union(t, nat::as.neuronlist(f))

  # Save
  googledrive_upload_neuronlistfh(z,
                                  team_drive = hemibrainr_team_drive(),
                                  file_name = sprintf("flywire_neurons_%s%s.rds", brain, flip),
                                  folder = "flywire_neurons",
                                  subfolder = brain)
}

#' @rdname flywire_neurons
#' @export
flywire_nblast_update <- function(x = NULL,
                                  nblast = c("hemibrain-flywire",
                                             "flywire",
                                             "flywire-mirror"),
                                  ...){
  # neuronlist
  if(!nat::is.neuronlist(x)){
    stop("x must be a neuronlist!")
  }

  # Get old NBLAST
  nblast = match.arg(nblast)
  old.nblast = hemibrain_nblast(nblast)

  # Get neurons to be updated
  message("Loading flywire neurons ...")
  old.nblast = old.nblast[,setdiff(colnames(old.nblast),x)]
  f = flywire_neurons(brain= "JRCFIB2018F")
  if(is.null(x)){
    x = setdiff(names(f),colnames(nblast))
  }
  f = f[names(f) %in% as.character(x)]
  f = nat::dotprops(f, ...)
  m = flywire_neurons(brain= "JRCFIB2018F", mirror = TRUE)
  m = m[names(m) %in% as.character(x)]
  m = nat::dotprops(m, ...)
  names(m) = paste0(names(m),"_m")
  d = union(f, m)

  ### NBLAST native
  message("NBLASTing ...")
  if(nblast=="hemibrain-flywire"){
    # Get hemibrain neurons
    message("Loading hemibrain neurons ...")
    all.neurons.flow.microns.dps = hemibrain_neurons(brain= "JRCFIB2018F", dotprops = TRUE)

    # NBLAST!
    message("NBLASTing hemibrain vs flywire")
    hemibrain.flywire.mean.1 = nat.nblast::nblast(query = d,
                                                  target = all.neurons.flow.microns.dps,
                                                  .parallel=TRUE,
                                                  normalised = TRUE)
    hemibrain.flywire.mean.2 = nat.nblast::nblast(query = all.neurons.flow.microns.dps,
                                                  target = d,
                                                  .parallel=TRUE,
                                                  normalised = TRUE)
    hemibrain.flywire.mean = (hemibrain.flywire.mean.1+t(hemibrain.flywire.mean.2))/2

    ## Average native and mirrored
    colnames(hemibrain.flywire.mean) = gsub("_m$","",colnames(hemibrain.flywire.mean))
    rownames(hemibrain.flywire.mean) = gsub("_m$","",rownames(hemibrain.flywire.mean))
    hemibrain.flywire.mean = collapse_matrix_by_names(hemibrain.flywire.mean, FUN = max)
    hemibrain.flywire.mean[hemibrain.flywire.mean<-0.5]=-0.5
    hemibrain.flywire.mean=round(hemibrain.flywire.mean, digits=3)
    hemibrain.flywire.mean = merge(hemibrain.flywire.mean, old.nblast, all.x = TRUE, all.y = TRUE)

    # Save NBLAST on hemibrain googledrive
    googledrive_upload_nblast(hemibrain.flywire.mean)

  }

  # FlyWire
  if(nblast=="flywire-mirror"){
    # Get hemibrain neurons
    message("Loading hemibrain neurons ...")
    all.fw = flywire_neurons()
    all.fw.dps = nat::dotprops(all.fw, .parallel = TRUE, OmitFailures = TRUE)
    rm("all.fw")
    flywire.mirror.mean.1 = nat.nblast::nblast(query = d,
                                                         target = all.fw.dps,
                                                         .parallel=TRUE,
                                                         normalised = TRUE)
    flywire.mirror.mean.2 = nat.nblast::nblast(query = all.fw.dps,
                                                         target = d,
                                                         .parallel=TRUE,
                                                         normalised = TRUE)
    flywire.mirror.mean = (flywire.mirror.mean.1+t(flywire.mirror.mean.2))/2
    colnames(flywire.mirror.mean) = gsub("_m$","",colnames(flywire.mirror.mean))
    rownames(flywire.mirror.mean) = gsub("_m$","",rownames(flywire.mirror.mean))
    flywire.mirror.mean[flywire.mirror.mean<-0.5]=-0.5
    flywire.mirror.mean=round(flywire.mirror.mean, digits=3)
    flywire.mirror.mean = merge(flywire.mirror.mean, old.nblast, all.x = TRUE, all.y = TRUE)
    googledrive_upload_nblast(flywire.mirror.mean)
  }

  # Just FlyWire left-right
  if(nblast=="flywire"){
    # Get hemibrain neurons
    message("Loading flywire neurons ...")
    all.fw = flywire_neurons()
    all.fw.dps = nat::dotprops(all.fw, .parallel = TRUE, OmitFailures = TRUE)
    rm("all.fw")
    flywire.mirror.mean = nat.nblast::nblast(query = f,
                                                         target = all.fw.dps,
                                                         .parallel=TRUE,
                                                         normalised = TRUE)
    flywire.mean=round(flywire.mean, digits=3)
    flywire.mean[flywire.mean<-0.5]=-0.5
    flywire.mean = merge(flywire.mean, old.nblast, all.x = TRUE, all.y = TRUE)
    googledrive_upload_nblast(flywire.mean)
  }

}

#' @rdname flywire_neurons
#' @export
flywire_basics <- function(flywire.neurons, ...){
  if(!nat::is.neuronlist(flywire.neurons)){
    stop("flywire.neurons must be a neuronlist")
  }

  # Get xyz for primary branch points
  simp = nat::nlapply(flywire.neurons,nat::simplify_neuron,n=1)
  branchpoints = sapply(simp, function(y) nat::xyzmatrix(y)[ifelse(length(nat::branchpoints(y)),nat::branchpoints(y),max(nat::endpoints(y))),])
  branchpoints = t(branchpoints)
  flywire.nm.xyz = apply(branchpoints, 1, paste_coords)

  # Convert to voxel space
  branchpoints.flywire.raw = scale(branchpoints, scale = c(4, 4, 40), center = FALSE)
  flywire.xyz = apply(branchpoints.flywire.raw, 1, paste_coords)

  # Flywire svid
  svids=fafbseg::flywire_xyz2id(nat::xyzmatrix(branchpoints.flywire.raw), root=FALSE, rawcoords = TRUE)

  # Get FAFBv14 nm coordinates
  FAFB.xyz = nat.templatebrains::xform_brain(branchpoints, sample = "FlyWire", reference = "FAFB14", ...)
  FAFB.xyz = apply(FAFB.xyz, 1, paste_coords)
  #FAFB.xyz = ""

  # Add meta data
  flywire.neurons[,"flywire.id"] =  names(flywire.neurons)
  flywire.neurons[,"flywire.xyz"] = flywire.xyz
  flywire.neurons[,"flywire.svid"] = svids
  flywire.neurons[,"FAFB.xyz"] = FAFB.xyz
  flywire.neurons[,"dataset"] = "flywire"
  flywire.neurons[,"id"] = NULL

  # Add IDs
  flywire.neurons = add_field_seq(flywire.neurons,flywire.neurons[,"flywire.id"],field="flywire.id")

  # return
  flywire.neurons

}

# Add neuron to request Google sheets
#' @rdname flywire_neurons
#' @export
flywire_request <- function(request,
                            gsheet = FALSE,
                            selected_sheet = options()$flywire_flagged_gsheet,
                            sheet = "flywire",
                            ...){
  if(!requireNamespace("fafbseg", quietly = TRUE)) {
    stop("Please install fafbseg using:\n", call. = FALSE,
         "remotes::install_github('natverse/fafbseg')")
  }
  # What kind of request is it?
  if(gsheet){
    xyz = c()
    type = "googlesheet"
    for(sheet in request){
      tabs = hemibrainr:::gsheet_manipulation(FUN = googlesheets4::sheet_names,
                                              ss = sheet,
                                              return = TRUE)
      for(tab in tabs){
        fts = try(flywire_tracing_sheet(ws = tab,
                                    selected_sheet=sheet,
                                    Verbose = FALSE), silent = TRUE)
        xyz = unique(c(xyz,fts$flywire.id))
      }
      if(!length(xyz)){
        stop("No flywire positions to add")
      }
    }
  }else{
    type = if(nat::is.neuronlist(request)){
      "neuronlist"
    }else if(is.data.frame(request)|is.matrix(request)){
      "xyz"
    }else{
      if(nrow(nat::xyzmatrix(request))){
        "xyz"
      }else{
        "ids"
      }
    }
  }
  message("Request is ", type)

  # Get coordinates
  if(type=='ids'){
    request = fafbseg::skeletor(request, ...)
  }
  if(nat::is.neuronlist(request)){
    fb = flywire_basics(request)
    xyz = do.call(rbind, lapply(fb[,"flywire.xyz"], function(y) strsplit(y,",| |;|:")))
  }else if(type!="googlesheet"){
    xyz = as.data.frame(nat::xyzmatrix(request), stringsAsFactors =FALSE)
  }
  fw.xyz = paste_coords(nat::xyzmatrix(xyz))
  fw.xyz = fw.xyz[fw.xyz!="(NA,NA,NA)"]

  # Add to Google sheet
  gs = try(flywire_tracing_sheet(ws=tab,open=FALSE,selected_sheet=selected_sheet), silent = FALSE)
  if(class(gs)!="try-error"){
    fw.xyz = setdiff(fw.xyz,gs$flywire.xyz)
    update = data.frame(User = "flywire", flywire.xyz = fw.xyz)
    for(col in setdiff(colnames(gs),colnames(update))){
      update[[col]] = NA
    }
    update = update[,colnames(gs)]
  }else{
    update = data.frame(User = "flywire", flywire.xyz = fw.xyz)
  }
  batches = split(1:nrow(update), ceiling(seq_along(1:nrow(update))/500))
  for(i in batches){
    gsheet_manipulation(FUN = googlesheets4::sheet_append,
                        data = update[min(i):max(i),],
                        ss = selected_sheet,
                        sheet = sheet)
  }
  message("FlyWire positions added")
}

#' Read precomputed flywire data from the hemibrainr Google Drive
#'
#' @description Read precomputed data available on the hemibrain Google Team
#'   Drive. (see \code{\link{hemibrainr_set_drive}}) and (see
#'   \code{\link{hemibrainr_rclone}}). This includes body IDs for all flywire
#'   neurons ((\code{flywire_ids})), and user contributions towards their
#'   creation (\code{flywire_contributions}), as well as Flywire related NBLAST
#'   scores retrieved using \code{\link{hemibrain_nblast}}.
#'
#' @inheritParams hemibrainr_googledrive_data
#'
#' @return a \code{data.frame}. Depending on which synapse function was called, it can contain the columns:
#'
#' \itemize{
#'
#'   \item{"flywire.xyz"} { - coordinates of a point in the neuron in flywire voxel space. XYZ, separated by a semicolon.}
#'
#'   \item{"flywire.id"}{ - the unique ID associated with this flywire neuron. This ID changes every time a neuron is, even slightly, modified. So it is an unstable identifier.
#'   This is why \code{flywire.xyz} is sometimes used.}
#'
#'   \item{"fw.x"}{ - the x coordinate of a point in the flywire neuron, in flywire voxel space..}
#'
#'   \item{"fw.y"}{ - the y coordinate of a point in the flywire neuron, in flywire voxel space..}
#'
#'   \item{"fw.z"}{ - the z coordinate of a point in the flywire neuron, in flywire voxel space..}
#'
#'   \item{"user_name"}{ - the name of the user who made the number of edits given in this row.}
#'
#'   \item{"edits"}{ - the number of edits (merges, splits, etc.) made by a user for the given \code{flywire.id}.}
#'
#'   \item{"proportion"}{ - the proportion of total edits for this neuron, that the given user made.}
#'
#'   \item{"dataset"}{ - the dataset this neuron is from, i.e. flywire.}
#'
#'}
#'
#' @examples
#' \donttest{
#' \dontrun{
#'
#' # All flywire IDs for neurons that have a split precomputed
#' fw.ids = flywire_ids()
#'
#' # For these flywire IDs, their meya data:
#' fw.meta = flywire_meta()
#'
#' # For flywire IDs, which users contributed what:
#' fw.edits = flywire_contributions()
#'
#' }}
#' @seealso \code{\link{hemibrain_splitpoints}},
#'   \code{\link{hemibrain_flow_centrality}},
#'   \code{\link{hemibrainr_googledrive_data}},
#'   \code{\link{hemibrain_metrics}}
#' @name flywire_googledrive_data
#' @aliases flywire_meta
#' @export
flywire_meta <-function(local = FALSE, folder = "flywire_neurons/", sql = FALSE, ...){
  savedir = good_savedir(local = local)
  if(sql){
    find_gsql(savedir = savedir, tab = "flywire_meta", sql.db = "flywire_data.sqlite", folder = folder, ...)
  }else{
    gfile = find_gfile(savedir = savedir, file = "flywire_meta", folder = folder)
    gcsv = suppressWarnings(readr::read_csv(gfile, col_types = sql_col_types))
    gcsv
  }
}

#' @rdname flywire_googledrive_data
#' @export
flywire_failed <-function(local = FALSE, folder = "flywire_neurons/", sql = FALSE, ...){
  savedir = good_savedir(local = local)
  if(sql){
    find_gsql(savedir = savedir, tab = "flywire_failed", sql.db = "flywire_data.sqlite", folder = folder, ...)
  }else{
    gfile = find_gfile(savedir = savedir, file = "flywire_failed", folder = folder)
    gcsv = suppressWarnings(readr::read_csv(gfile, col_types = sql_col_types))
    gcsv
  }
}

# hidden
sql_col_types = readr::cols(.default = "c",
                 edits = "i",
                 total.edits = "i",
                 pre = "i",
                 post = "i",
                 upstream = "i",
                 downstream = "i",
                 voxels = "i",
                 layer = "n",
                 ct.layer = "n",
                 orig.soma = "?",
                 soma = "?",
                 soma.edit = "?",
                 truncated = "?",
                 splittable = "?",
                 checked = "?",
                 edited.cable = "?",
                 skeletonization = "?",
                 time = "?",
                 total.outputs = "n",
                 axon.outputs = "n",
                 dend.outputs = "n",
                 total.inputs = "n",
                 axon.inputs = "n",
                 dend.inputs = "n",
                 axon.outputs = "n",
                 dend.outputs = "n",
                 axon.inputs = "n",
                 dend.inputs = "n",
                 total.outputs.density = "n",
                 total.inputs.density = "n",
                 axon.outputs.density = "n",
                 dend.outputs.density = "n",
                 axon.inputs.density = "n",
                 dend.inputs.density = "n",
                 total.length = "n",
                 axon.length = "n",
                 dend.length = "n",
                 pd.length = "n",
                 length = "?",
                 segregation_index = "n",
                 X = "n",
                 Y = "n",
                 Z = "n",
                 x = "n",
                 x = "n",
                 z = "n",
                 fw.x ="n",
                 fw.y = "n",
                 fw.z = "n",
                 cable.length = "n",
                 proportion = "n",
                 count = "n",
                 weight = "n",
                 norm = "n",
                 synapses = "n",
                 syn = "n",
                 syns = "n",
                 n.syn = "n",
                 nodes = "n",
                 endpoints = "n",
                 segments = "n",
                 branchpoints = "n",
                 root = "i")

#' @rdname flywire_googledrive_data
#' @export
flywire_contributions <-function(local = FALSE, folder = "flywire_neurons/", sql = TRUE, ...){
  savedir = good_savedir(local = local)
  if(sql){
    find_gsql(savedir = savedir, tab = "flywire_edits", sql.db = "flywire_data.sqlite", folder = folder, ...)
  }else{
    gfile = find_gfile(savedir = savedir, file = "flywire_edits", folder = folder)
    gcsv = suppressWarnings(readr::read_csv(gfile, col_types = sql_col_types))
    gcsv
  }
}

#' @rdname flywire_googledrive_data
#' @export
flywire_ids <-function(local = FALSE, folder = "flywire_neurons/", sql = FALSE, ...){
  savedir = good_savedir(local = local)
  if(sql){
    find_gsql(savedir = savedir, tab = "flywire_ids", sql.db = "flywire_data.sqlite", folder = folder, ...)
  }else{
    gfile = find_gfile(savedir = savedir, file = "flywire_ids", folder = folder)
    gcsv = suppressWarnings(readr::read_csv(gfile, col_types = sql_col_types))
    as.character(gcsv$x)
  }
}

#' Update the flywire.id column in a set of google sheets based on flywire xyz positions
#'
#' @description This function retrieves flywire IDs based on xyz positions in flywire voxel space, from a set of google sheets.
#' It also writes the updated flywire IDs to the same google sheets. This is often helpful because flywire IDs are inherently unstable, they change every time
#' a neuron is modified even slightly. Users can record 'stable' points in a neuron that identify it, e.g. a single xyz position in the cell body fibre, or at the soma, and
#' then use this function to update and get the correct flywire ID whenever they wish.
#'
#' @param selected_sheets character vector. the google sheet(s) to update. Each entry is a unique google sheet ID. You can find these in a sheet's URL.
#' If \code{NULL} then defaults to \code{option('hemibrainr_gsheets')}.
#' @param chosen.columns as well as writing column updates to the specified google sheets, this function returns a \code{data.frame} built from all given sheets and their
#' individual tabs, that have been updated. This argument specifies which column you want returned. Filled with NAs if it does not exist.
#' @param ws character vector, tabs (i.e. work sheets) on the google sheet to query/read. This works with regex, so you only need to give the name partially.
#' If set to \code{NULL} for \code{flywire_tracing_sheets}, the whole google sheet is read and all tabs are combined using \code{plyr::rbind.fill}.
#' @param regex logical, use \code{ws} with regex.
#' @param match logical. If \code{TRUE}, hemibrain matches given.
#' @param meta meta data for flywire neurons, e.g. as retreived using \code{\link{flywire_meta}}. Used to efficiently input \code{flywire.xyz} column if only a \code{flywire.id} entry has been given.
#' Only works if that id is also in this provided \code{data.frame}, \code{meta}.
#' @inheritParams matches_update
#'
#' @details For this function to work, the specified google sheet(s) must have either the column \code{flywire.xyz},
#' which gives the xyz position of points in a format that can be read by \code{nat::xyzmatrix}, for example \code{"(135767,79463,5284)"} or \code{"(135767;79463;5284)"}.
#' If this is missing, then the columns: \code{fw.x}, \code{fw.y}, \code{fw.z} must be specified. The xyz positions must be in FlyWire voxel space, which is what you get if you use the
#' copy location tool in the flywire.ai web-interface.
#'
#' The logic of the update procedure is:, find the \code{flywire.xyz} column.
#' If that does not exist, find: \code{fw.x}, \code{fw.y}, \code{fw.z}, and use that to create a \code{flywire.xyz} column.
#' We use \code{flywire.xyz} if both are given, and there is a mismatch.
#' For each row, a \code{flywire.id} is then found based on these points, using \code{fafbseg::flywire_xyz2id} (using the argument \code{rawcoords = TRUE}).
#' The google sheet columns \code{flywire.id},\code{flywire.xyz}, \code{fw.x}, \code{fw.y}, \code{fw.z} are then updated if they exist in the original google sheet.
#' If they do not, they are not updated. The function returns a \code{data.frame} combining all tabs of all googlesheets specified, but returning only the columns
#' specified by the argument \code{chosen.columns}.
#'
#' @return a \code{data.frame} with columns from the given google sheet(s), specified using the argument \code{chosen.columns}.
#'
#' @examples
#' \donttest{
#' \dontrun{
#'
#' # Update flywire.ids in the sheet:
#' ### https://docs.google.com/spreadsheets/d/
#' ### 1rzG1MuZYacM-vbW7100aK8HeA-BY6dWAVXQ7TB6E2cQ
#' ### /edit#gid=0
#' ## And return
#' fw.gsheet.meta = flywire_ids_update(selected_sheets =
#' "1rzG1MuZYacM-vbW7100aK8HeA-BY6dWAVXQ7TB6E2cQ")
#'
#' }}
#' @seealso \code{\link{flywire_ids}},
#'   \code{\link{flywire_meta}},
#'   \code{\link{flywire_neurons}}
#' @name flywire_ids_update
#' @export
flywire_ids_update <- function(selected_sheets = NULL, # "1rzG1MuZYacM-vbW7100aK8HeA-BY6dWAVXQ7TB6E2cQ"
                               chosen.columns = c('flywire.xyz', "flywire.svid",
                                                  "flywire.id", "skid",
                                                  "FAFB.xyz", "cell.type", "side",
                                                  "ItoLee_Hemilineage", "Hartenstein_Hemilineage",
                                                  "status"),
                               ws = NULL,
                               regex = FALSE,
                               meta = NULL,
                               match = FALSE,
                               matching_sheet = options()$hemibrainr_matching_gsheet,
                               priority = c("FAFB","hemibrain")){
  if(is.null(selected_sheets)){
    selected_sheets = getOption("hemibrainr_gsheets", stop("Please set option('hemibrainr_gsheets')"))
  }
  fw.columns = c("flywire.id","fw.x","fw.y","fw.z","flywire.xyz","flywire.svid")
  gs = data.frame(stringsAsFactors = FALSE)
  for(selected_sheet in selected_sheets){
    ## Read Google sheets and extract glywire neuron positions
    message("### Working on google sheet: ", selected_sheet)
    if(is.null(ws)){
      tabs = gsheet_manipulation(FUN = googlesheets4::sheet_names,
                                 ss = selected_sheet,
                                 return = TRUE)
    }else{
      if(regex){
        ws=regex_tab_names(regex=ws,selected_sheet=selected_sheet)$name
      }
      tabs=ws
    }
    pb = progress::progress_bar$new(
      format = "  updating tab :what [:bar] :percent eta: :eta",
      clear = FALSE, total = length(tabs))
    for(tab in tabs){
      pb$tick(tokens = list(what = tab))
      gs.t = gs.t.current = flywire_tracing_sheet(ws = tab, selected_sheet=selected_sheet, Verbose = FALSE)
      used.cols = colnames(gs.t.current)
      if(nrow(gs.t)&&ncol(gs.t)&&sum(grepl("fw.x|flywire.xyz",used.cols))>0){
        # Separate x,y,z positions
        for(fw.c in fw.columns){
          if(is.null(gs.t[[fw.c]])){
            gs.t[[fw.c]] = NA
          }
        }
        good.xyz = sapply(gs.t$flywire.xyz,function(x) length(tryCatch(nat::xyzmatrix(x),error = function(e) NA))==3)
        gs.t[!good.xyz,c("flywire.xyz")] = apply(gs.t[!good.xyz,c("fw.x","fw.y",'fw.z')],1,paste_coords)
        if(!is.null(meta)){
          justids = gs.t$flywire.xyz==paste_coords(matrix(NA,ncol=3))&!is.na(gs.t$flywire.id)
          replacement.xyz = meta[match(gs.t[justids,"flywire.id"],meta$flywire.id),]$flywire.xyz
          gs.t[justids,"flywire.xyz"] = replacement.xyz
          needs.updating = is.na(replacement.xyz)
          please.update = gs.t$flywire.id[justids][needs.updating]
          if(length(please.update)){
            updated.ids = tryCatch(sapply(please.update, fafbseg::flywire_latestid), error = function(e) NULL)
            if(!is.null(updated.ids)){
              gs.t$flywire.id[justids][needs.updating] = updated.ids
              good.xyz = sapply(gs.t$flywire.xyz,function(x) length(tryCatch(nat::xyzmatrix(x),error = function(e) NA))==3)
              gs.t[!good.xyz,c("flywire.xyz")] = apply(gs.t[!good.xyz,c("fw.x","fw.y",'fw.z')],1,paste_coords)
              justids = gs.t$flywire.xyz==paste_coords(matrix(NA,ncol=3))&!is.na(gs.t$flywire.id)
              replacement.xyz = meta[match(gs.t[justids,"flywire.id"],meta$flywire.id),]$flywire.xyz
              # replacement.xyz[is.na(replacement.xyz)] = "(NA,NA,NA)"
              gs.t[justids,"flywire.xyz"] = replacement.xyz
            }
          }
        }
        good.xyz = sapply(gs.t$flywire.xyz,function(x) length(tryCatch(nat::xyzmatrix(x),error = function(e) NA))==3)
        if(sum(good.xyz)){
          gs.t[good.xyz,c("fw.x","fw.y","fw.z")] = nat::xyzmatrix(gs.t[good.xyz,"flywire.xyz"])
        }
        bad.svids = (is.na(gs.t$flywire.svid)|gs.t$flywire.svid=="0")&!is.na(gs.t$flywire.xyz)
        if(sum(bad.svids)>0){
          gs.t$flywire.svid[bad.svids] = fafbseg::flywire_xyz2id(nat::xyzmatrix(gs.t$flywire.xyz[bad.svids]),
                                                                 root=FALSE,
                                                                 rawcoords = TRUE)
        }
        if(!all(is.na(gs.t$flywire.xyz))){
          # Get flywire IDs from these positions
          bbx = matrix(c(5100, 1440, 16, 59200, 29600, 7062),ncol=3,byrow = TRUE)
          bbx = nat::boundingbox(scale(bbx, scale = 1/c(4, 4, 40), center = FALSE))
          gs.t$fw.x = as.numeric(gs.t$fw.x)
          gs.t$fw.y = as.numeric(gs.t$fw.y)
          gs.t$fw.z = as.numeric(gs.t$fw.z)
            pos = gs.t[apply(gs.t, 1, function(row) sum(is.na(row[c("fw.x","fw.y",'fw.z')]))==0),]
            p = nat::pointsinside(pos[,c("fw.x","fw.y",'fw.z')],bbx)
            pos = pos[p,]
            pos = pos[!is.na(pos$fw.x),]
            if(nrow(pos)){
              foreach.ids = fafbseg::flywire_xyz2id(pos[,c("fw.x","fw.y",'fw.z')], rawcoords = TRUE)
              names(foreach.ids) = pos[,"flywire.xyz"]
            }else{
              foreach.ids = NULL
            }
          fids = unlist(foreach.ids)
          fids[is.na(fids)|is.nan(fids)] = "0"
          replacement = fids[match(gs.t$flywire.xyz ,names(fids))]
          coordsmissing = gs.t$flywire.xyz==paste_coords(matrix(NA,ncol=3))
          coordsmissing[is.na(gs.t$flywire.xyz)] = TRUE
          coordsmissing[is.na(replacement)] = TRUE
          replacement[coordsmissing] = gs.t$flywire.id[coordsmissing]
          gs.t$flywire.id = replacement
          gs.t$flywire.xyz[gs.t$flywire.xyz==paste_coords(matrix(NA,ncol=3))] = NA
          if(nrow(gs.t)!=nrow(gs.t.current)){
            stop("Sheet processing corruption.")
          }
          write.cols = intersect(fw.columns,used.cols)
          if(length(fids) & length(write.cols)){
            gsheet_update_cols(
              write.cols = write.cols,
              gs=gs.t[,used.cols],
              selected_sheet = selected_sheet,
              sheet = tab,
              Verbose = FALSE)
          }
        }
        # Now continue processing
        gs.t = gs.t[,colnames(gs.t)%in%chosen.columns]
        for(col in chosen.columns){
          if(is.null(gs.t[[col]])){
            gs.t[[col]] = NA
          }
        }
        if(nrow(gs.t)){
          gs.t$gsheet = selected_sheet
          gs = plyr::rbind.fill(gs.t,gs)
        }
      }
    }
  }
  # Make this unique, but keep row with most information
  master = gs[!is.na(gs$flywire.xyz),]
  master$filled = apply(master, 1, function(r) sum(!is.na(r)))
  master = master[order(master$filled,decreasing = TRUE),]
  master = master[!duplicated(master$flywire.xyz),]
  rownames(master) = master$flywire.xyz
  if(match){
    matches = hemibrain_matches(selected_file = matching_sheet, priority = priority)
    matches = subset(matches, matches$match.dataset == "hemibrain" &  matches$dataset == "flywire")
    master$hemibrain.match = matches[match(master$flywire.id, matches$id), "match"]
    master$hemibrain.match.quality = matches[match(master$flywire.id, matches$id),"quality"]
    master$FAFB.hemisphere.match = matches[match(master$flywire.id, matches$id), "FAFB.hemisphere.match"]
    master$FAFB.hemisphere.match.quality = matches[match(master$flywire.id, matches$id),"FAFB.hemisphere.match.quality"]
    master$cell.type = matches[match(master$flywire.id, matches$id), "cell.type"]
  }
  master = subset(master, !is.na(master$flywire.id))
  master$filled = NULL
  master
}

#' Update neuron match information on a google sheet
#'
#' @description This function retrieves neuron matches (\code{\link{hemibrain_matches}}) from a master-matching google sheet. If then
#' can update other google sheets, specified by the user to update neuron-match information.
#' Just columns giving the match, match quality and cell type are updated.
#'
#' @param selected_sheets character vector. the google sheet(s) to update. Each entry is a unique google sheet ID. You can find these in a sheet's URL.
#' If \code{NULL} then defaults to \code{option('hemibrainr_gsheets')}.
#' @param matching_sheet the master matching sheet. Cannot be in \code{selected_sheets}. This sheet will be processed using \code{\link{hemibrain_matches}} to produce a
#' \code{data.frame} describing hemibrain-FAFB and FAFB-hemibrain matches.
#' @param priority whether to use FAFB->hemibrain matches (FAFB) or hemibrain->FAFB matches (hemibrain) in order to ascribe
#' cell type names to FAFB neurons. In both cases, cell type names are attached to hemibrain bodyids, and propagated to their FAFB matches.
#' @param id the ID specifying unique neuron identity for each row, for each tab of the google sheets specified with \code{selected_sheets}. When matches are added
#' they are for the neuron specified in the \code{id} column of each worksheet.
#' @param match.field which match to record. E.g. if \code{id} is \code{"flywire.id"} you may want to add the hemibrain match from the master matching sheets. To do this
#' set \code{match.field} to \code{"hemibrain"}. Then, if there is a \code{"hemibrain.match"} and/or \code{"hemibrain.match.quality"} column, it will be updated. You
#' may also want to know the hemispheric match within FAFB, in which case \code{"FAFB.hemisphere"} could be used.
#' @param chosen.columns as well as writing column updates to the specified google sheets, this function returns a \code{data.frame} built from all given sheets and their
#' individual tabs, that have been updated. This argument specifies which column you want returned. Filled with NAs if it does not exist.
#'
#' @details For this function to work, the specified google sheet(s) must have either the column specified with the argument \code{id},
#' which gives the ID of a neuron from the FAFB (CATMAID or flywire) data set, or the hemibrain data set.
#' Each tab of each specified google sheet is read in turn. The data set match (specified with the argument \code{match.field}). A match is recovered
#' using \code{hemibrain_matching(selected_file = matching_sheet, priority=priority)}. If the right columns exist in the google sheet, i.e. 'hemibrain.match'
#' and 'hemibrain.match.quality' for a hemibrain match, then this column in the google sheet is wiped and updated. If not, no update takes place.
#'
#' @return a \code{data.frame} with columns from the given google sheet(s), specified using the argument \code{chosen.columns}.
#'
#' @examples
#' \donttest{
#' \dontrun{
#'
#' # Update flywire.ids in the sheet:
#' ## https://docs.google.com/spreadsheets/d/
#' ## 1spGSuhUX6Hhn-8HH0U_ArIWUuPpMBFNjIjeSSh_MFVY
#' ## /edit#gid=603762040
#' lineage.gsheet.meta = matches_update(id = "flywire.id",
#' match.field = "hemibrain",
#' selected_sheets = "1spGSuhUX6Hhn-8HH0U_ArIWUuPpMBFNjIjeSSh_MFVY")
#'
#' }}
#' @seealso \code{\link{hemibrain_matches}},
#'   \code{\link{hemibrain_matched}},
#'   \code{\link{hemibrain_matching}}
#' @name matches_update
#' @export
matches_update <- function(matching_sheet = options()$hemibrainr_matching_gsheet,
                           priority = c("FAFB","hemibrain"),
                           selected_sheets,
                           id = c("flywire.id","bodyid","skid","id"),
                           match.field = c("hemibrain","CATMAID","flywire","LM","FAFB.hemisphere"),
                           chosen.columns = c("cell.type", "flywire.xyz", "side",
                                              "ItoLee_Hemilineage", "Hartenstein_Hemilineage")){
  id = match.arg(id)
  match.field = match.arg(match.field)
  priority = match.arg(priority)
  matches = hemibrain_matches(selected_file = matching_sheet, priority = priority)
  matches = subset(matches, ! matches$quality %in% "NBLAST")
  if(sum(selected_sheets%in%c(options()$hemibrainr_matching_gsheet,matching_sheet))>1){
    stop("selected_sheets should not indicate master matching sheets.
         I.e. in options()$hemibrainr_matching_gsheet.")
  }
  # iterate through sheets and tabs
  if(is.null(selected_sheets)){
    selected_sheets = getOption("hemibrainr_gsheets", stop("Please set option('hemibrainr_gsheets')"))
  }
  gs = data.frame(stringsAsFactors = FALSE)
  for(selected_sheet in selected_sheets){
    ## Read Google sheets and extract glywire neuron positions
    tabs = gsheet_manipulation(FUN = googlesheets4::sheet_names,
                               ss = selected_sheet,
                               return = TRUE)
    for(tab in tabs){
      gs.t = gs.t.current = flywire_tracing_sheet(ws = tab, selected_sheet=selected_sheet)
      used.cols = colnames(gs.t.current)
      field = match.field
      if(field=="CATMAID"){ # for historic reasons. Chould fix.
        field = "FAFB"
      }
      match.column = paste0(field,".match")
      quality.column = paste0(field,".match.quality")
      if(!match.column%in%used.cols){
        match.column = paste0(field,"_match")
        if(!quality.column%in%used.cols){
          quality.column = paste0(field,"_match.quality")
        }
      }
      chosen.columns = unique(c(chosen.columns,match.column,quality.column,"cell.type"))
      in.sheet = intersect(gs.t[[id]],matches$id)
      if(nrow(gs.t) & length(in.sheet)){
        matches.sel = subset(matches, matches$id%in%in.sheet & matches$match.dataset == match.field)
        if(match.field=="hemibrain"){
          gs.t[[match.column]] = matches.sel$match[match(gs.t[[id]], matches.sel$id)]
          gs.t[[quality.column]] = matches.sel$quality[match(gs.t[[id]], matches.sel$id)]
        }else if(match.field == "FAFB.hemisphere"){
          gs.t[[match.column]] = matches.sel$FAFB.hemisphere.match[match(gs.t[[id]], matches.sel$id)]
          gs.t[[quality.column]] = matches.sel$FAFB.hemisphere.match.quality[match(gs.t[[id]], matches.sel$id)]
        }else if(match.field == "LM"){
          gs.t[[match.column]] = matches.sel$LM.match[match(gs.t[[id]], matches.sel$id)]
          gs.t[[quality.column]] = matches.sel$LM.match.quality[match(gs.t[[id]], matches.sel$id)]
        }else if(match.field == "CATMAID"){
          gs.t[[match.column]] = matches.sel$match[match(gs.t[[id]], matches.sel$id)]
          gs.t[[quality.column]] = matches.sel$quality[match(gs.t[[id]], matches.sel$id)]
        }
        gs.t[["cell.type"]] = matches.sel$cell.type[match(gs.t[[id]], matches.sel$id)]
        if(nrow(gs.t)!=nrow(gs.t.current)){
          stop("Sheet processing corruption.")
        }
        write.cols = intersect(c(match.column,quality.column,"cell.type"),used.cols)
        if(nrow(matches.sel) & length(write.cols)){
          gsheet_update_cols(
            write.cols = write.cols,
            gs=gs.t[,used.cols],
            selected_sheet = selected_sheet,
            sheet = tab)
        }
        gs.t = gs.t[,colnames(gs.t)%in%chosen.columns]
        for(col in chosen.columns){
          if(is.null(gs.t[[col]])){
            gs.t[[col]] = NA
          }
        }
        gs = plyr::rbind.fill(gs,gs.t)
      }
    }
  }
  gs
}


