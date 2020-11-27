#' Download a large set of well-traced skeletonised neurons from FlyWire
#'
#'@description Get a large number of skeletonised neurons from FlyWire. The neuronlist is saved as a
#'  \code{nat::neuronlistfh} object so certain neurons may be read from it
#'  without loading the entire, large \code{neuronlist} into memory. You will need access to the hemibrain Google Team Drive and
#'  have it mounted with Google filestream.The function \code{flywire_neurons_update} can be used to update the available data.
#'  If you want to flag flywire neurons that should be added to the Google drive, without doing this yourself, you can use
#'  \code{flywire_request}.
#'
#' @param brain the brainspace in which hemibrain neurons have been registered. Defaults to raw voxel space for the FlyWire project.
#' @param local \code{FALSE} or path. By default (\code{FALSE}) data is read from \code{options()$Drive_hemibrain_data}), but the user can specify an alternative path.
#' @param mirror logical, whether or not to read neurons that have been mirrored (i.e. flipped to the 'other' brain hemisphere).
#' @param x flywire IDs to update, for the saved Google drive \code{neuronlistfh} objects called with \code{flywire_neurons}.
#' @param request a neuronlist, matrix of x,y,z position or flywire ID to add to a
#' \href{https://docs.google.com/spreadsheets/d/1rzG1MuZYacM-vbW7100aK8HeA-BY6dWAVXQ7TB6E2cQ/edit#gid=0}{Google sheet} that records flywire positions
#' flagged to be processed into neuron skeletons that can be called by \code{flywire_neurons}.
#' @param nblast which flywire NBLAST to update on Google drive.
#' @param selected_file the Google sheet onto which to add new flywire coordinate. I.e. \href{https://docs.google.com/spreadsheets/d/1rzG1MuZYacM-vbW7100aK8HeA-BY6dWAVXQ7TB6E2cQ/edit#gid=0}{Google sheet}.
#' @param sheet the tab onto which to add your requests.
#' @param ... Additional arguments passed to \code{nat::nlapply}.and/or \code{fafbseg::skeletor}.
#'
#' @examples
#' \donttest{
#' \dontrun{
#' fw.neurons = flywire_neurons()
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
flywire_neurons <- function(local = FALSE,
                            brain = c("FlyWire", "JRCFIB2018Fraw","JRCFIB2018F","FAFB","FAFB14","JFRC2", "JFRC2013","JRC2018F","FCWB"),
                            mirror = FALSE,
                            ...){
  brain = match.arg(brain)
  if(brain == "JRCFIB2018Fraw"){
    brain = "JRCFIB2018F"
    scale = TRUE
  }else if (brain %in% c("FAFB","FAFB14")){
    brain = "FlyWire"
    scale = FALSE
  }else{
    scale = FALSE
  }

  # Get Google drive folder
  savedir = good_savedir(local = local)
  neuronsdir = file.path(savedir,"flywire_neurons/")
  fhdir = file.path(neuronsdir,brain,"/")

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
  }else{
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


# hidden
flywire_basics <- function(x){
  if(!nat::is.neuronlist(x)){
    stop("x must be a neuronlist")
  }

  # Get xyz for primary branch points
  simp = nat::nlapply(x,nat::simplify_neuron,n=1)
  branchpoints = sapply(simp, function(y) nat::xyzmatrix(y)[ifelse(length(nat::branchpoints(y)),nat::branchpoints(y),max(nat::endpoints(y))),])
  branchpoints = t(branchpoints)
  flywire.xyz = apply(branchpoints, 1, paste_coords)

  # Get FAFBv14 nm coordinates
  # roots.flywire.raw = scale(roots, scale = 1/c(4, 4, 40), center = FALSE)
  #FAFB.xyz = nat.templatebrains::xform_brain(roots.flywire.raw, sample = "FlyWire", reference = "FAFB14")
  #FAFB.xyz = apply(FAFB.xyz, 1, paste_coords)
  FAFB.xyz = ""

  # Add
  x[,"flywire.id"] =  names(x)
  x[,"flywire.xyz"] = flywire.xyz
  x[,"FAFB.xyz"] = FAFB.xyz
  x[,"dataset"] = "flywire"
  x[,"id"] = NULL
  x

}

# Add neuron to request Google sheets
#' @rdname flywire_neurons
#' @export
flywire_request <- function(request,
                            selected_file = options()$flywire_flagged_gsheet,
                            sheet = "flywire",
                            ...){
  if(!requireNamespace("fafbseg", quietly = TRUE)) {
    stop("Please install fafbseg using:\n", call. = FALSE,
         "remotes::install_github('natverse/fafbseg')")
  }
  # What kind of request is it?
  type = if(nat::is.neuronlist(request)){
    "neuronlist"
  }else if(is.data.frame(request)|is.matrix(request)){
    if(nrow(nat::xyzmatrix(request))){
      "xyz"
    }else{
      "ids"
    }
  }else{
    "ids"
  }
  message("Request is ", type)

  # Get coordinates
  if(type=='ids'){
    request = fafbseg::skeletor(request, ...)
  }
  if(nat::is.neuronlist(request)){
    fb = flywire_basics(request)
    xyz = do.call(rbind, lapply(fb[,"flywire.xyz"], function(y) strsplit(y,",| |;|:")))
  }else{
    xyz = as.data.frame(nat::xyzmatrix(request), stringsAsFactors =FALSE)
  }

  # Add to Google sheet
  batches = split(1:nrow(xyz), ceiling(seq_along(1:nrow(xyz))/500))
  for(i in batches){
    gsheet_manipulation(FUN = googlesheets4::sheet_append,
                        data = xyz[min(i):max(i),],
                        ss = selected_file,
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
flywire_meta <-function(local = FALSE, folder = "flywire_neurons/", sql = TRUE, ...){
  savedir = good_savedir(local = local)
  if(sql){
    find_gsql(savedir = savedir, tab = "flywire_meta", sql.db = "flywire_data.sqlite", folder = folder, ...)
  }else{
    gfile = find_gfile(savedir = savedir, file = "flywire_meta", folder = folder)
    gcsv = readr::read_csv(gfile)
    gcsv
  }
}

#' @rdname flywire_googledrive_data
#' @export
flywire_contributions <-function(local = FALSE, folder = "flywire_neurons/", sql = TRUE, ...){
  savedir = good_savedir(local = local)
  if(sql){
    find_gsql(savedir = savedir, tab = "flywire_edits", sql.db = "flywire_data.sqlite", folder = folder, ...)
  }else{
    gfile = find_gfile(savedir = savedir, file = "flywire_edits", folder = folder)
    gcsv = readr::read_csv(gfile)
    gcsv
  }
}

#' @rdname flywire_googledrive_data
#' @export
flywire_ids <-function(local = FALSE, folder = "flywire_neurons/", sql = TRUE, ...){
  savedir = good_savedir(local = local)
  if(sql){
    find_gsql(savedir = savedir, tab = "flywire_ids", sql.db = "flywire_data.sqlite", folder = folder, ...)
  }else{
    gfile = find_gfile(savedir = savedir, file = "flywire_ids", folder = folder)
    gcsv = readr::read_csv(gfile)
    as.character(gcsv$x)
  }
}

#' Update the flywire.id column in a set of google sheets based on flywire xyz positions
#'
#' @description This function retreives flywire IDs based on xyz positions in flywire voxel space, from a set of google sheets.
#' It also writes the updated flywire IDs to the same google sheets. This is often helpful because flywire IDs are inherently unstable, they change every time
#' a neuron is modified even slightly. Users can record 'stable' points in a neuron that identify it, e.g. a single xyz position in the cell body fibre, or at the soma, and
#' then use this function to update and get the correct flywire ID whenever they wish.
#'
#' @param selected_sheets character vector. the google sheet(s) to update. Each entry is a unique google sheet ID. You can find these in a sheet's URL.
#' If \code{NULL} then defaults to \code{option('hemibrainr_gsheets')}.
#' @param chosen.columns as well as writing column updates to the specified google sheets, this function returns a \code{data.frame} built fromm all given sheets and their
#' individual tabs, that have been updated. This argument specifies which column you want returned. Filled with NAs if it does not exist.
#' @param numCores if run in parallel, the number of cores to use. This is not necessary unless you have >10k points and want to see if you can get a speed up.
#' @param max.tries maximum number of attempts to write to/read from the google sheets before aborting. Sometimes attempts fail due to sporadic connections or API issues.
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
#' ## https://docs.google.com/spreadsheets/d/1rzG1MuZYacM-vbW7100aK8HeA-BY6dWAVXQ7TB6E2cQ/edit#gid=0
#' ## And return
#' fw.gsheet.meta = flywire_ids_update(selected_sheets = "1rzG1MuZYacM-vbW7100aK8HeA-BY6dWAVXQ7TB6E2cQ")
#'
#' }}
#' @seealso \code{\link{flywire_ids}},
#'   \code{\link{flywire_meta}},
#'   \code{\link{flywire_neurons}}
#' @name flywire_ids_update
#' @export
flywire_ids_update <- function(selected_sheets = NULL,
                               chosen.columns = c("fw.x","fw.y",'fw.z', 'flywire.xyz',
                                                  "flywire.id", "skid",
                                                  "FAFB.xyz", "cell.type", "side",
                                                  "ItoLee_Hemilineage", "Hartenstein_Hemilineage",
                                                  "hemibrain_match"),
                               numCores = 1,
                               max.tries = 10){
  if(is.null(selected_sheets)){
    selected_sheets = getOption("hemibrainr_gsheets", stop("Please set option('hemibrainr_gsheets')"))
  }
  fw.columns = c("flywire.id","fw.x","fw.y","fw.z","flywire.xyz")
  gs = data.frame(stringsAsFactors = FALSE)
  for(selected_sheet in selected_sheets){
    ## Read Google sheets and extract glywire neuron positions
    tabs = gsheet_manipulation(FUN = googlesheets4::sheet_names,
                               ss = selected_sheet,
                               return = TRUE)
    for(tab in tabs){
      gs.t = gs.t.current = gsheet_manipulation(FUN = googlesheets4::read_sheet,
                                                wait = 20,
                                                ss = selected_sheet,
                                                sheet = tab,
                                                guess_max = 3000,
                                                return = TRUE)
      used.cols = colnames(gs.t.current)
      if(nrow(gs.t)&&ncol(gs.t)&&sum(grepl("fw.x|flywire.xyz",used.cols))>0){
        # Separate x,y,z positions
        for(fw.c in fw.columns){
          if(is.null(gs.t[[fw.c]])){
            gs.t[[fw.c]] = NA
          }
        }
        good.xyz = sapply(gs.t$flywire.xyz,function(x) length(tryCatch(nat::xyzmatrix(x),error = function(e) NA))==3)
        if(sum(good.xyz)){
          gs.t[good.xyz,c("fw.x","fw.y","fw.z")] = nat::xyzmatrix(gs.t[good.xyz,"flywire.xyz"])
        }
        gs.t[!good.xyz,c("flywire.xyz")] = apply(gs.t[!good.xyz,c("fw.x","fw.y",'fw.z')],1,paste_coords)
        if(sum(is.na(gs.t$flywire.xyz))>0){
          # Get flywire IDs from these positions
          bbx = matrix(c(5100, 1440, 16, 59200, 29600, 7062),ncol=3,byrow = TRUE)
          bbx = nat::boundingbox(scale(bbx, scale = 1/c(4, 4, 40), center = FALSE))
          gs.t$fw.x = as.numeric(gs.t$fw.x)
          gs.t$fw.y = as.numeric(gs.t$fw.y)
          gs.t$fw.z = as.numeric(gs.t$fw.z)
          if(numCores>1){
            batch = 1
            batches = split(1:nrow(gs.t), round(seq(from = 1, to = numCores, length.out = nrow(gs.t))))
            foreach.ids <- foreach::foreach (batch = 1:length(batches)) %dopar% {
              pos = gs.t[batches[[batch]],]
              pos = pos[apply(pos, 1, function(row) sum(is.na(row[c("fw.x","fw.y",'fw.z')]))==0),]
              p = nat::pointsinside(pos[,c("fw.x","fw.y",'fw.z')],bbx)
              pos = pos[p,]
              pos = pos[!is.na(pos$fw.x),]
              if(nrow(pos)){
                tries = 0
                try.again = TRUE
                while(try.again){
                  tries = tries+1
                  i <- tryCatch(fafbseg::flywire_xyz2id(pos[,c("fw.x","fw.y",'fw.z')], rawcoords = TRUE),
                                error = function(e){cat(as.character(e));rep("failed",nrow(pos))})
                  names(i) = pos$flywire.xyz
                  i[is.na(i)|is.nan(i)] = "failed"
                  if(all(i%in%c("failed"))&tries<max.tries){
                    try.again = TRUE
                  }else{
                    try.again = FALSE
                  }
                }
                i
              }
            }
          }else{
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
          }
          fids = unlist(foreach.ids)
          fids[is.na(fids)|is.nan(fids)] = "0"
          gs.t$flywire.id = fids[match(gs.t$flywire.xyz ,names(fids))]
          if(nrow(gs.t)!=nrow(gs.t.current)){
            stop("Sheet processing corruption.")
          }
          write.cols = intersect(fw.columns,used.cols)
          if(length(fids) & length(write.cols)){
            for(column in write.cols){
              letter = LETTERS[match(column,colnames(gs.t))]
              range = paste0(letter,2,":",letter,nrow(gs.t)+1)
              gsheet_manipulation(FUN = googlesheets4::range_write,
                                  ss = selected_sheet,
                                  range = range,
                                  data = as.data.frame(gs.t[,column], stringsAsFactors = FALSE),
                                  sheet = tab,
                                  col_names = FALSE)
            }
          }
        }
        # Now continue processing
        gs.t = gs.t[,colnames(gs.t)%in%chosen.columns]
        for(col in chosen.columns){
          if(is.null(gs.t[[col]])){
            gs.t[[col]] = NA
          }
        }
        gs = plyr::rbind.fill(gs.t,gs)
      }
    }
  }
  # Make this unique, but keep row with most information
  master = gs
  master = master[!is.na(master$fw.x),]
  master = master[!is.na(master$fw.y),]
  master = master[!is.na(master$fw.z),]
  master$filled = apply(master, 1, function(r) sum(!is.na(r)))
  master = master[order(master$filled,decreasing = TRUE),]
  master = master[!duplicated(master$flywire.xyz),]
  rownames(master) = master$flywire.xyz
  master = subset(master, !is.na(master$flywire.id))
  master$filled = NULL
  master
}

#' Update neuron match information on a google sheet
#'
#' @description This function retreives neuron matches (\code{\link{hemibrain_matches}}) from a master-matching google sheet. If then
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
#' @param chosen.columns as well as writing column updates to the specified google sheets, this function returns a \code{data.frame} built fromm all given sheets and their
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
      gs.t = gs.t.current = gsheet_manipulation(FUN = googlesheets4::read_sheet,
                                                wait = 20,
                                                ss = selected_sheet,
                                                sheet = tab,
                                                guess_max = 3000,
                                                return = TRUE)
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
          for(column in write.cols){
            letter = LETTERS[match(column,colnames(gs.t))]
            range = paste0(letter,2,":",letter,nrow(gs.t)+1)
            gsheet_manipulation(FUN = googlesheets4::range_write,
                                ss = selected_sheet,
                                range = range,
                                data = as.data.frame(gs.t[,column], stringsAsFactors = FALSE),
                                sheet = tab,
                                col_names = FALSE)
          }
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


