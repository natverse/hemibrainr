#' Download a large set of well-traced skeletonised neurons from FlyWire
#'
#'@description Get a large number of skeletonised neurons from FlyWire. The neuronlist is saved as a
#'  \code{nat::neuronlistfh} object so certain neurons may be read from it
#'  without loading the entire, large \code{neuronlist} into memory. You will need access to the hemibrain Google Team Drive and
#'  have it mounted with Google filestream.
#'  If you want to flag flywire neurons that should be added to the Google drive, without doing this yourself, you can use
#'  \code{flywire_request}. The \code{flywire_basics} function will calculate some useful meta-data, namely
#' a point in the primary neurite tract that can be used as a stable reference for this neuron.
#'
#' @param brain the brainspace in which hemibrain neurons have been registered. Defaults to raw voxel space for the FlyWire project.
#' @param x flywire IDs for desired. If left as \code{NULL}, all flywire neurons that can be summoned from the hemibrainr Google drive are summoned.
#' @param WithConnectors logical, if \code{TRUE} the neuron fetched have predicted synapses located within each neuron at \code{neuron$connectors} and
#' have a predicted axon/dendrite split at \code{neuron$d$Label} for each point in its skeleton, created using \code{\link{flow_centrality}}.
#' @param local \code{FALSE} or path. By default (\code{FALSE}) data is read from \code{options()$Drive_hemibrain_data}), but the user can specify an alternative path.
#' @param mirror logical, whether or not to read neurons that have been mirrored (i.e. flipped to the 'other' brain hemisphere).
#' @param zip logical. If \code{TRUE} then \code{nat::neuronlistz} is used to resad neurons from a \code{.zip} archive. Otherwise, \code{nat::neuronlistfh} is used to read neurons from a \code{.rds} file with a linked 'data' folder.
#' The \code{.zip} method is preffered. Both methods load a dynamic neuronlist, which only loads skeleton data into memroy at the point where it is processed.
#' @param flywire.neurons a \code{neuronlist} of flywire neurons in FlyWire space. The \code{flywire_basics} function will calculate some useful meta-data, namely
#' a point in the primary neurite tract that can be used as a stable reference for this neuron.
#' @param request a neuronlist, matrix of x,y,z position or flywire ID to add to a
#' \href{https://docs.google.com/spreadsheets/d/1rzG1MuZYacM-vbW7100aK8HeA-BY6dWAVXQ7TB6E2cQ/edit#gid=0}{Google sheet} that records flywire positions
#' flagged to be processed into neuron skeletons that can be called by \code{flywire_neurons}.
#' @param selected_sheet the Google sheet onto which to add new flywire coordinate. I.e. \href{https://docs.google.com/spreadsheets/d/1rzG1MuZYacM-vbW7100aK8HeA-BY6dWAVXQ7TB6E2cQ/edit#gid=0}{Google sheet}.
#' @param sheet the tab onto which to add your requests.
#' @param gsheet logical, whether or not the request is are googlesheet keys. If they are, every item in their \code{flywire.xyz} columns is added to
#' the sheet specified by \code{selected_sheet}.
#' @param type When using neurons with \code{flywire_neurons} from the Google drive a neuronlistfh object of neuron skeletons (default),
#' the saved \code{.swc} files, a \code{nat::dotprops} objects or a \code{nat::dotprops object} cut to the hemibrain volume. Not all brainspaces supported.
#' @param User character, the user who is added new data.
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
#' # Now get some split neurons with synapses
#' fw.split = flywire_neurons(WithConnectors=TRUE)
#' nopen3d()
#' nlscan_split(fw.split[sample(length(fw.split),15)])
#'
#' # Get neuron dotprops object in hemibrain brainspace,
#' ## ready for NBLAST with hemibrain neurons
#' fw.dps.cut = flywire_neurons(type="cut_dotprops", brain = "JRCFIB2018F")
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
                            WithConnectors = FALSE,
                            local = FALSE,
                            brain = c("FlyWire", "JRCFIB2018Fraw","JRCFIB2018F","FAFB","FAFB14","JFRC2", "JFRC2013","JRC2018F","FCWB"),
                            mirror = FALSE,
                            type = c("neurons", "swc","dotprops","cut_dotprops"),
                            zip = FALSE,
                            ...){
  brain = match.arg(brain)
  type = match.arg(type)
  pattern = ifelse(zip,"zip","rds")
  if(type=="neurons"){
    type=""
  }
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

  # Get synapses from private drive
  if(type%in%c("dotprops","cut_dotprops")){
    if(brain!="JRCFIB2018F"){
      brain = "JRCFIB2018F"
      fhdir = file.path(neuronsdir,brain,"/")
      warning("When type = 'dotprops' only JRCFIB2018F brainspace supported from hemibrain google drive")
    }
  }else{
    if(WithConnectors){
      type = ""; local = FALSE; mirror = FALSE; brain = "FlyWire"
      if(brain!="FlyWire"){
        warning("When WithConnectors = TRUE only FlyWire brainspace supported from hemibrain google drive")
      }
    }
  }

  # Exists?
  if(!file.exists(fhdir)){
    stop("Cannot find file: ", fhdir)
  }

  # Read
  if(type=="swc"){
    if(brain!="FlyWire"){
      brain = "FlyWire"
      warning("Only native flywire neurons in FlyWire space supported when swc = TRUE")
    }
    if(zip){
      warning("zip set to FALSE, when swc is TRUE")
      zip = FALSE
    }
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
    filelist = list.files(path = fhdir, pattern = pattern, full.names = TRUE)
    if(type=="cut_dotprops"){
      filelist = filelist[grepl(type,filelist)]
    }else if (WithConnectors){
      filelist = filelist[grepl('flow',filelist)]
    }else{
      filelist = filelist[!grepl("cut_dotprops",filelist)]
      filelist = filelist[grepl(type,filelist)]
    }
    filelist = filelist[grepl("mirror",filelist)==mirror]
    filelist = sort(filelist,decreasing = TRUE)
    if(length(filelist)){
      fh.file = filelist[1]
      if(zip){
        neurons.flow.fh = nat::neuronlistz(fh.file)
      }else{
        neurons.flow.fh = nat::read.neuronlistfh(fh.file)
      }
      test = tryCatch(neurons.fh[[1]], error = function(e){
        warning(e)
        try(file.remove(paste0(attributes(neurons.fh)$db@datafile,"___LOCK")), silent = TRUE)
      })
      attr(neurons.fh,"df") = neurons.fh[,]
    }else{
      warning(pattern, " file for neuronlist not found at: ", fhdir)
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
                            User = "hemibrainr",
                            ...){
  if(!requireNamespace("fafbseg", quietly = TRUE)) {
    stop("Please install fafbseg using:\n", call. = FALSE,
         "remotes::install_github('natverse/fafbseg')")
  }
  # What kind of request is it?
  if(gsheet){
    xyz = c()
    type = "googlesheet"
    for(req in request){
      tabs = hemibrainr:::gsheet_manipulation(FUN = googlesheets4::sheet_names,
                                              ss = req,
                                              return = TRUE)
      for(tab in tabs){
        fts = try(flywire_tracing_sheet(ws = tab,
                                        selected_sheet=req,
                                        Verbose = FALSE), silent = TRUE)
        xyz = unique(c(xyz,fts$flywire.xyz))
      }
      if(!length(xyz)){
        stop("No flywire positions to add")
      }
    }
    xyz=xyz[!is.na(xyz)]
    xyz = nat::xyzmatrix(xyz)
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
    xyz = do.call(rbind, lapply(fb[,"flywire.xyz"], nat::xyzmatrix))
  }else if(type!="googlesheet"){
    xyz = as.data.frame(nat::xyzmatrix(request), stringsAsFactors =FALSE)
  }
  fw.xyz = apply(xyz,1,paste_coords)
  fw.xyz = fw.xyz[fw.xyz!="(NA,NA,NA)"]

  # Add to Google sheet
  gs = try(flywire_tracing_sheet(ws=sheet,open=FALSE,selected_sheet=selected_sheet,Verbose=FALSE), silent = FALSE)
  fw.xyz = setdiff(fw.xyz,gs$flywire.xyz)
  if(is.null(fw.xyz)||!length(fw.xyz)){
    warning("No new flywire positions to add")
    return(NULL)
  }
  if(class(gs)!="try-error"){
    fw.xyz = setdiff(fw.xyz,gs$flywire.xyz)
    update = data.frame(User = "flywire", flywire.xyz = fw.xyz)
    for(col in setdiff(colnames(gs),colnames(update))){
      if(col=="status"){
        update[[col]] = "unassessed"
      }else if(col=="workflow"){
        update[[col]] = "trace"
      }else if(col=="added_by"){
        update[[col]] = User
      }else if(col=="User"){
        update[[col]] = User
      }else if(col=="User"){
        update[[col]] = randomwords(n=nrow(update),words= 2,collapse = "_")
      }else{
        update[[col]] = NA
      }
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
  message(sheet, ": ", nrow(update), " FlyWire positions added")
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
#' @param simplified flywire synapses, where nearby Buhmann predicted synapses (often multiple per real synapses) are algorithmically combined into approximate single synapses
#' using \code{hemibrainr:::flywire_synapse_simplify}.
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
                            y = "n",
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
                            root = "i",
                            top.p = "n",
                            top.nt = "c",
                            scores = 'n',
                            cleft_scores = 'n',
                            top.nt = "c",
                            gaba  = "n",
                            acetylcholine  = "n",
                            glutamate  = "n",
                            octopamine  = "n",
                            serotonin  = "n",
                            dopamine  = "n",
                            prepost = "i")

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
#' @param Verbose logical, whether or not to supply you with messages.
#' @param retry integer, sometimes \code{fafbseg::flywire_xyz2id} fails due to a server timeout. This is the number of times to re-arrempt failed calls before accepting defeat.
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
                               priority = c("FAFB","hemibrain"),
                               Verbose = TRUE,
                               retry = 1){
  if(is.null(selected_sheets)){
    selected_sheets = getOption("hemibrainr_gsheets", stop("Please set option('hemibrainr_gsheets')"))
  }
  fw.columns = c("flywire.id","fw.x","fw.y","fw.z","flywire.xyz","flywire.svid")
  tracing.list = list()
  for(selected_sheet in selected_sheets){
    ## Read Google sheets and extract glywire neuron positions
    if(Verbose) { message("### Working on google sheet: ", selected_sheet) }
    if(is.null(ws)){
      tabs = gsheet_manipulation(FUN = googlesheets4::sheet_names,
                                 ss = selected_sheet,
                                 return = TRUE,
                                 Verbose = FALSE)
    }else{
      if(regex){
        ws=regex_tab_names(regex=ws,selected_sheet=selected_sheet)$name
      }
      tabs=ws
    }
    pb = progress::progress_bar$new(
      format = "  updating tab :what [:bar] :percent eta: :eta",
      clear = !Verbose, total = length(tabs))
    for(tab in tabs){
      # Update progress bar
      pb$tick(tokens = list(what = tab))
      # Get google sheet tab as data frame
      gs.t = gs.t.current = flywire_tracing_sheet(ws = tab, selected_sheet=selected_sheet, Verbose = FALSE)
      # Get original columns for later
      used.cols = colnames(gs.t.current)
      ## if valid data frame with flywire columnes
      if(nrow(gs.t)&&ncol(gs.t)&&sum(fw.columns%in%used.cols)>0){
        # Separate x,y,z positions
        for(fw.c in fw.columns){
          if(is.null(gs.t[[fw.c]])){
            gs.t[[fw.c]] = NA
          }
        }
        # Get seemingly good xyz flywire positions
        good.xyz = sapply(gs.t$flywire.xyz,function(x) length(tryCatch(nat::xyzmatrix(x),error = function(e) NA))==3)
        # If coordinates separately given, combine
        gs.t[!good.xyz,c("flywire.xyz")] = apply(gs.t[!good.xyz,c("fw.x","fw.y",'fw.z')],1,paste_coords)
        good.xyz = sapply(gs.t$flywire.xyz,function(x) length(tryCatch(nat::xyzmatrix(x),error = function(e) NA))==3)
        # Provide separate columns for x,y,z positions
        if(sum(good.xyz)){
          gs.t[good.xyz,c("fw.x","fw.y","fw.z")] = nat::xyzmatrix(gs.t[good.xyz,"flywire.xyz"])
        }
        # If flywire.svid missing, update
        bad.svids = (is.na(gs.t$flywire.svid)|gs.t$flywire.svid=="0")&!is.na(gs.t$flywire.xyz)
        if(sum(bad.svids)>0 & "flywire.svid"%in%colnames(gs.t)){
          gs.t$flywire.svid[bad.svids] = tryCatch(fafbseg::flywire_xyz2id(nat::xyzmatrix(gs.t$flywire.xyz[bad.svids]),
                                                                 root=FALSE,
                                                                 rawcoords = TRUE),
                                                  error = function(e){
                                                    NA
                                                    warning(e)
                                                  })
        }
        fids = NULL
        if(!all(is.na(gs.t$flywire.xyz))){
          # Get the ids that need to be updated
          fwids.old = gs.t$flywire.id
          fwids.old.not.good = is.na(gs.t$flywire.id)|gs.t$flywire.id%in%c("0","NA",""," ","\n")|!grepl("^[0-9]{1,}$", gs.t$flywire.id)
          fwids.old[fwids.old.not.good] = 0
          latest = fafbseg::flywire_islatest(fwids.old)
          fwids.need.a.look = (fwids.old.not.good+!latest)>0
          gs.n = gs.t[fwids.need.a.look,]
          if(sum(fwids.need.a.look)>0){
            # Get flywire IDs from these positions
            bbx = matrix(c(5100, 1440, 16, 59200, 29600, 7062),ncol=3,byrow = TRUE)
            bbx = nat::boundingbox(scale(bbx, scale = 1/c(4, 4, 40), center = FALSE))
            gs.n$fw.x = as.numeric(gs.n$fw.x)
            gs.n$fw.y = as.numeric(gs.n$fw.y)
            gs.n$fw.z = as.numeric(gs.n$fw.z)
            pos = gs.n[apply(gs.n, 1, function(row) sum(is.na(row[c("fw.x","fw.y",'fw.z')]))==0),]
            p = nat::pointsinside(pos[,c("fw.x","fw.y",'fw.z')],bbx)
            pos = pos[p,]
            pos = pos[!is.na(pos$fw.x),]
            if(nrow(pos)){
              foreach.ids = try(fafbseg::flywire_xyz2id(pos[,c("fw.x","fw.y",'fw.z')], rawcoords = TRUE), silent = TRUE)
              sleep = 10
              if(class(foreach.ids)=="try-error" & retry>0){
                retry = as.integer(retry)
                while(retry>0&class(foreach.ids)=="try-error"){
                  Sys.sleep(sleep)
                  retry=retry-1
                  foreach.ids = try(fafbseg::flywire_xyz2id(pos[,c("fw.x","fw.y",'fw.z')], rawcoords = TRUE), silent = TRUE)
                  sleep = sleep+sleep
                }
              }
              if(class(foreach.ids)=="try-error"){
                stop(paste0("fafbseg::flywire_xyz2id could not be used for ", tab," in sheet ", selected_sheet))
              }
              names(foreach.ids) = pos[,"flywire.xyz"]
            }else{
              foreach.ids = NULL
            }
            fids = unlist(foreach.ids)
            fids[is.na(fids)|is.nan(fids)] = "0"
            replacement = fids[match(gs.n$flywire.xyz ,names(fids))]
            coordsmissing = gs.n$flywire.xyz==paste_coords(matrix(NA,ncol=3))
            coordsmissing[is.na(gs.n$flywire.xyz)] = TRUE
            coordsmissing[is.na(replacement)] = TRUE
            replacement[coordsmissing] = gs.n$flywire.id[coordsmissing]
            gs.t$flywire.id[fwids.need.a.look] = replacement
          }
        }
        # If ID is given and no xyz
        bad.xyz = (is.na(gs.t$flywire.xyz)|gs.t$flywire.xyz=="0")&!is.na(gs.t$flywire.id)
        if(sum(bad.xyz)>0){
          ids.fresh = try(fafbseg::flywire_latestid(gs.t$flywire.id[bad.xyz]),silent = TRUE)
          if(!"try-error"%in%class(ids.fresh)){
            gs.t$flywire.id[bad.xyz] = ids.fresh
          }else{
            warning(ids.fresh)
          }
        }
        # If only skid given try and guess flywire.id
        if("skid"%in%used.cols){
          justskids=((gs.t$flywire.xyz==paste_coords(matrix(NA,ncol=3))|is.na(gs.t$flywire.xyz))
                     &(is.na(gs.t$flywire.id)|gs.t$flywire.id==0)
                     &!is.na(gs.t$skid))
          justskids[is.na(justskids)] = FALSE
          if(sum(justskids)>0){
            if(Verbose) message("Geting flywire IDs for skids")
            replacement.ids = unlist(pbapply::pbsapply(gs.t[justskids,"skid"], function(x)
              tryCatch(suppress(fafb14_to_flywire_ids_timed(x, only.biggest = TRUE))$flywire.id,error=function(e){message(e);NA})))
            if(inherits(replacement.ids,"try-error")){
              warning(replacement.ids)
            }else{
              gs.t[justskids,"flywire.id"] = replacement.ids
            }
          }
        }
        # update with flywire.xyz if possible
        if(!is.null(meta)){
          justids=(gs.t$flywire.xyz==paste_coords(matrix(NA,ncol=3))
                  |is.na(gs.t$flywire.xyz)
                  &!is.na(gs.t$flywire.id))
          justids[is.na(justids)] = FALSE
          if(sum(justids)>0){
            replacement.xyz = meta[match(gs.t[justids,"flywire.id"],meta$flywire.id),]$flywire.xyz
            gs.t[justids,"flywire.xyz"] = replacement.xyz
          }
        }
        # Change 0 to NA
        good.xyz = sapply(gs.t$flywire.xyz,function(x) length(tryCatch(nat::xyzmatrix(x),error = function(e) NA))==3)
        gs.t$flywire.xyz[good.xyz] = try(apply(nat::xyzmatrix(gs.t$flywire.xyz[good.xyz]),1,paste_coords), silent = TRUE)
        gs.t$flywire.xyz[gs.t$flywire.xyz==paste_coords(matrix(NA,ncol=3))] = NA
        gs.t$flywire.id[gs.t$flywire.id==0]=NA
        # Write to google sheet
        if(nrow(gs.t)!=nrow(gs.t.current)){
          stop("Sheet processing corruption. Flywire.id update failed.")
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
        # Now continue processing
        gs.t = gs.t[,colnames(gs.t)%in%chosen.columns]
        for(col in chosen.columns){
          if(is.null(gs.t[[col]])){
            gs.t[[col]] = NA
          }
        }
        if(nrow(gs.t)){
          gs.t$gsheet = selected_sheet
          entry = paste0(selected_sheet,"_",tab)
          tracing.list[[entry]] = gs.t
        }
      }
    }
  }
  # Make this unique, but keep row with most information
  master = do.call(plyr::rbind.fill, tracing.list)
  master = master[!is.na(master$flywire.xyz),]
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
  tracing.list = list()
  for(selected_sheet in selected_sheets){
    ## Read Google sheets and extract glywire neuron positions
    tabs = gsheet_manipulation(FUN = googlesheets4::sheet_names,
                               ss = selected_sheet,
                               return = TRUE,
                               Verbose = FALSE)
    pb = progress::progress_bar$new(
      format = "  updating :what [:bar] :percent eta: :eta",
      clear = FALSE, total = length(tabs))
    for(tab in tabs){
      pb$tick(tokens = list(what = tab))
      gs.t = gs.t.current = flywire_tracing_sheet(ws = tab, selected_sheet=selected_sheet, Verbose = FALSE)
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
          suppressMessages(gsheet_update_cols(
            write.cols = write.cols,
            gs=gs.t[,used.cols],
            selected_sheet = selected_sheet,
            sheet = tab,
            Verbose = FALSE))
        }
        gs.t = gs.t[,colnames(gs.t)%in%chosen.columns]
        for(col in chosen.columns){
          if(is.null(gs.t[[col]])){
            gs.t[[col]] = NA
          }
        }
        entry = paste(selected_sheet,tab,collapse ="_")
        tracing.list[[entry]] = gs.t
      }
    }
  }
  do.call(plyr::rbind.fill, tracing.list)
}


