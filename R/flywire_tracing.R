# lineage google sheets

#' Retrieve or go to a specific google sheet describing flywire neurons.
#'
#' @description These functions can help you
#' either go to or read the right tab (\code{flywire_tracing_sheet}) on a large google sheet, or read the whole sheet and all of its tabs ((\code{flywire_tracing_sheets}),
#' note, memoised for 30 minutes).
#' You can also see if a flywire neuron of interest in in a google sheet (\code{flywire_in(fw.meta = flywire_tracing_sheets())}),
#' and see if a flywire neuron of interest is among the data stored on the hemibrainr drive, i.e. (\code{flywire_in(fw.meta = flywire_meta())}).
#'
#' @param ws character vector, tabs (i.e. work sheets) on the google sheet to query/read. This works with regex, if \code{regex==TRUE}, so you only need to give the name partially.
#' If set to \code{NULL} for \code{flywire_tracing_sheets}, the whole google sheet is read and all tabs are combined using \code{plyr::rbind.fill}.
#' @param regex logical, use \code{ws} with regex.
#' @param open logial, if \code{TRUE} the relevant google sheet tab is opened in your browser. Else, a \code{data.frame} of the tab is returned.
#' @param selected_sheet character vector indicating the a flywire tracing google sheet. This defaults to the master 'lineage' tracing sheet used by the Drosphila Connectomics Group to
#' store annotations on flywire neurons in different developmental lineages. Note: You may not have access.
#' @param query a vector of data indicating flywire neurons. These can be a vector of flywire XYZ positions (readable by \code{nat::xyzmatrix}),
#' flywire root IDs or flywire supervoxel IDs. These are used to get the most up to date root IDs from FlyWire, which are then matched to the
#' results from \code{flywire_tracing_sheets}.
#' @param query.type  wtherh the qurey is a vector of xyz positions, flywire supervoxel IDs or flywire IDs ("flywire.xyz"/"flywire.svid"/"flywire.id").
#' @param fw.meta a \code{data.frame} of meta data on flywire neurons, e.g. as produced by \code{flywire_tracing_sheets}.or read by \code{flywire_meta}.
#' @param cloudvolume.url URL for CloudVolume to fetch segmentation image data. The default value of NULL choose
#' @param Verbose logical, if \code{TRUE} then \code{hemibrainr} communicates what it has found.
#' @param ... methods passed to \code{googlesheets4} functions.
#'.
#' @details flywire tracing google sheets should have the columns: flywire.xyz (a single cardinal position for each neuron in raw FlyWire voxel space),
#'flywire.svid (the supervoxel ID for this position) and flywire.id (the flywire root ID for the neuron that can be found at this position). Additional, informative columns
#'are encouraged. A sheet in this format can be updated using \code{\link{flywire_ids_update}}.
#'This updates the 'unstable' flywire.id column based on the stable flywire.xyz column, and also provides
#'stable flywire.svid column if it is empty.
#'
#' @return a \code{data.frame} with columns from the right tab, chosen using \code{regex}, from the given google sheet(s) \code{selected_sheet},
#' specified using the argument \code{chosen.columns}.
#'
#' @examples
#' \donttest{
#' \dontrun{
#'
#' # Quick look at what we have on the google sheet
#' sheetlist=googlesheets4::sheet_properties(
#' "1QyuHFdqz705OSxXNynD9moIsLvZGjjBjylx5sGZP2Yg")
#' View(sheetlist)
#'
#' # Go to specific tab of lineage google sheet
#' flywire_tracing_sheet("AOTUv1_v", open = TRUE)
#' flywire_tracing_sheet("CREa1_v.*_r", open = TRUE)
#'
#' # Get this tab as a data frame
#' flywire.lineage.sheet = flywire_tracing_sheet("CREa1_v.*_l")
#'
#' # Get the whole sheet as a data frame
#' flywire.lineage.sheet = flywire_tracing_sheets()
#' ## Note, the memoisation lasts for 30 minutes.
#' ## I.e. this takes a long time to run but is saved for
#' ## 30 mins on yuor computer. So if the function is called again
#' ## with the same arguments, the sheets are not read again.
#'
#' # See if we have a certain flywire neuron in our sheets
#' in.lin = flywire_in(query = "111437.5,21121.25,2661",
#' query.type = "flywire.xyz", Verbose = TRUE,
#' fw.meta = flywire_tracing_sheets())
#'
#' # And in the meta data for neurons on the google drive
#' in.meta = flywire_in(query = "111437.5,21121.25,2661",
#' query.type = "flywire.xyz", Verbose = TRUE)
#'
#' }}
#' @seealso \code{\link{flywire_ids}},
#'   \code{\link{flywire_meta}},
#'   \code{\link{flywire_neurons}}
#' @name flywire_tracing_sheet
#' @export
flywire_tracing_sheet <- function(ws,
                                  regex = FALSE,
                                  open=FALSE,
                                  selected_sheet = options()$flywire_lineages_gsheet,
                                  Verbose = TRUE,
                                  ...) {
  if(open&interactive()){
    u= sprintf("https://docs.google.com/spreadsheets/d/%s/edit", selected_sheet)
    sel=regex_tab_names(regex=ws,selected_sheet=selected_sheet,...)
    if(nrow(sel)>1) {
      warning("Multiple matches. Keeping first: ", paste(sel$name, collapse = ','))
      sel=sel[1,,drop=F]
    }
    uu=paste0(u, "#gid=", sel$id)
    utils::browseURL(uu)
  }else{
    if(regex){
      sel=regex_tab_names(regex=ws,selected_sheet=selected_sheet,...)
      if(nrow(sel)>1) {
        warning("Multiple matches. Keeping first: ", paste(sel$name, collapse = ','))
        sel=sel[1,,drop=F]
      }
      ws = sel$name
    }
    gs = data.frame(strstringsAsFactors = FALSE)
    for(w in ws){
      gm = gsheet_manipulation(FUN = googlesheets4::read_sheet,
                          wait = 20,
                          ss = selected_sheet,
                          guess_max = 3000,
                          sheet = w,
                          return = TRUE,
                          Verbose = Verbose,
                          ...)
      gs = plyr::rbind.fill(gs, gm)
    }
    gs
  }
}

# hidden
regex_tab_names <- function(regex, selected_sheet, ...){
  u=sprintf("https://docs.google.com/spreadsheets/d/%s/edit", selected_sheet)
  props=sheet_properties.memo(u)
  m=grepl(regex, props$name, ...)
  sel=props[m,]
  if(nrow(sel)==0){
    stop("No matches for ", regex)
  }
  sel
}

#' @export
#' @rdname flywire_tracing_sheet
flywire_tracing_sheets <- memoise::memoise(function(ws = NULL,
                                                    selected_sheet = options()$flywire_lineages_gsheet){
  flywire_tracing_sheets.now(ws = ws, selected_sheet=selected_sheet)
},  ~memoise::timeout(30*60))

# hidden
flywire_tracing_sheets.now <- function(ws = NULL,
                                      selected_sheet = options()$flywire_lineages_gsheet){
  gs.lineages  = data.frame(stringsAsFactors = FALSE)
  if(!is.null(ws)){
    sel=regex_tab_names(regex=ws,selected_sheet=selected_sheet)
    tabs=sel$name
  }else{
    tabs=gsheet_manipulation(FUN = googlesheets4::sheet_names,
                                                    ss = selected_sheet,
                                                    return = TRUE)
  }
  pb = progress::progress_bar$new(
    format = "  downloading :what [:bar] :percent eta: :eta",
    clear = FALSE, total = length(tabs))
  for(tab in tabs){
    pb$tick(tokens = list(what = tab))
    gs.lin = flywire_tracing_sheet(ws = tab, selected_sheet=selected_sheet, Verbose = FALSE)
    gs.lin$ws = tab
    gs.lineages = plyr::rbind.fill(gs.lineages, gs.lin)
  }
  gs.lineages
}

#' @export
#' @rdname flywire_tracing_sheet
flywire_in <- function(query,
                             query.type = c("flywire.id","flywire.xyz","flywire.svid"),
                             ws = NULL,
                             fw.meta = flywire_meta(),
                             cloudvolume.url = NULL,
                             Verbose = TRUE){
  # Match arguments
  query.type = match.arg(query.type)
  if(Verbose){
    message(sprintf("Looking for %s queries of type %s",length(query),query.type))
  }
  # Get SVIDS from sheets
  sheet.svids = if(is.null(fw.meta$flywire.svid)){
    fafbseg::flywire_xyz2id(fw.meta$flywire.xyz, rawcoords = TRUE, root=FALSE)
  }else{
    fw.meta$flywire.svid
  }
  # Get query as a SVID
  if(query.type=="flywire.xyz"){
    flywire.ids = fafbseg::flywire_xyz2id(nat::xyzmatrix(query), rawcoords = TRUE, cloudvolume.url=cloudvolume.url)
  }else if(query.type=="flywire.svid"){
    flywire.ids = fafbseg::flywire_rootid(query, cloudvolume.url=cloudvolume.url)
  }
  flywire.svids = lapply(flywire.ids, fafbseg::flywire_leaves, cloudvolume.url=cloudvolume.url)
  names(flywire.svids) = query
  df = data.frame(stringsAsFactors = FALSE)
  for(fi in query){
    svids = unlist(flywire.svids[fi])
    idx = which(sheet.svids%in%svids)
    if(!length(idx)||is.na(idx)){
      next
    }else{
      df = plyr::rbind.fill(df, cbind(query=query, fw.meta[idx,]))
    }
  }
  if(Verbose){
    message(sprintf("%s of %s queries found:",sum(query%in%df$query),length(query)))
    key.cols = c("query", "flywire.xyz", "flywire.svid","flywire.id", "status", "cell.type", "ItoLee_Hemilineage",
                 "side", "hemibrain.match", "hemibrain.match.quality", "tab", "ws")
    print(knitr::kable(df[,intersect(key.cols,colnames(df))]))
  }
  df
}

#' @param tab the tab on \code{selected_sheet} to update.
#' @param update a \code{data.frame} to use to update columns or add rows to the google sheet tab specified. This item must have at least
#' the same number of entries as the google sheet tab.
#' @param write.cols character vector, the columns from update that will be used to replace the column in the google sheet.
#' @param by character, the column by which to join the old google sheet tab and \code{update} using \code{dplyr::inner_join}. This column
#' must be in both the google sheet tab and \code{update}.
#' @param return logical, whether or not to return the updated googlesheet tab as a \code{data.frame}.
#' @name flywire_tracing_sheet
#' @export
flywire_tracing_update <- function(tab,
                                   selected_sheet,
                                   update,
                                   write.cols = colnames(update),
                                   by = "flywire.id",
                                   Verbose = TRUE,
                                   return = FALSE){
  # Read sheet
  gs = try(flywire_tracing_sheet(ws=tab,open=FALSE,selected_sheet=selected_sheet, Verbose=Verbose), silent = TRUE)

  # Add new tab if needed
  if(class(gs)=="try-error"){
    tabs = gsheet_manipulation(FUN = googlesheets4::sheet_names,
                                            ss = selected_sheet,
                                            return = TRUE)
    if(!tab%in%tabs){
      warning("Adding new tab ", tab)
      gs.added = gsheet_manipulation(FUN = googlesheets4::sheet_add,
                                     ss = selected_sheet,
                                     sheet = tab)
      gs.added = gsheet_manipulation(FUN = googlesheets4::sheet_write,
                                     data = update,
                                     ss = selected_sheet,
                                     sheet = tab)
    }
    return(invisible())
  }

  # Make sure write.cols in update
  update = unlist_df(update)
  missing.cols = setdiff(write.cols,colnames(update))
  if(length(missing.cols)){
    stop("write.cols missing from update")
  }
  missing.ids = setdiff(gs[[by]],update[[by]])
  if(length(missing.cols)){
    warning("entries from ", by," missing in update")
  }

  # Append missing rows
  if(nrow(update)>nrow(gs)){
    missing.rows = subset(update, !updated[[by]]%in%gs[[by]])
    if(Verbose){
      message("Adding ", nrow(update)," new rows")
    }
    missing.cols = setdiff(colnames(gs),colnames(update))
    for(mc in missing.cols){
      missing.rows[[mc]] = NA
    }
    missing.rows = missing.rows[,colnames(gs)]
    gsheet_manipulation(FUN = googlesheets4::sheet_append,
                        data = missing.rows,
                        ss = selected_sheet,
                        sheet = tab,
                        Verbose = Verbose)
    gs = try(flywire_tracing_sheet(ws=tab,open=FALSE,selected_sheet=selected_sheet, Verbose=Verbose), silent = TRUE)
  }

  # Merge gs and update
  gs[,write.cols] = update[match(gs[[by]],update[[by]]),write.cols]
  if(!nrow(update)){
    stop("nothing to update")
  }

  #  Update columns
  if(!identical(gs,update)){
    gsheet_update_cols(
      write.cols = write.cols,
      gs=gs,
      selected_sheet=selected_sheet,
      sheet = tab,
      Verbose = Verbose)
  }

  # Return a data.frame that represents new sheet
  if(return){
    try(flywire_tracing_sheet(ws=tab,open=FALSE,selected_sheet=selected_sheet), silent = FALSE)
  }else{
    invisible()
  }
}

# hidden, caches result for 5min in current session
sheet_properties.memo <- memoise::memoise(googlesheets4::sheet_properties, ~memoise::timeout(5*60))

