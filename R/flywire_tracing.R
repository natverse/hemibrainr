# lineage google sheets

#' @title Retrieve or go to a specific google sheet describing flywire neurons.
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
#' @param selected_sheet character indicating the a flywire tracing google sheet. This defaults to the master 'lineage' tracing sheet used by the Drosphila Connectomics Group to
#' store annotations on flywire neurons in different developmental lineages. Note: You may not have access.
#' @param selected_sheets a character vector of multiple google sheet IDs, as for \code{selected_sheet}.
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
                                  selected_sheet = options()$flywire_lineages_gsheets,
                                  Verbose = TRUE,
                                  ...) {
  data = data.frame()
  for(ss in selected_sheet){
    if(open&interactive()){
      u= sprintf("https://docs.google.com/spreadsheets/d/%s/edit", ss)
      sel=regex_tab_names(regex=ws,selected_sheet=ss,...)
      if(nrow(sel)>1) {
        warning("Multiple matches. Keeping first: ", paste(sel$name, collapse = ','))
        sel=sel[1,,drop=F]
      }
      uu=paste0(u, "#gid=", sel$id)
      utils::browseURL(uu)
    }else{
      if(regex){
        sel=regex_tab_names(regex=ws,selected_sheet=ss,...)
        if(nrow(sel)>1) {
          warning("Multiple matches. Keeping first: ", paste(sel$name, collapse = ','))
          sel=sel[1,,drop=F]
        }
        ws = sel$name
      }
      gs = data.frame(stringsAsFactors = FALSE)
      for(w in ws){
        gm = gsheet_manipulation(FUN = googlesheets4::read_sheet,
                                 wait = 20,
                                 ss = ss,
                                 guess_max = 3000,
                                 sheet = w,
                                 return = TRUE,
                                 Verbose = Verbose,
                                 ...)
        gs = plyr::rbind.fill(gs, gm)
      }
      data = plyr::rbind.fill(data, gs)
    }
  }
  data
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
                                                    selected_sheet = options()$flywire_lineages_gsheets){
  flywire_tracing_sheets.now(ws = ws, selected_sheet=selected_sheet)
},  ~memoise::timeout(30*60))

# hidden
flywire_tracing_sheets.now <- function(ws = NULL,
                                      selected_sheet = options()$flywire_lineages_gsheets){

  data = data.frame()
  for(ss in selected_sheet){
    if(!is.null(ws)){
      sel=regex_tab_names(regex=ws,selected_sheet=ss)
      tabs=sel$name
    }else{
      tabs=gsheet_manipulation(FUN = googlesheets4::sheet_names,
                               ss = ss,
                               return = TRUE)
    }
    pb = progress::progress_bar$new(
      format = "  downloading :what [:bar] :percent eta: :eta",
      clear = FALSE, total = length(tabs))
    tracing.list = list()
    for(tab in tabs){
      pb$tick(tokens = list(what = paste0(ss,": ", tab)))
      gs.lin = flywire_tracing_sheet(ws = tab, selected_sheet=ss, Verbose = FALSE)
      if(sum(nrow(gs.lin))){
        gs.lin$ws = tab
        tracing.list[[tab]] = gs.lin
      }
    }
    df = do.call(plyr::rbind.fill,tracing.list)
    data = plyr::rbind.fill(data,df)
  }
  data
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
                                            return = TRUE,
                               Verbose = Verbose)
    if(!tab%in%tabs){
      warning("Adding new tab ", tab)
      gs.added = gsheet_manipulation(FUN = googlesheets4::sheet_add,
                                     ss = selected_sheet,
                                     sheet = tab,
                                     Verbose = Verbose)
      gs.added = gsheet_manipulation(FUN = googlesheets4::sheet_write,
                                     data = update,
                                     ss = selected_sheet,
                                     sheet = tab,
                                     Verbose = Verbose)
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
    missing.rows = subset(update, !update[[by]]%in%gs[[by]])
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
    if(isTRUE(nrow(gs)>0)){
      gsheet_update_cols(
        write.cols = write.cols,
        gs=gs,
        selected_sheet=selected_sheet,
        sheet = tab,
        Verbose = Verbose)
    }
  }

  # Return a data.frame that represents new sheet
  if(return){
    try(flywire_tracing_sheet(ws=tab,open=FALSE,selected_sheet=selected_sheet, Verbose = Verbose), silent = FALSE)
  }else{
    invisible()
  }
}

# Standardise flywire sheet
#' @param whimsy logical. If TRUE then a column with randomly generated 'whimsical' names
#' is generated (unless there is already a name in this column).
#' @param field character, the the colum in the google sheet that you wish to have unique,
#' when \code{remove.duplicates} is \code{TRUE}.
#' @param remove.duplicates logical, whether or not to remvoe duplicate rows from google sheet. Currently,
#' this will also trigger \code{reorder} even if it is set to \code{FALSE}.
#' @param reorder logical, if \code{TRUE} then the googlesheet is reordered by cell type, connection weight and users.
#' @name flywire_tracing_sheet
#' @export
flywire_tracing_standardise <- function(ws = NULL,
                                        regex = FALSE,
                                        field = "flywire.id",
                                        selected_sheets = options()$flywire_lineages_gsheets,
                                        Verbose = TRUE,
                                        whimsy = FALSE,
                                        reorder = TRUE,
                                        remove.duplicates = TRUE){
  if(length(selected_sheets)>1){
    for(selected_sheet in selected_sheets){
      if(!is.null(ws)){ warning("examining all tabs on each sheet")}
      flywire_tracing_standardise(ws = NULL,
                                  field = field,
                                  regex=FALSE,
                                  whimsy=whimsy,
                                  selected_sheets=selected_sheet,
                                  Verbose=Verbose,
                                  reorder = reorder,
                                  remove.duplicates = remove.duplicates)
    }
    return(invisible())
  }
  if(is.null(ws)){
    if(Verbose){ message("### Working on google sheet: ", selected_sheets) }
    tabs = gsheet_manipulation(FUN = googlesheets4::sheet_names,
                                            ss = selected_sheets,
                                            return = TRUE,
                                            Verbose = FALSE)
    tabs = tabs[tabs!="readme"]
    pb = progress::progress_bar$new(
      format = "  standardising :what [:bar] :percent eta: :eta",
      clear = !Verbose, total = length(tabs))
    for(tab in tabs){
      try(pb$tick(tokens = list(what = tab)), silent = TRUE)
      flywire_tracing_standardise(ws=tab,
                                  field = field,
                                  regex=FALSE,
                                  whimsy=whimsy,
                                  selected_sheets=selected_sheets,
                                  Verbose=Verbose,
                                  reorder = reorder,
                                  remove.duplicates = remove.duplicates)
    }
    return(invisible())
  }else{
    tab = ws
  }
  gs = update = flywire_tracing_sheet(ws=tab,regex=regex,open=FALSE,selected_sheet=selected_sheets,Verbose=Verbose)
  write.cols = colnames(gs)
  if(!is.null(gs$status)){
    update$status = standard_statuses(update$status)
    write.cols = union("status",write.cols)
  }
  if(!is.null(gs$workflow)){
    update$workflow = standard_workflow(update$workflow)
    write.cols = union("workflow",write.cols)
  }
  if(whimsy){
    write.cols = union(write.cols,"whimsy")
    if(is.null(update$whimsy)){
      update$whimsy  = randomwords(n=nrow(update),words= 2,collapse = "_")
    }else{
      whismy.na = is.na(update$whimsy)
      if(sum(whismy.na)){
        update$whimsy[whismy.na] = randomwords(n=sum(whismy.na),words= 2,collapse = "_")
      }
    }
    update$whimsy = gsub(" ","_",update$whimsy)
  }
  bys = colnames(update)[grepl("\\_by$",colnames(update))]
  if(length(bys)){
    dupes = duplicated(update[[field]])
    for(by in bys){
      update[[by]] = gsub("DUPLICATED_","",update[[by]])
      update[[by]][dupes] = paste0("DUPLICATED_", update[[by]][dupes])
    }
  }
  if(isTRUE(nrow(update)>0)){
    if(reorder||remove.duplicates){
      if(!identical(gs,update)|reorder){
        gsheet_reorder(gs=update,
                       tab=tab,
                       selected_sheet=selected_sheets,
                       field = field,
                       remove.duplicates = remove.duplicates,
                       Verbose = Verbose)
      }
    }else{
      if(!identical(gs,update)){
        gsheet_update_cols(
          write.cols = write.cols,
          gs = update,
          selected_sheet=selected_sheets,
          sheet = tab,
          Verbose = Verbose)
      }
    }
  }
}

#' @param main_sheet a google sheet of flywire neurons with columns: \code{workflow}, \code{whimsy}, \code{flywire.id}. See \code{\link{standard_workflow}}.
#' This sheet is ead, and a separate tab in \code{target_sheet} is made for each workflow.
#' @param target_sheet a sheet to which to add 'workflow' tabs. Workflow tabs are lists of flywire neurons up/downstream of neurons
#' entered into the \code{main_sheet}. Accepted workflows: inputs, outputs, matches.
#' @param transmitters logical, if \code{TRUE} then inputs/outputs workflows include transmitter predictions for partners.
#' @param work.flows a workflow for a flywire neuron. Either we want to trace upstream of it, downstream of it, or match it to a hemibrain neuron.
#' @inheritParams fafbseg::flywire_ntpred
#' @inheritParams fafbseg::flywire_partner_summary
#' @name flywire_tracing_sheet
#' @export
flywire_deploy_workflows <-function(ws = "flywire",
                                      main_sheet, # "1nVEkC-WBcRMODhkKAp5KW5OWIRHzdxpY2ipFOV7r7k4"
                                      target_sheet = main_sheet, # "1WI7ri9yHkCGXDZ68PM5PnwAL6mw2keW7QYcxtz9Fwtw"
                                      regex = FALSE,
                                      threshold = 10,
                                      cleft.threshold = 100,
                                      transmitters = FALSE,
                                      local = NULL,
                                      cloudvolume.url = NULL,
                                      Verbose = TRUE,
                                      work.flows = c("inputs","outputs","matches","synapses")){
  if(length(setdiff(work.flows,c("inputs","outputs","matches","synapses")))){
    stop("The only supported workflows are: inputs, outputs and matches")
  }
  gs = flywire_tracing_sheet(ws=ws,regex=regex,open=FALSE,selected_sheet=main_sheet,Verbose=Verbose)
  if(!all(c("whimsy","workflow")%in%colnames(gs))){
    stop("Please give a column named 'whimsy' with human-memorable names for neurons\n Make sure there is a flywire.id column with valid entires.")
  }
  gs = gs[!is.na(gs$whimsy),]
  gs = gs[!duplicated(gs$whimsy),]
  gs$workflow = standard_workflow(gs$workflow)
  gs$status = standard_statuses(gs$status)
  to.work = subset(gs, grepl(paste(work.flows,collapse="|"),gs$workflow))
  if(!nrow(to.work)){
    message("No neurons in need of a tracing workflow")
    return(invisible())
  }
  tasks = sapply(strsplit(to.work$workflow,"/"), function(x) x[x%in%work.flows])
  indices = unlist(sapply(seq_along(tasks), function(x) rep(x,length((tasks[[x]])))))
  main = to.work[indices,]
  main$workflow = unlist(tasks)
  main$tab = workflow.tabs = paste(main$whimsy,main$workflow,sep="_")
  tabs = gsheet_manipulation(FUN = googlesheets4::sheet_names,
                                          ss = target_sheet,
                                          return = TRUE)
  workflow.tabs.missing = setdiff(workflow.tabs,tabs)
  workflow.tabs.there= intersect(workflow.tabs,tabs)
  pbm = progress::progress_bar$new(
    format = "  adding new tab :what [:bar] :percent eta: :eta",
    clear = FALSE, total = length(tabs))
  for(missing in workflow.tabs.missing){
    try(pbm$tick(tokens = list(what = missing)), silent = TRUE)
    please.add = subset(main, main$tab == missing)
    fw = flywire_workflow(flywire.id = please.add$flywire.id,
                          status = please.add$status,
                          ws=missing,
                          threshold = threshold,
                          cleft.threshold = cleft.threshold,
                          transmitters=transmitters,
                          local = local,
                          cloudvolume.url = cloudvolume.url,
                          Verbose = Verbose)
    gs.added = gsheet_manipulation(FUN = googlesheets4::sheet_add,
                                                ss = target_sheet,
                                                sheet = missing,
                                                Verbose = Verbose)
    gs.added = gsheet_manipulation(FUN = googlesheets4::sheet_write,
                                                data = fw,
                                                ss = target_sheet,
                                                sheet = missing,
                                                Verbose = Verbose)
  }
  if(length(workflow.tabs.there)){
    flywire_update_workflow(
      main = main,
      ws=workflow.tabs.there,
                            target_sheet=target_sheet,
                            Verbose = Verbose,
                            threshold = threshold,
                            cleft.threshold = cleft.threshold,
                            transmitters=transmitters,
                            local = local,
                            cloudvolume.url = cloudvolume.url)
  }
}

# hidden
flywire_update_workflow <-function(main,
                                   ws,
                                   target_sheet, # "1WI7ri9yHkCGXDZ68PM5PnwAL6mw2keW7QYcxtz9Fwtw"
                                   Verbose = TRUE,
                                   threshold = 10,
                                   cleft.threshold = 100,
                                   transmitters = FALSE,
                                   local = NULL,
                                   cloudvolume.url = NULL){
  if(length(ws)>1){
    pb = progress::progress_bar$new(
      format = "  updating tab :what [:bar] :percent eta: :eta",
      clear = FALSE, total = length(ws))
    for(item in ws){
      pb$tick(tokens = list(what = item))
      Sys.sleep(5) # to prevent server overload
      flywire_update_workflow(main = subset(main, main$tab == item),
                                         ws=item,
                                         target_sheet=target_sheet,
                                         Verbose = FALSE,
                                         threshold = threshold,
                                         cleft.threshold = cleft.threshold,
                                         transmitters=transmitters,
                                         local = local,
                                         cloudvolume.url = cloudvolume.url)
    }
  }else{
    gs = update = flywire_tracing_sheet(ws=ws,regex=FALSE,open=FALSE,selected_sheet=target_sheet,Verbose=Verbose)
    if(!is.null(gs$status)){
      gs$status = standard_statuses(gs$status)
      id = ifelse("pre_id"%in%colnames(gs),"pre_id","post_id")
      gs[[id]][is.na(gs[[id]])] = "0"
      gs[[id]] = tryCatch(fafbseg::flywire_latestid(gs[[id]]), error = function(e){
        warning(e)
        unlist(sapply(gs[[id]], function(x) tryCatch(fafbseg::flywire_latestid(x), error = function(e) NA)))
      })
      if(!nrow(gs)){
        warning("Workflow sheet has no rows: ", ws)
        return(invisible())
      }
      fw = flywire_workflow(flywire.id = main$flywire.id,
                            status = main$status,
                            ws=ws,
                            threshold = threshold,
                            cleft.threshold = cleft.threshold,
                            transmitters=transmitters,
                            local = local,
                            cloudvolume.url = cloudvolume.url,
                            Verbose = Verbose)
      shared.cols = setdiff(intersect(colnames(fw),colnames(gs)),id)
      for(sc in shared.cols){
        fw[[sc]] = gs[[sc]][match(fw[[id]],gs[[id]])]
      }
    }
    gsheet_manipulation(FUN = googlesheets4::sheet_write,
                                     data = fw,
                                     ss = target_sheet,
                                     sheet = ws,
                                     Verbose = Verbose)
  }
}

# hidden
#' @import fafbseg
flywire_workflow <- function(flywire.id,
                             status = "unassessed",
                              ws,
                              threshold = 10,
                              cleft.threshold = 100,
                             transmitters = FALSE,
                             nblast = NULL,
                              local = NULL,
                              cloudvolume.url = NULL,
                             Verbose = TRUE,
                             max.hits = 50,
                             db = flywire_neurons(WithConnectors = TRUE)){
  if(length(flywire.id)>1){
    stop("Only one flywire.id at a time please")
  }
  match = grepl("match",ws)
  syns = grepl("synapses",ws)
  if(match){
    if(is.null(nblast)){
      nblast = hemibrain_nblast("hemibrain-flywire")
    }
    entries = tab.entries[,colnames(nblast)%in%flywire.id]
    entries = sort(entries[entries>0],decreasing = TRUE)
    entries = entries[1:min(max.hits,length(entries))]
    tab.entries = as.data.frame(bodyid = names(entries),
                                nblast = unname(entries),
                                quality = "none",
                                note = NA,
                                stringsAsFactors = FALSE)
    main = data.frame(flywire.id, nblast = "main", quality = "none", note = gsub("_.*","",ws), stringsAsFactors = FALSE)
    colnames(main) = colnames(tab.entries)
    tab.entries = as.data.frame(tab.entries, stringsAsFactors = FALSE)
    plyr::rbind.fill(main, tab.entries)
  }else if (syns){
    tab.entries = flywire_synapse_annotations(flywire.id,
                                            partners = NULL,
                                            db = db,
                                            keep.dist.nm = NULL,
                                            cleft_scores.thresh = 50,
                                            sample = NULL,
                                            write.csv = FALSE)
    tab.entries
  }else{
    partners = ifelse(grepl("outputs",ws),"outputs","inputs")
    tab.entries = fafbseg::flywire_partner_summary(flywire.id,
                                                    partners = partners,
                                                    threshold = threshold,
                                                    remove_autapses = TRUE,
                                                    cleft.threshold = cleft.threshold,
                                                    details = FALSE,
                                                    roots = TRUE,
                                                    Verbose = Verbose,
                                                    local = local,
                                                    cloudvolume.url=cloudvolume.url)[,1:2]
    if(transmitters){
      ntpredictions=try(fafbseg:::ntpredictions_tbl(local=local),silent=TRUE)
      if(is.null(ntpredictions)){
        warning("Cannot find transmitter predictions")
      } else {
        nts.all = list()
        for(i in tab.entries[,1]){
          nt = fafbseg::flywire_ntpred(i,
                                       local = local,
                                       cleft.threshold=cleft.threshold,
                                       cloudvolume.url=cloudvolume.url)
          top.nt = names(sort(table(nt$top.nt),decreasing = TRUE))[1]
          nts.all[[i]] = data.frame(id = i, top.nt = top.nt)
        }
        nts = do.call(plyr::rbind.fill, nts.all)
        tab.entries$top.nt = nts$top.nt[match(tab.entries[,1],nts$id)]
      }
    }
    tab.entries$status = "unassessed"
    tab.entries$note = NA
    main = data.frame(flywire.id,weight = "main", status = status, note = gsub("_.*","",ws), stringsAsFactors = FALSE)
    colnames(main) = colnames(tab.entries)
    tab.entries = as.data.frame(tab.entries, stringsAsFactors = FALSE)
    plyr::rbind.fill(main, tab.entries)
  }
}

# hidden, caches result for  5min in current session
sheet_properties.memo <- memoise::memoise(googlesheets4::sheet_properties, ~memoise::timeout(5*60))

# Find descending neurons in flywire
#' @param side logical. The side of the brain these neurons are manually annotated, as being on.
#' @param gsheets character, google sheet codes for google sheets independently containing flywire IDs for right and left DNs.
#' @param chosen.columns character, which columns to read from the chosen google sheets..
#' @name flywire_tracing_sheet
#' @export
flywire_dns <- function(side = c("both","right","left"),
                        gsheets = c(left = "1Gq_-L1tpvuSxYs5O-_PggiI2qdF52xoTG64WzSs9KTU", right = "10T0JE6nVSz_uUdoHGOpV2odO_k75-arRPKdOiBXlS80"),
                        chosen.columns = c("flywire.id","flywire_id","flywire.xyz","Tags"),
                        ...){
  side = match.arg(side)
  if(side=="both"){
    side = c("right","left")
  }
  gs = data.frame()
  for(si in side){
    selected_sheet = gsheets[si]
    gm = gsheet_manipulation(FUN = googlesheets4::read_sheet,
                             wait = 20,
                             ss = selected_sheet,
                             guess_max = 3000,
                             sheet = "DNs",
                             return = TRUE,
                             Verbose = FALSE,
                             ...)
    gm = gm[,colnames(gm)%in%chosen.columns]
    gm$side = si
    gm$ws = "DNs"
    if(!is.null(gm$flywire.id)){
      gm$flywire.id = as.character(gm$flywire.id)
    }
    gs = plyr::rbind.fill(gs, gm)
  }
  gs
}

#' Generate a CSV of neuron synapses to import to a flywire annotation layer
#'
#' @param fw.ids character vector, a vector of valid flywire IDs
#' @param partners character vector else not used i \code{NULL}. A vector of
#'   valid flywire IDs for postsynaptic neurons to keep. Synapses to other
#'   postsynaptic targets are filtered out.
#' @param db a \code{neuronlist} of flywire neurons with synapses attached
#' @param keep.dist.nm numeric, minimum distance in nm that one synapse can be
#'   from another. One random synapse is chosen from groups of proximal
#'   synapses. To deactivate, enter \code{NULL}.
#' @param cleft_scores.thresh numeric, cleft_score threshold for synapse
#'   inclusion
#' @param sample numeric, the number of synapses to choose from both the
#'   dendrite and axon of each neuron in \code{fw.ids}
#' @param write.csv logical, whether or not to write a \code{.csv} output file,
#'   ready for import into flywire. One for each neuron, named by
#'   \code{csv.path}.
#' @param csv.path character, the path to which to save \code{.csv} files.
#' @export
flywire_synapse_annotations <- function(fw.ids,
                                        partners = NULL,
                                        db = flywire_neurons(WithConnectors = TRUE),
                                        keep.dist.nm = 500,
                                        cleft_scores.thresh = 75,
                                        sample = 250,
                                        write.csv = TRUE,
                                        csv.path = getwd()){
  flywire.scans = data.frame(stringsAsFactors = FALSE)
  for(fw.id in fw.ids){
    if(!fw.id%in%names(db)){
      warning(fw.id," missing from db")
      next
    }
    fw.neurons.syn.ac.syns = hemibrain_extract_synapses(db[fw.id], .parallel = TRUE, OmitFailures = TRUE)
    if(is.null(partners)){
      fw.neurons.syn.ac.syns <- fw.neurons.syn.ac.syns %>%
        dplyr::filter(partner %in% partners)
    }
    if(!nrow(fw.neurons.syn.ac.syns)){
      warning(fw.id," no valid synapses after partner filtering")
      next
    }
    if(is.null(fw.neurons.syn.ac.syns$Label)){
      fw.neurons.syn.ac.syns$Label = "unknown"
    }
    synister.synapse.sample <- fw.neurons.syn.ac.syns %>%
      dplyr::filter(prepost==0, cleft_scores > cleft_scores.thresh, Label %in% c("axon","dendrite", "2","3","unknown")) %>%
      dplyr::mutate(bin.cleft_scores= plyr::round_any(cleft_scores, 10)) %>%
      dplyr::group_by(Label) %>%
      dplyr::sample_n(size=sample, replace = TRUE) %>%
      dplyr::collect()
    if(!nrow(synister.synapse.sample)){
      warning(fw.id," no valid synapses after filtering")
      next
    }
    if(!is.null(keep.dist.nm)){
      d = t(apply(nat::xyzmatrix(synister.synapse.sample),1,function(r) r*c(4,4,40)))
      close = nabor::knn(d,d,k=10,radius=keep.dist.nm) # remove synapses within 500 nm of each other
      near = close$nn.idx
      keep = c()
      remove = c()
      for(r in 1:nrow(near)){
        remove = unique(c(remove,near[r,-1]))
        if(!r%in%remove){
          keep = c(keep, r)
        }
      }
      synister.synapse.sample = synister.synapse.sample[keep,]
    }
    synister.synapse.sample = as.data.frame(synister.synapse.sample, stringsAsFactors = FALSE)
    synister.synapse.sample = synister.synapse.sample[!duplicated(synister.synapse.sample$offset),]
    synister.synapse.sample$`Coordinate 1` = apply(nat::xyzmatrix(synister.synapse.sample),1,function(x) paste_coords(x/c(4,4,40)))
    flywire.scan = data.frame(`Coordinate 1` = synister.synapse.sample$`Coordinate 1`,
                              `Coordinate 2` = "",
                              `Ellipsoid Dimensions` = "",
                              tags = "",
                              Description = synister.synapse.sample$top.nt,
                              `Segment IDs` = "",
                              `Parent ID` = "",
                              Type = "Point",
                              ID = "",
                              offset = synister.synapse.sample$offset,
                              scores = synister.synapse.sample$scores,
                              cleft_scores = synister.synapse.sample$cleft_scores,
                              Label = synister.synapse.sample$Label,
                              flywire.id = fw.id)
    colnames(flywire.scan) = gsub("\\."," ",colnames(flywire.scan))
    flywire.scan$`Coordinate 1` = as.character(flywire.scan$`Coordinate 1`)
    if(write.csv){
      readr::write_excel_csv(flywire.scan, file = file.path(csv.path, paste0("flywire_",fw.id,"_synapse_annotations.csv")))
    }
    flywire.scan = plyr::rbind.fill(flywire.scans,flywire.scan)
  }
  flywire.scans
}
