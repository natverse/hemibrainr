#####################################################################################
################################ Save Manual Split ##################################
#####################################################################################

# # Settings
# bodyids = NULL
# phases = "I"
# db = NULL
# check_thresh = 1
# batch_size = 10
# by.type = TRUE
# assignments = TRUE
# brain = hemibrainr::hemibrain.surf
# update_regularly = TRUE
# motivate = TRUE
# clean = TRUE
# prioritise = TRUE
# selected_file = "1YjkVjokXL4p4Q6BR-rGGGKWecXU370D1YMc1mgUYr8E"
# selected_col = "#fadadd"

# hidden
## Set up Google Drive file as database
setup_splitcheck_sheet <-function(selected_file = "1YjkVjokXL4p4Q6BR-rGGGKWecXU370D1YMc1mgUYr8E"){
  ### Set up
  roots = subset(hemibrainr::hemibrain_all_splitpoints, hemibrainr::hemibrain_all_splitpoints$point == "root")
  lhn.gs = unlist_df(googlesheets4::read_sheet(ss = "1dH4d3-9aWLvoA8U3bz94JWMXOtW1krzLPnFc7tkgZnM"))
  rownames(lhn.gs) = lhn.gs$bodyId
  i = intersect(roots$bodyid,lhn.gs$bodyId)
  ### Underlying data.frame for database
  roots$type = hemibrainr::hemibrain_metrics[as.character(roots$bodyid),"type"]
  roots$type[is.na(roots$type)] = paste0("unknown_",roots$bodyid[is.na(roots$type)] )
  roots$soma.edit = !hemibrainr::hemibrain_metrics[as.character(roots$bodyid),"soma"]
  #roots[as.character(i),]$soma = unlist(nullToNA(lhn.gs[as.character(i),]$hasSoma))
  roots$cut = hemibrainr::hemibrain_metrics[as.character(roots$bodyid),"cropped"]
  roots[as.character(i),]$cut = unlist(nullToNA(lhn.gs[as.character(i),]$cropped2))
  roots$truncated = roots$cut
  roots$manual_edit = FALSE
  roots$splittable = TRUE
  roots$skeletonization = "good"
  roots$checked = 0
  roots$user = "flyconnectome"
  roots$time = Sys.time()
  roots$note = "none"
  ### Format for manual splitpoints sheet
  manual = roots
  ### Get ids for classes
  orn.ids = class2ids("ORN")
  hrn.ids = class2ids("HRN")
  pn.ids = class2ids("PN")
  upn.ids = class2ids("uPN")
  mpn.ids = class2ids("mPN")
  vppn.ids = class2ids("VPPN")
  alln.ids = class2ids("ALLN")
  dan.ids = class2ids("DAN")
  mbon.ids = class2ids("MBON")
  ### Prioritise the neurons we care about most
  roots$priority = 0
  roots$priority[roots$bodyid%in%hemibrainr::ton.ids] = 2
  roots$priority[roots$bodyid%in%hemibrainr::alln.ids] = 2
  roots$priority[roots$bodyid%in%hemibrainr::dan.ids] = 2
  roots$priority[roots$bodyid%in%hemibrainr::pn.ids] = 2
  roots$priority[roots$bodyid%in%hemibrainr::upn.ids] = 2
  roots$priority[roots$bodyid%in%hemibrainr::mpn.ids] = 2
  roots$priority[roots$bodyid%in%hemibrainr::mbon.ids] = 2
  roots$priority[roots$bodyid%in%hemibrainr::vppn.ids] = 1
  roots$priority[roots$bodyid%in%hemibrainr::orn.ids] = 1
  roots$priority[roots$bodyid%in%hemibrainr::hrn.ids] = 1
  roots$priority[roots$bodyid%in%lhn.gs$bodyId] = 3
  roots = roots[order(roots$bodyid, decreasing = TRUE),]
  roots = roots[order(roots$priority, decreasing = TRUE),]
  ### Get ready to record cable edits
  roots$edited.cable = 0
  ### Record original values
  roots$orig.soma = hemibrainr::hemibrain_metrics[as.character(roots$bodyid),"soma"]
  roots$orig.cut = hemibrainr::hemibrain_metrics[as.character(roots$bodyid),"cropped"]
  ### Assign users
  rownames(roots) = roots$bodyid
  ton.batches = split(hemibrainr::ton.ids, ceiling(seq_along(1:length(hemibrainr::ton.ids))/2800))
  roots[roots$bodyid%in%c(hemibrainr::upn.ids,hemibrainr::mpn.ids, hemibrainr::orn.ids),"user"] = "ND"
  roots[roots$bodyid%in%c(hemibrainr::dan.ids),"user"] = "GD"
  roots[roots$bodyid%in%c(hemibrainr::vppn.ids,hemibrainr::hrn.ids),"user"] = "RT"
  roots[roots$bodyid%in%c(hemibrainr::mbon.ids),"user"] = "MWP"
  roots[roots$bodyid%in%c(hemibrainr::alln.ids),"user"] = "TS"
  roots[roots$bodyid%in%c(hemibrainr::lhn.ids),"user"] = "AJ"
  roots[roots$bodyid%in%ton.batches[[1]],"user"] = "IT"
  roots[roots$bodyid%in%ton.batches[[2]],"user"] = "JH"
  roots = subset(roots, !is.na(roots$bodyid))
  ### Write to Google Sheet
  googlesheets4::write_sheet(roots[0,],
                             ss = selected_file,
                             sheet = "roots")
  batches = split(1:nrow(roots), ceiling(seq_along(1:nrow(roots))/5000))
  for(i in batches){
    gsheet_manipulation(FUN = googlesheets4::sheet_append,
                        data = roots[min(i):max(i),],
                        ss = selected_file,
                        sheet = "roots")
  }
  googlesheets4::write_sheet(manual[0,],
                             ss = selected_file,
                             sheet = "manual")
}

#' @param gs A Google sheet object, as read by \code{googlesheets4::read_sheet(selected_file)}.
#' @examples
#' \donttest{
#' \dontrun{
#' # Assign mPN review to Nik Drummond (ND)
#' hemibrain_task_update(bodyids = mpn.ids, column = "user", update = "ND")
#'
#' # Assign members of the FlyConnectome group to work on their favourite neurons:
#' ton.batches = split(hemibrainr::ton.ids,
#'   ceiling(seq_along(1:length(hemibrainr::ton.ids))/2700))
#' hemibrain_task_update(bodyids =
#'   c(hemibrainr::upn.ids, hemibrainr::mpn.ids, hemibrainr::orn.ids),
#'   column = "user", update = "ND")
#' hemibrain_task_update(bodyids = hemibrainr::dan.ids, column = "user", update = "GD")
#' hemibrain_task_update(bodyids = c(hemibrainr::vppn.ids,hemibrainr::hrn.ids),
#'   column = "user", update = "RT")
#' hemibrain_task_update(bodyids = hemibrainr::alln.ids, column = "user", update = "TS")
#' hemibrain_task_update(bodyids = hemibrainr::mbon.ids, column = "user", update = "MWP")
#' hemibrain_task_update(bodyids = ton.batches[[1]], column = "user", update = "AJ")
#' hemibrain_task_update(bodyids = ton.batches[[2]], column = "user", update = "IT")
#' hemibrain_task_update(bodyids = ton.batches[[3]], column = "user", update = "JH")
#' }}
#' @export
#' @rdname hemibrain_adjust_saved_split
hemibrain_task_update <- function(bodyids,
                                 update,
                                 gs = NULL,
                                 selected_file = "1YjkVjokXL4p4Q6BR-rGGGKWecXU370D1YMc1mgUYr8E",
                                 column = c("soma.edit","cut", "truncated", "manual_edit", "splittable", "checked", "user",
                                              "time", "note", "priority", "edited.cable", "skeletonization")
                                 ){
  column = match.arg(column)
  message("Updating task field: ", column)
  if(is.null(gs)){
    gs = googlesheets4::read_sheet(ss = selected_file, sheet = "roots")
    gs = unlist_df(gs)
  }
  bodyids = bodyids[bodyids%in%unlist(gs$bodyid)]
  rows = match(bodyids, gs$bodyid)
  rows = rows[!is.na(rows)]
  letter = LETTERS[match(column,colnames(gs))]
  if(length(update)!=length(rows)){
    update = rep(update[1], length(rows))
  }
  for(r in 1:length(rows)){
    row = rows[r]
    range = paste0(letter,row+1)
    gsheet_manipulation(FUN = googlesheets4::range_write,
                        ss = selected_file,
                        range = range,
                        data = data.frame(update[r]),
                        sheet = "roots",
                        col_names = FALSE)
  }
  message("Task updated! ")
}

#' Manually review split neurons and save your results to a Google Sheet for everyone to use
#'
#' @description A 3D interactive pipeline for reviewing our \code{\link{hemibrain_precomputed_splitpoints}}. Each user can
#' manually edit neurons if needs be, add whether or not they have a soma or are cropped, and make additional notes. The purpose of this
#' function if for interested users to work together to curate axon-dendrite splits suggested by the \code{\link{flow_centrality}} algorithm.
#' Currently, this function is built to review our most promising set of split points (\code{\link{hemibrain_splitpoints_polypre_centrifugal_synapses}}).
#' As the Google sheet gets updated, new manual splits will be re-saved into this package, as \code{hemibrain_manual_splitpoints} and implemented as default for those neurons
#' when using \code{\link{hemibrain_flow_centrality}} and \code{\link{hemibrain_read_neurons}}.For this function to work, you may have to add
#' the \href{https://docs.google.com/spreadsheets/d/1YjkVjokXL4p4Q6BR-rGGGKWecXU370D1YMc1mgUYr8E/edit?usp=sharing}{Google Sheet database} to your Google Drive.
#'
#' @inheritParams nlscan_split
#' @param brain which brain to plot while splitting neuron. If \code{NULL} no brain is plotted. Otherwise
#' it is recommended to use \code{hemibrain.surf} (which is in raw voxel space).
#' @param bodyids bodyids which specify individual neurons, as used in the relevant neuPrint project
#' @param phases which Phases of the split check process you want to engage with. Phase I involves choosing which neurons need to be manually split,
#' making notes on them and marking if they are truncated by by the volume or not Phase II involves manually splitting those neurons that have been flagged as
#' as in need to splitting. And Phase III involves reviewing the manual splits that have been made, and re-splitting neurons if needs be. The default, complete,
#' is to go through all three phases with sequential batches of neurons, of \code{batch_size}.
#' @param db a \code{neuronlist} object to act as a database of neuron skeletons. If you want to provide all
#' hemibrain neurons that could possibly be split, then you can use \code{\link{hemibrain_download_neurons}} to download
#' all hemibrain neurons as a \code{nat::neuronlistfh} file. All neurons can be read
#' using \code{db = hemibrain_read_neurons(savedir=TRUE)}. This can then be given as the argument to \code{db}.
#' @param check_thresh when checking skeleton splits and correcting them, we do not want to look at skeletons that have already been checked by others many times.
#' If a skeleton has already been checked by this many users or more, it will not be given to the user to edit.
#' @param batch_size when reviewing splits for hemibrain neurons, we examine neurons (or putative cell types, see parameter \code{by.type}) in batches of this size. We then go through two phases,
#' choosing neurons that are incorrect in some way (wrongly split, cropped) and then editing them manually, before saving them to a Google Sheet.
#' You may not have access to this Google Sheet, in which case you will not be able to contribute manual splits for neurons to this package. Contact us
#' if you would like access.
#' @param by.type logical. If \code{TRUE} then putative cell types (as assigned as the type field in neuPrint) of neurons are presented. The number of cell types pulled
#' are of size \code{batch_size}.
#' @param assignments logical. If TRUE then the user is invited to examine only neurons that they are flagged for them. These are entries into the 'roots' tab of the Google Sheet,
#' where the user column is the same as the initials the user enters.
#' @param update_regularly whether or not to read from the Google Sheet database regularly, to avoid duplicating manual splits. Should be set to \code{TRUE} unless
#' reading the Google Sheet introduces a prohibitive time delay.
#' @param motivate whether or not to plot 2D motivationals from \href{(https://inspirobot.me/)}{InspiroBot} to keep to entertained while splitting neurons.
#' @param prioritise whether or not to look at higher priority neurons first. Initial, this means olfactory system neurons.
#' @param clean whether or not to set synapse-less branches to \code{Label = 0}.
#' @param update a vector of \code{length(bodyids)} or \code{1}, which will replace cells in the column specified by \code{column} for given bodyids (\code{bodyids}).
#' E.g. If I wanted to changed the given bodyids on the task sheet to checked = TRUE, I would use \code{update = TRUE}.
#' @param column the column in the split-check task's Google Sheet (specified by \code{selected_file}) to be updated with the value of the argument \code{update}.
#' E.g. If I wanted to changed the given bodyids on the task sheet to checked = TRUE, I would use \code{column = "check"}.
#' @param ... additional arguments passed to plotting functions, such as \code{plot3d_split} and \code{nlscan_split}
#' @return Updates a Google Sheet which records whether neurons have been checked in one tab ("roots") and records actual manually made splits in a second tab ("manual").
#' @details the pipeline works in three phases and you choose to take a batch_size of X neurons/types (defaults to 10) through these stages:
#' \itemize{
#' \item{Phase I}{ - choose which neurons are okay and which need modifying /
#' are truncated / are unsplittable / have a bad soma, etc.}
#' \item{Phase II}{ - manually edit those neurons that need modifying. Top tip:
#' You can ‘cycle’ through branches and update them to a certain cable type one
#' at a time. This works by removing the primary neurite and linker,
#' so it is good to define these first for a neuron, and then cycle
#' through branches as a quick way to assign axon/dendrite.}
#' \item{Phase III}{ - Review manually modified neurons.}
#' }
#' If you have authentication issues with the Google sheet, look \href{(https://googlesheets4.tidyverse.org/articles/articles/auth.html)}{here} first.
#' Video tutorial: \url{https://drive.google.com/open?id=1FB7m1QFKApysO74maHAx-yg2CQWkx49C}
#' @examples
#' \donttest{
#' \dontrun{
#' # Plot splits to check them, correcting any errors, and saving your changes
#' nat::nopen3d()
#' hemibrain_adjust_saved_split()
#'
#' # In order to read neurons more quickly from a local source / Google Drive
#' db = hemibrain_neurons()
#' nat::nopen3d()
#' hemibrain_adjust_saved_split(db=db)
#'
#' # We can also try to correct somas en masse
#' nat::nopen3d()
#' hemibrain_adjust_saved_somas()
#' }}
#' @export
#' @rdname hemibrain_adjust_saved_split
hemibrain_adjust_saved_split <- function(bodyids = NULL,
                                         phases = c("complete", "I", "II", "III"),
                                         db = NULL,
                                         check_thresh = 1,
                                         batch_size = 10,
                                         by.type = TRUE,
                                         assignments = TRUE,
                                         brain = hemibrainr::hemibrain.surf,
                                         update_regularly = TRUE,
                                         motivate = TRUE,
                                         clean = TRUE,
                                         prioritise = TRUE,
                                         selected_file = "1YjkVjokXL4p4Q6BR-rGGGKWecXU370D1YMc1mgUYr8E",
                                         selected_col = "#fadadd",
                                         ...){
  ### Get Google sheet, assess what is to be done
  phases = match.arg(phases)
  if(motivate){plot_inspirobot()}
  initials = readline(prompt = "What are your initials? ")
  ### Get Google Sheet data
  gs = gsheet_manipulation(FUN = googlesheets4::read_sheet,
                      ss = selected_file,
                      sheet = "roots",
                      return = TRUE)
  if(is.null(gs)){
    stop("Google sheet database could not be read: ", selected_file)
  }
  manual = gsheet_manipulation(FUN = googlesheets4::read_sheet,
                           ss = selected_file,
                           sheet = "manual",
                           return = TRUE)
  manual = manual[!is.na(manual$bodyid),]
  say_hello(greet = initials)
  ### Process data
  message("
          ###Colours###
          By default:
          axon = orange,
          dendrite = blue,
          linker = green,
          primary neurite = purple,
          input = navy,
          output = red
          ###Colours###
          ")
  ### Get chosen bodyIDs
  undone.ids <- check_undoneids(bodyids = bodyids,
                                check_thresh = check_thresh,
                                gs = gs,
                                phases = phases,
                                assignments = assignments,
                                initials = initials,
                                manual = manual,
                                prioritise = prioritise,
                                by.type = by.type)
  ### Make batches
  if(by.type){
    message("You will be presented with ", batch_size, " cell types at a time (which may be more neurons)")
    batches = type_batches(undone.ids, gs = gs, batch_size = batch_size)
  }else{
    batches = split(undone.ids, ceiling(seq_along(undone.ids)/batch_size))
  }
  if(!length(batches)){
    stop("None of the given bodyIDs have been checked by less than ",check_thresh," user(s)")
  }else{
    for(batch in batches){
      exit = FALSE
      ### Get bodyids
      if(by.type){
        batch = gs$bodyid[gs$type%in%batch]
        batch = intersect(batch, undone.ids)
        if(is.issue(batch)){
          next
        }
        message("Neurons in batch: ", length(batch))
      }
      ### What priority are we at
      ranks = purify(gs[match(batch,gs$bodyid),"priority"])
      message("Highest priority neuron under consideration: ", max(ranks, na.rm=TRUE))
      ### Read batch
      message("Reading batch of ", length(batch)," hemibrain neurons")
      someneuronlist = pipeline_read_neurons(batch = batch, db = db, clean = clean, motivate = motivate)
      ### Get vectors that we will need to update
      edits = replace_with_none(purify(gs[match(names(someneuronlist), gs$bodyid),"manual_edit"]))
      notes = replace_with_none(purify(gs[match(names(someneuronlist), gs$bodyid),"note"]))
      truncated = replace_with_none(purify(gs[match(names(someneuronlist), gs$bodyid),"truncated"]))
      splittable = replace_with_none(purify(gs[match(names(someneuronlist), gs$bodyid),"splittable"]))
      cuts = replace_with_none(purify(gs[match(names(someneuronlist), gs$bodyid),"cut"]))
      skels = replace_with_none(purify(gs[match(names(someneuronlist), gs$bodyid),"skeletonization"]))
      soma.edits = replace_with_none(purify(gs[match(names(someneuronlist), gs$bodyid),"soma.edit"]))
      someneuronlist = hemibrain_settags(someneuronlist,
                                         manual_edit = edits,
                                         note = notes,
                                         truncated = truncated,
                                         splittable = splittable,
                                         cut = cuts,
                                         skeletonization = skels,
                                         soma.edit = soma.edits)
      message("Current entries:")
      print(knitr::kable(hemibrain_seetags(someneuronlist)))
      ### Let's Go. For. It.
      satisfied = FALSE
      while(!satisfied){
        ### Select neurons to edit and make notes on neurons
        if(phases %in% c("complete", "I")){
          message("Phase I : select neurons that need to be edited manually")
          phaseI <- splitcheck_phaseI(someneuronlist = someneuronlist,
                                      brain = brain,
                                      selected_col = selected_col)
          selected = phaseI[["selected"]]
          someneuronlist = phaseI[["someneuronlist"]]
          exit = phaseI[["exit"]]
          if(length(selected)){
            someneuronlist[as.character(selected)] <- hemibrain_settags(someneuronlist[as.character(selected)],manual_edit = rep(TRUE,length(selected)))
          }
          mes <- NULL
        }else{
          selected <- as.character(batch)
        }
        ### Manually fix selected skeletons
        if(phases %in% c("complete", "II")){
          message("Phase II : manually edit ", length(selected), " neurons")
          mes <- splitcheck_phaseII(someneuronlist = someneuronlist,
                                    selected = selected,
                                    brain = brain,
                                    manual = manual,
                                    update_regularly = update_regularly,
                                    selected_file = selected_file)
        }else if(phases == "I"){
          if(!is.issue(selected)){
            someneuronlist[selected] = hemibrain_settags(someneuronlist[selected],manual_edit = rep(TRUE,length(selected)))
          }
        }
        if(phases == "III" | length(mes)){
          ### Review the changes
          message("Phase III : review ", length(mes), " neurons")
          if(phases =="III"){
            message("Select neurons you may want to edit again")
          }
          mes <- splitcheck_phaseIII(mes = mes,
                                     selected = selected,
                                     manual = manual,
                                     motivate = motivate,
                                     brain = brain)
          phases = "complete"
        }
        ### Continue?
        update = prepare_update(someneuronlist=someneuronlist,gs=gs,initials=initials)
        satisfied = !hemibrain_choice(prompt = "Final check: Do you want to run the edit process again? yes/no
                                     ")
      }
      ### Save manual split to Google Sheet
      if(motivate){plot_inspirobot()}
      if(length(mes)){
        someneuronlist[names(mes)] = mes
        mes.sp = hemibrain_splitpoints(x = mes)
        mes.i = match(mes.sp$bodyid, names(someneuronlist))
        mes.sp$type = gs$type[match(mes.sp$bodyid,gs$bodyid)]
        mes.sp$checked = TRUE
        mes.sp$user = initials
        mes.sp$time = Sys.time()
        mes.tags = hemibrain_seetags(mes)
        mes.sp = merge(mes.sp,
                       mes.tags[,setdiff(c("bodyid",colnames(mes.tags)),colnames(mes.sp))],
                       all.x = TRUE, all.y = FALSE)
        if(nrow(mes.sp)){
          message("Saving new manual splits ...")
          success = FALSE
          gsheet_manipulation(FUN = googlesheets4::sheet_append,
                              data = mes.sp[,colnames(manual)],
                              ss = selected_file,
                              sheet = "manual")
          message("Saved new manual splits to Google sheet! ")
        }
      }
      ### Update our checks on Google Sheet
      message("Updating task completion ...")
      if(phases == "I" & update_regularly & sample(1:10,1)==10 ){
        message("Checking status of task by reading Google Sheet ...")
        message("Note: this can take ~ a minute")
        if(motivate){plot_inspirobot()}
        gs.safe = gs
        gs = gsheet_manipulation(FUN = googlesheets4::read_sheet,
                                 ss = selected_file,
                                 sheet = "roots",
                                 return = TRUE)
        if(is.null(gs)){
          gs = gs.safe
        }
      }
      update = prepare_update(someneuronlist=someneuronlist,gs=gs,initials=initials)
      rows = as.numeric(rownames(update))
      for(r in rows){
        range = paste0("H",r,":Q",r)
        up = update[as.character(r),intersect(colnames(gs),colnames(update))]
        if(sum(is.na(up))>1){
          message("Erroneuous NAs generated for row ", r, ", dropping this update")
          print(up)
          next
        }
        gsheet_manipulation(FUN = googlesheets4::range_write,
                            ss = selected_file,
                            range = range,
                            data = up,
                            sheet = "roots",
                            col_names = FALSE)
      }
      message("Task updated! ")
      say_encouragement(greet = initials)
      if(exit){
        break
      }
    }
  }
  message("That's all folks!")
}

# hidden
prepare_update <- function(someneuronlist,
                           gs,
                           initials = NULL,
                           print = TRUE){
  rows = match(names(someneuronlist), purify(gs$bodyid))
  update = hemibrain_seetags(someneuronlist)
  checked = suppress(as.numeric(purify(gs[rows,]$checked)))
  if(!is.null(initials)){
    users = unlist(gs[rows,]$user)
    users[users=="flyconnectome"] = ""
    users = gsub(initials,"",users)
    users= paste(users,initials, sep = "/")
    users = gsub("^/","",users)
    checked[is.na(checked)] = 0
    update2 = data.frame(
      checked = checked + 1,
      user = users,
      time = Sys.time())
    update = cbind(update,update2)
  }
  update[is.na(update)] = ""
  update = update3 = update[,colnames(update)!="bodyid"]
  rownames(update) = rows + 1
  if(print){
    message("Your updates: ")
    print(knitr::kable(update3))
  }
  update
}

# hidden
must_be <- function(prompt, answers){
  answer = "no answer yet fool"
  while(!answer %in% answers){
    answer = readline(prompt)
  }
  answer
}

# hidden
gsheet_manipulation <- function(FUN, ..., return = FALSE){
  sleep = 10
  success = FALSE
  while(!success){
    g = tryCatch(FUN(...),
                 error = function(e){
                   message(e)
                   return(NULL)
                 })
    if(!is.null(g)){
      success = TRUE
    }else{
      Sys.sleep(sleep)
      message("Google sheet read/write failure(s), re-trying in ", sleep," seconds ...")
      sleep = sleep + 10
      if(sleep > 600){
        slep <- 600
      }
    }
  }
  if(return){
    if(is.data.frame(g)){
      g = unlist_df(g)
    }
    return(g)
  }
}

# hidden
splitcheck_phaseI <- function(someneuronlist,
                              brain = hemibrainr::hemibrain.surf,
                              selected_col = "#fadadd"){
  message("
          ###################################PhaseI###################################
          Please select (s) neurons that need to be manually edited.
          Important: Please also mark neurons that have a bad/no soma (r),
          that are slightly truncated (t), seriously cropped (k), unsplittable (u),
          or have a bad skeletonization (z). You can also make custom notes (e).
          ###################################PhaseI####################################")
  neurons = names(someneuronlist)
  nams = as.data.frame(someneuronlist, stringsAsFactors = FALSE)$name
  frames <- length(neurons)
  selected <- character()
  i <- 1
  chc <- NULL
  rgl::bg3d(color = "white")
  reset3d(brain=brain)
  WithConnectors = TRUE
  exit = FALSE
  while(TRUE) {
    if (i > length(neurons) || i < 1){
      reset3d(brain=brain)
      if(length(selected)){
        rgl::plot3d(someneuronlist[selected], col = hemibrain_bright_colour_ramp(length(selected)))
      }
      end = hemibrain_choice("Done selecting neurons to edit (shown)? yes/no ")
      reset3d(brain=brain)
      if(end){
        break
      }else{
        i <- 1
      }
    }
    message("
########selection########
Please make notes on neurons
and select neurons for further
manual editing
########selection########")
    n <- neurons[i]
    cat("Current neuron:", n, " ", nams[i]," (", i, "/", length(neurons),")",ifelse(n %in% selected,"SELECTED","NOT SELECTED"), "\n")
    print(someneuronlist[n,], row.names = FALSE, right = TRUE)
    note = someneuronlist[[i]]$tags$note
    if(!is.issue(note)){
      message("Note for this neuron: ",note)
    }
    if(n %in% selected){
      rgl::bg3d(color = selected_col)
    }else{
      rgl::bg3d(color = "white")
    }
    plot3d_split(someneuronlist[i], soma.alpha = 0.5, WithConnectors = WithConnectors)
    chc <- must_be(prompt = "
Return to continue, b to go back, s to (de)select (needs edit), y to toggle synapses
n to make a note: truncated, cropped, bad soma, unsplittable, bad skeletoniztion, custom
c to cancel (with selection) and 'exit' to exit with current edits and end  ",
                   answers = c("","b","s","n","c","y","exit"))
    if (chc == "c") {
      if (is.null(chc) || chc == "c")
        chc = must_be(prompt = sprintf("Selection is: %s . Continue (c) or go back and make selection (b)? ",paste(ifelse(length(selected),selected,"none"),collapse=", ")), answers = c("c","b"))
      if (is.null(chc) || chc == "c"){
        break
      }
    }
    if (chc == "exit") {
      chc = must_be(prompt = paste0("Exit this pipeline and save current changes for the first", i ," neurons in this batch (exit)? Or go back and make selection (b)? "), answers = c("exit","b"))
      if (is.null(chc) || chc == "exit"){
        someneuronlist = someneuronlist[1:i]
        someneuronlist = nat::union(someneuronlist, someneuronlist[as.character(selected)])
        exit = TRUE
        break
      }
    }
    if (chc == "s") {
      if (n %in% selected) {
        message("Deselected: ", n)
        selected <- setdiff(selected, n)
      }
      else selected <- union(selected, n)
    }
    if (chc == "b"){
      i <- i - 1
    }else if (chc =="y"){
      WithConnectors = !WithConnectors
      message("Plot synapses: ", WithConnectors)
    }else if (chc == "n"){
      someneuronlist[[i]] = hemibrain_makenote(x=someneuronlist[[i]])
    }else{
      i <- i + 1
    }
    reset3d(brain=brain)
  }
  list(selected = selected,
       someneuronlist = someneuronlist,
       exit = exit)
}

# hidden
splitcheck_phaseII <- function(selected,
                               someneuronlist,
                               manual,
                               selected_file,
                               update_regularly = TRUE,
                               brain = hemibrainr::hemibrain.surf,
                               ...){
  message("
          #####################PhaseII#####################
          Please manually edit neurons that are not
          correctly split. Tip: mark out primary neurite
          and linker first, then cycle through branches (9)
          #####################PhaseII#####################")
  if(length(selected)){
    reset3d(brain=brain)
    rgl::bg3d(color = "white")
    if(update_regularly){
      manual.try = tryCatch(googlesheets4::read_sheet(ss = selected_file,
                                                      sheet = "manual"),
                            error = function(e) NULL)
      if(!is.null(manual.try)){
        manual = manual.try
        message("Manual splits by users: ")
        print(table(manual[!duplicated(manual$bodyid),]$user))
      }
    }
    if(any(selected%in%manual$bodyid)){
      reject = intersect(selected, as.character(manual$bodyid))
      message("The following bodyids have been manually split already. Please delete from
                Google Sheet (https://docs.google.com/spreadsheets/d/1YjkVjokXL4p4Q6BR-rGGGKWecXU370D1YMc1mgUYr8E/edit?usp=sharing) before adding them again: ", paste(reject,collapse=", "))
      utils::browseURL("https://docs.google.com/spreadsheets/d/1YjkVjokXL4p4Q6BR-rGGGKWecXU370D1YMc1mgUYr8E/edit?usp=sharing",
                       browser = getOption("browser"),
                       encodeIfNeeded = FALSE)
      selected = setdiff(selected, as.character(manual$bodyid))
    }
    me = nat::as.neuronlist(manually_assign_labels(someneuronlist[as.character(selected)]))
    edited = unlist(sapply(me, function(x) x$tags$manual_edit))
    mes = me[edited]
    message("Manually edited neurons: ", length(mes), " neurons")
    mes
  }else{
    message("No neurons selected for manual editing")
    NULL
  }
}

# hidden
splitcheck_phaseIII <- function(mes = NULL,
                                selected = NULL,
                                manual = NULL,
                                phases = c("complete", "I", "II", "III"),
                                motivate = TRUE,
                                brain = hemibrainr::hemibrain.surf,
                                ...){
  phases = match.arg(phases)
  message("
          ###################PhaseIII###################
          Please select neurons that are correctly split
          By default, all neurons are selected.
          ###################PhaseIII###################")
  if(is.null(mes)&phases=="III"){
    selected = intersect(selected, manual$bodyid)
    if(!is.issue(selected)){
      if(motivate){plot_inspirobot()}
      message("Reading ",length(selected)," neurons to check manual split")
      sel = hemibrain_read_neurons(as.character(selected), microns = FALSE, ...)
      mes = hemibrain_use_splitpoints(sel, df = manual)
    }
  }
  if(length(mes)){
    if(phases!="III"){
      save = hemibrain_choice(prompt = "Are you satisfied with your edits (y)? Or do you want to review them (n)? ")
      mes.sp = hemibrain_splitpoints(x = mes)
    }else{
      save = FALSE
      mes.sp = subset(manual, manual$bodyid %in% names(mes))
    }
    new.select = names(mes)
    while(!save){
      reset3d(brain=brain)
      message("Select (s) which neurons to save")
      mes2 = hemibrain_use_splitpoints(x = mes, df = mes.sp)
      message("Note: You must have a neuron selected (pink background) in order to save it.")
      new.select = nlscan_split(mes2, already_selected = names(mes2))
      message("Manually edited neurons: ", length(mes2), " neurons")
      if(!length(new.select)){
        new.select = "none"
      }
      message("Manual edits to be saved for: ",paste(new.select, collapse =", "))
      save = hemibrain_choice(prompt = "Are you satisfied with your edits (y)? Or do you want to review them again (n)? ")
    }
    mes = mes[names(mes)%in%new.select]
  }
  mes
}

# hidden
check_undoneids <- function(bodyids,
                            check_thresh = 1,
                            gs,
                            phases,
                            assignments,
                            initials,
                            manual,
                            prioritise,
                            by.type){
  ### Find undone
  lookedat = unlist(gs$checked)
  lookedat[lookedat=="FALSE"] = 0
  lookedat[lookedat=="TRUE"] = 1
  lookedat = as.numeric(lookedat)
  lookedat[is.na(lookedat)] = 0
  undone = gs[as.numeric(lookedat)<check_thresh,]
  if(phases=="II"){
    undone = gs[!gs$bodyid %in% unique(manual$bodyid) & gs$manual_edit,]
  }else if(phases=="III"){
    undone = gs[gs$bodyid %in% unique(manual$bodyid),]
  }
  ### Get unseen IDs
  undone.ids = unique(undone$bodyid)
  done = gs[as.numeric(lookedat)>=check_thresh,]
  message(length(unique(done$bodyid)), "/",length(unique(gs$bodyid))," have been checked by at least ",check_thresh," user")
  users = unlist(done[,"user"])
  print(table(users))
  message("Manual splits by users: ")
  print(table(manual[!duplicated(manual$bodyid),]$user))
  ### Choose particular IDs if selected
  if(!is.null(bodyids)){
    undone.ids = intersect(bodyids, undone.ids)
    nids = setdiff(bodyids, undone.ids)
    if(length(nids)){
      message("Some of the given IDs have already been examined by >=", check_thresh ," users: ", paste(nids,collapse=", "))
      message("Examining ", length(undone.ids)," ids")
    }
  }
  ### Get assigned IDs
  if(assignments){
    gs.b = gs[!is.na(gs$user),]
    gs.a = gs.b[gs.b$user==initials,]
    undone.ids.a = intersect(undone.ids,gs.a$bodyid)
    if(is.issue(undone.ids.a)){
      warning("There are no unchecked neurons to select assigned to ", initials)
    }else{
      message("Choosing only those neuron assigned to ", initials)
      undone.ids = undone.ids.a
    }
  }
  if(phases %in% c("II","III")){
    undone = subset(gs, gs$mes.tags %in% c("TRUE",TRUE))
    if(phases == "II"){
      undone.ids = setdiff(unique(undone$bodyid), manual$bodyid)
    }
    if(!length(undone.ids)){
      stop("There are no neurons marked for manual splitting (in Phase I) that have not had a manual split saved.")
    }
  }
  ### Prioritise
  if(prioritise){
    message("Prioritising 1st-3rd order olfactory ",ifelse(by.type,"neurons","putative cell types"))
    undone.ids = sample(undone.ids, length(undone.ids), replace = FALSE)
    ranks = as.numeric(gs[match(undone.ids,gs$bodyid),"priority"])
    undone.ids = undone.ids[order(ranks,decreasing = TRUE)]
    message("Priorities: ")
    print(table(ranks))
  }else{
    message("Choosing un-checked ",ifelse(by.type,"neurons","putative cell types")," at random")
    undone.ids = sample(undone.ids, length(undone.ids), replace = FALSE)
  }
  ### Do we have IDs?
  if(!length(undone.ids)){
    stop("There are no unchecked neurons to select for Phase:", phases)
  }
  undone.ids = as.character(undone.ids)
  undone.ids = undone.ids[!is.na(undone.ids)]
  undone.ids
}

# hidden
hemibrain_perfectstart <-function(x, ...) UseMethod("hemibrain_perfectstart")
hemibrain_perfectstart.neuron <- function(x, ...){
  x$tags$manual_edit = FALSE
  x$tags$cut = FALSE
  x$tags$soma.edit = "automatic"
  x$tags$truncated = FALSE
  x$tags$splittable = TRUE
  x$tags$skeletonization = "good"
  x$tags$soma.edit = "automatic"
  x$tags$note = "none"
  x
}
hemibrain_perfectstart.neuronlist<-function(x, ...){
  nat::nlapply(x,hemibrain_perfectstart.neuron, ...)
}

# hidden
hemibrain_makenote <- function(x){
  message("
########make notes########
Make notes on this neuron
########make notes########")
  while(TRUE){
    chc <- must_be(prompt = "t to mark as truncated, k to mark as seriously cropped, u to mark as unsplittable,
r to mark as missing soma, z to mark poor skeletonization, e to make a custom note
ENTER to continue (with notes made), c to cancel (without notes made).
                   ",
                   answers = c("t","c","e","u","r","k","z",""))
    x.safe = x
    if (chc == "e"){
      x$tags$note <- readline(prompt = "Add your note here: yes/no ")
      message("Note is: ", x$tags$note)
    }else if (chc == "t"){
      x$tags$truncated <- hemibrain_choice(prompt = "Is this neuron (even slightly) truncated?: yes/no ")
      message("Truncation status is: ",  x$tags$truncated)
    }else if (chc == "u"){
      x$tags$splittable <- hemibrain_choice(prompt = "Can this neuron be split? yes/no ")
      message("Splittable status is: ", x$tags$splittable)
    }else if (chc == "z"){
      dec<- must_be(prompt = "Is this neuron appropriately skeletonized? good(g)/minor errors(m)/major errors(e) ",
                          answers = c("g","m", "e"))
      dec = ifelse(dec == "g","good",dec)
      dec = ifelse(dec == "m","minor errors","major errors")
      x$tags$skeletonization <- dec
      message("Skeletonization status is: ", x$tags$skeletonization)
    }else if (chc == "r"){
      x$tags$soma.edit <- !hemibrain_choice(prompt = "Is this neuron missing its soma/obvious tract to soma? yes/no ")
      if(x$tags$soma.edit){
        reset <- hemibrain_choice(prompt = "Does the soma need to be manually re-set? yes/no ")
        x$tags$soma.edit <- ifelse(reset,"manual",TRUE)
      }
      message("Soma status is: ", x$tags$soma.edit)
    }else if (chc == "k"){
      x$tags$cut <- hemibrain_iscropped()
      if(x$tags$cut %in% c(TRUE,"TRUE","MB")){
        x$tags$truncated <- TRUE
      }
      if (x$tags$cut %in% c(TRUE,"TRUE")){
        x$tags$splittable <- FALSE
      }
      message("Cut status is: ", x$tags$cut)
      message("Truncation status is: ", x$tags$truncated)
      message("Splittable status is: ", x$tags$splittable)
    }
    if(chc == "c"){
      x = x.safe
      message("Notes not saved")
      break
    }else if (chc == ""){
      message("Notes saved")
      break
    }
  }
  x
}

# hidden
type_batches <- function(ids, gs, batch_size = 10){
  pcts = unique(purify(gs[match(ids, gs$bodyid),"type"]))
  pcts = pcts[order(pcts)]
  batches = split(pcts, ceiling(seq_along(pcts)/batch_size))
  batches
}
