#####################################################################################
################################ Save Manual Split ##################################
#####################################################################################

# hidden
## Set up Google Drive file as database
setup_splitcheck_sheet <-function(ss = "1YjkVjokXL4p4Q6BR-rGGGKWecXU370D1YMc1mgUYr8E"){
  ### Set up
  roots = subset(hemibrainr::hemibrain_all_splitpoints, point == "root")
  lhn.gs = as.data.frame(googlesheets4::read_sheet(ss = "1dH4d3-9aWLvoA8U3bz94JWMXOtW1krzLPnFc7tkgZnM"))
  rownames(lhn.gs) = lhn.gs$bodyId
  i = intersect(roots$bodyid,lhn.gs$bodyId)
  ### Underlying data.frame for database
  roots$type = hemibrainr::hemibrain_metrics[as.character(roots$bodyid),"type"]
  roots$type[is.na(roots$type)] = paste0("unknown_",roots$bodyid[is.na(roots$type)] )
  roots$soma = hemibrainr::hemibrain_metrics[as.character(roots$bodyid),"soma"]
  roots[as.character(i),]$soma = unlist(nullToNA(lhn.gs[as.character(i),]$hasSoma))
  roots$cut = hemibrainr::hemibrain_metrics[as.character(roots$bodyid),"cropped"]
  roots[as.character(i),]$cut = unlist(nullToNA(lhn.gs[as.character(i),]$cropped2))
  roots$truncated = roots$cut
  roots$split = "automatic_synapses"
  roots$splittable = TRUE
  roots$checked = FALSE
  roots$user = "flyconnectome"
  roots$time = Sys.time()
  roots$note = ""
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
  roots$priority[roots$bodyid%in%lhn.gs$bodyId] = 3
  roots$priority[roots$bodyid%in%alln.ids] = 2
  roots$priority[roots$bodyid%in%dan.ids] = 2
  roots$priority[roots$bodyid%in%pn.ids] = 2
  roots$priority[roots$bodyid%in%upn.ids] = 2
  roots$priority[roots$bodyid%in%mpn.ids] = 2
  roots$priority[roots$bodyid%in%mbon.ids] = 2
  roots$priority[roots$bodyid%in%vppn.ids] = 1
  roots$priority[roots$bodyid%in%orn.ids] = 1
  roots$priority[roots$bodyid%in%hrn.ids] = 1
  roots = roots[order(roots$bodyid, decreasing = TRUE),]
  roots = roots[order(roots$priority, decreasing = TRUE),]
  ### Get ready to record cable edits
  roots$edited.cable = 0
  ### Record original values
  roots$orig.soma = hemibrainr::hemibrain_metrics[as.character(roots$bodyid),"soma"]
  roots$orig.cut = hemibrainr::hemibrain_metrics[as.character(roots$bodyid),"cropped"]
  ### Write to Google Sheet
  googlesheets4::write_sheet(roots[0,],
                             ss = ss,
                             sheet = "roots")
  batches = split(1:nrow(roots), ceiling(seq_along(1:nrow(roots))/1000))
  for(i in batches){
    gsheet_manipulation(FUN = googlesheets4::sheets_append,
                        data = roots[min(i):max(i),],
                        ss = selected_file,
                        sheet = "roots")
  }
  googlesheets4::write_sheet(manual[0,],
                             ss = ss,
                             sheet = "manual")
  ### Some assignments for tracers in the FlyConnectome group
  hemibrain_task_update(bodyids = c(upn.ids,mpn.ids), column = "user", update = "ND")
  hemibrain_task_update(bodyids = dan.ids, column = "user", update = "GD")
  hemibrain_task_update(bodyids = c(vppn.ids,hrn.ids), column = "user", update = "ECM")
  hemibrain_task_update(bodyids = alln.ids, column = "user", update = "TS")
  hemibrain_task_update(bodyids = mbon.ids, column = "user", update = "MWP")
  hemibrain_task_update(bodyids = ton.ids, column = "user", update = "AJ")
  hemibrain_task_update(bodyids = orn.ids, column = "user", update = "MMC")
}

#' @examples
#' \donttest{
#' \dontrun{
#' # Assign mPN review to Nik Drummond (ND)
#' hemibrain_task_update(bodyids = mpn.ids, column = "user", update = "ND")
#' }}
#' @export
#' @rdname hemibrain_adjust_saved_split
hemibrain_task_update <- function(bodyids,
                                  update,
                                 selected_file = "1YjkVjokXL4p4Q6BR-rGGGKWecXU370D1YMc1mgUYr8E",
                                 column = c("soma","cut", "truncated", "split", "splittable", "checked", "user",
                                              "time", "note", "priority", "edited.cable")
                                 ){
  column = match.arg(column)
  message("Updating task field: ", column)
  gs = googlesheets4::read_sheet(ss = selected_file, sheet = "roots")
  rows = match(bodyids, gs$bodyid)
  rows = rows[!is.na(rows)]
  letter = LETTERS[match(column,colnames(gs))]
  range = paste0(letter,min(rows)+1,":",letter,max(rows)+1)
  if(length(update)!=length(rows)){
    update = rep(update, length(rows))
  }
  gsheet_manipulation(FUN = googlesheets4::sheets_edit,
                        ss = selected_file,
                        range = range,
                        data = data.frame(update),
                        sheet = "roots",
                        col_names = FALSE)
  message("Task updated! ")
}

#' Manually review split neurons and save your results to a Google Sheet for everyone to use
#'
#' @description A 3D interactive pipeline for reviewing our \code{\link{hemibrain_precomputed_splitpoints}}. Each user can
#' manually edit neurons if needs be, add whether or not they have a soma or are cropped, and make additional notes. The purpose of this
#' function if for interested users to work together to curate axon-dendrite splits suggested by the \code{\link{flow_centrality}} algorithm.
#' Currently, this function is built to review our most promising set of split points (\code{\link{hemibrain_splitpoints_polypre_centrifugal_synapses}}).
#' As the google sheet gets updated, new manual splits will be re-saved into this package, as \code{hemibrain_manual_splitpoints} and implemented as default for those neurons
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
#'
#' \item Phase I - choose which neurons are okay and which need modifying /
#' are truncated / are unsplittable / have a bad soma, etc.}
#'
#' \item Phase II - manually edit those neurons that need modifying. Top tip:
#' You can ‘cycle’ through branches and update them to a certain cable type one
#' at a time. This works by removing the primary neurite and linker,
#' so it is good to define these first for a neuron, and then cycle
#' through branches as a quick way to assign axon/dendrite
#'
#' \item Phase III - Review manually modified neurons.
#' }
#'
#' Video tutorial: https://drive.google.com/open?id=1Cw-5bVgY5k-eaZAXQl6jgU1leYY7YbTM
#'
#' @examples
#' \donttest{
#' \dontrun{
#' # Plot splits to check them, correcting any errors, and saving your changes
#' nat::nopen3d()
#' hemibrain_adjust_saved_split()
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
  ### Get google sheet, assess what is to be done
  phases = match.arg(phases)
  if(motivate){plot_inspirobot()}
  initials = readline(prompt = "What are your initials? ")
  say_hello(greet = initials)
  ### Get Google Sheet data
  gs = googlesheets4::read_sheet(ss = selected_file, sheet = "roots")
  manual = googlesheets4::read_sheet(ss = selected_file, sheet = "manual")
  ### Process data
  undone = subset(gs, checked < check_thresh)
  undone.ids = unique(undone$bodyid)
  done = subset(gs, checked >= check_thresh)
  message(length(unique(done$bodyid)), "/",length(unique(gs$bodyid))," have been checked by at least ",check_thresh," user")
  users = gs[gs$checked,"user"]
  print(table(users))
  message("Manual splits by users: ")
  print(table(manual[!duplicated(manual$bodyid),]$user))
  message("
          ###Colours###
          By default: axon = orange, dendrite = blue,
          linker = green, primary neurite = purple,
          input = navy, output = red
          ###Colours###")
  ### Get chosen bodyIDs
  undone.ids <- check_undoneids(undone.ids = undone.ids,
                                bodyids = bodyids,
                                gs = gs,
                                phases = phases,
                                assignments = assignments,
                                initials = initials,
                                manual = manual)
  ### Prioritise
  if(prioritise){
    message("Prioritising 1st-3rd order olfactory ",ifelse(by.type,"neurons","putative cell types"))
    undone.ids = sample(undone.ids, length(undone.ids), replace = FALSE)
    ranks = gs[match(undone.ids,gs$bodyid),"priority"][[1]]
    undone.ids = undone.ids[order(ranks,decreasing = TRUE)]
    message("Priorities: ")
    print(table(ranks))
  }else{
    message("Choosing un-checked ",ifelse(by.type,"neurons","putative cell types")," at random")
    undone.ids = sample(undone.ids, length(undone.ids), replace = FALSE)
  }
  ### Make batches
  if(by.type){
    pcts = purify(gs[match(undone.ids, gs$bodyid),"type"])
    pcts = pcts[order(pcts)]
    batches = split(pcts, ceiling(seq_along(pcts)/batch_size))
  }else{
    batches = split(undone.ids, ceiling(seq_along(undone.ids)/batch_size))
  }
  if(!length(batches)){
    stop("None of the given bodyIDs have been checked by less than ",check_thresh," user(s)")
  }else{
    for(batch in batches){
      ### Get bodyids
      if(by.type){
        batch = purify(gs$bodyid[match(batch, gs$type)])
        message("Neurons in batch: ", length(batch))
      }
      ### What priority are we at
      ranks = gs[match(batch,gs$bodyid),"priority"][[1]]
      message("Highest priority neuron under consideration: ", max(ranks, na.rm=TRUE))
      ### Read batch
      message("Reading batch of ", batch_size," neurons from the hemibrain project")
      if(is.null(db)){
        if(motivate){plot_inspirobot()}
        message("Reading and manipulating neurons from neuPrint ...")
        someneuronlist = hemibrain_read_neurons(x = as.character(batch), microns = FALSE)
      }else{
        message("Reading locally saved neurons ...")
        someneuronlist = db[as.character(batch)]
      }
      if(clean){
        someneuronlist = hemibrain_clean_skeleton(someneuronlist, rval = "neuron")
      }
      ### Get vectors that we will need to update
      notes = as.character(gs[match(names(someneuronlist), gs$bodyid),"note"][[1]])
      truncated = as.character(gs[match(names(someneuronlist), gs$bodyid),"truncated"][[1]])
      somas = as.character(gs[match(names(someneuronlist), gs$bodyid),"soma"][[1]])
      splittable = as.character(gs[match(names(someneuronlist), gs$bodyid),"splittable"][[1]])
      ### Let's Go. For. It.
      satisfied = FALSE
      while(!satisfied){
        ### Select neurons to edit and make notes on neurons
        if(phases %in% c("complete", "I")){
          message("Phase I : select neurons that need to be edited manually")
          phaseI <- splitcheck_phaseI(someneuronlist = someneuronlist,
                                      notes = notes,
                                      truncated=truncated,
                                      somas=somas,
                                      brain = brain,
                                      splittable=splittable,
                                      selected_col = selected_col)
          selected = phaseI["selected"][[1]]
          notes = phaseI["notes"][[1]]
          truncated = phaseI["truncated"][[1]]
          somas = phaseI["somas"][[1]]
          splittable = phaseI["splittable"][[1]]
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
        }
        if(phases %in% c("complete", "III")){
          ### Review the changes
          message("Phase III : review ", length(mes), " neurons")
          mes <- splitcheck_phaseIII(mes = mes,
                                     selected = selected,
                                     manual = manual,
                                     motivate = motivate,
                                     brain = brain)
        }
        ### Continue?
        satisfied = hemibrain_choice(prompt = "Continue to saving changes (y), or would you like to run the edit process
        for neurons in this batch again (keeping curent changes) (n)? ")
        someneuronlist[names(mes)] = mes
      }
      ### Save manual split to Google Sheet
      if(motivate){plot_inspirobot()}
      if(length(mes)){
        mes.sp = hemibrain_splitpoints(x = mes)
        mes.i = match(mes.sp$bodyid, names(someneuronlist))
        mes.sp$type = gs$type[match(mes.sp$bodyid,gs$bodyid)]
        mes.sp$split = "manual"
        mes.sp$checked = TRUE
        mes.sp$user = initials
        mes.sp$cut = unlist(sapply(mes[as.character(mes.sp$bodyid)], function(x) nullToNA((x$tags$cropped))))
        mes.sp$truncated = truncated[mes.i]
        mes.sp$soma = unlist(sapply(mes[as.character(mes.sp$bodyid)], function(x) nullToNA((x$tags$soma))))
        mes.sp$splittable = splittable[mes.i]
        mes.sp$time = Sys.time()
        mes.sp$note = notes[mes.i]
        if(nrow(mes.sp)){
          message("Saving new manual splits ...")
          success = FALSE
          gsheet_manipulation(FUN = googlesheets4::sheets_append,
                              data = mes.sp[,colnames(manual)],
                              ss = selected_file,
                              sheet = "manual")
          message("Saved new manual splits to Google sheet! ")
        }
      }
      ### Update our checks on Google Sheet
      message("Updating task completion ...")
      if(phases == "I" & update_regularly){
        message("Checking status of task by reading Google Sheet ...")
        if(motivate){plot_inspirobot()}
        gs = googlesheets4::read_sheet(ss = selected_file, sheet = "roots")
      }
      rows = match(names(someneuronlist), gs$bodyid)
      splits = as.character(gs$split[rows])
      cuts = as.character(gs$cut[rows])
      if(length(mes)){
        which.mes = which(names(someneuronlist)%in%selected)
        splits[which.mes] = "manual"
        cuts[which.mes] = unlist(sapply(mes, function(x) x$tags$cropped))
        somas[which.mes] = unlist(sapply(mes, function(x) x$tags$soma))
      }
      somas[somas=="automatic"] = TRUE
      which.cuts = which(cuts%in%c(TRUE,"MB","TRUE"))
      truncated[which.cuts] = TRUE
      update = data.frame(
        soma = somas,
        cut = cuts,
        truncated = truncated,
        split = splits,
        splittable = splittable,
        checked = gs[rows,]$checked+1,
        user = rep(initials, length(rows)),
        time = Sys.time(),
        note = notes)
      update[is.na(update)] = ""
      rownames(update) = rows
      for(r in rows){
        range = paste0("H",r+1,":P",r+1)
        gsheet_manipulation(FUN = googlesheets4::sheets_edit,
                            ss = selected_file,
                            range = range,
                            data = update[as.character(r),intersect(colnames(gs),colnames(update))],
                            sheet = "roots",
                            col_names = FALSE)
      }
      message("Task updated! ")
      say_encouragement(greet = initials)
    }
  }
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
gsheet_manipulation <- function(FUN, ...){
  success = FALSE
  while(!success){
    g = tryCatch(FUN(...),
                 error = function(e) NULL)
    if(!is.null(g)){
      success = TRUE
    }else{
      Sys.sleep(5)
      message("Google sheet read failures, re-trying ...")
    }
  }
}

# hidden
splitcheck_phaseI <- function(someneuronlist,
                              notes,
                              truncated,
                              splittable,
                              somas,
                              brain = NULL,
                              selected_col = "#fadadd"
){
  if(nat::is.neuronlist(someneuronlist)) {
    neurons = names(someneuronlist)
    nams = as.data.frame(someneuronlist)$name
  }
  frames <- length(neurons)
  selected <- character()
  i <- 1
  chc <- NULL
  rgl::bg3d(color = "white")
  reset3d(brain=brain)
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
    n <- neurons[i]
    cat("Current neuron:", n, " ", nams[i]," (", i, "/", length(neurons),")",ifelse(n %in% selected,"SELECTED","NOT SELECTED"), "\n")
    print(someneuronlist[n,], row.names = FALSE, right = TRUE)
    note = notes[i]
    if(!is.issue(note)){
      message("Note for this neuron: ",note)
    }
    if(n %in% selected){
      rgl::bg3d(color = selected_col)
    }else{
      rgl::bg3d(color = "white")
    }
    plot3d_split(someneuronlist[i])
    chc <- must_be(prompt = "
          Return to continue, b to go back, s to select (needs edit / seriously cropped),
          t to mark as truncated, u to mark as unsplittable, r to mark as missing soma,
          c to cancel (with selection), e to make a note ",
                   answers = c("","b","s","t","c","e","u"))
    if (chc == "c") {
      if (is.null(chc) || chc == "c")
        chc = must_be(prompt = sprintf("Selection is: %s . Continue (c) or go back and make selection (b)? ",paste(ifelse(length(selected),selected,"none"),collapse=", ")), answers = c("c","b"))
      if (is.null(chc) || chc == "c"){
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
    }else if (chc == "e"){
      notes[i] <- readline(prompt = "Add your note here: yes/no ")
      message("Note is: ", notes[i])
    }else if (chc == "t"){
      truncated[i] <- readline(prompt = "Is this neuron (even slightly) truncated?: yes/no ")
      message("Truncation status is: ", truncated[i])
    }else if (chc == "u"){
      splittable[i] <- hemibrain_choice(prompt = "Can this neuron be split? yes/no ")
      message("Splittable status is: ", splittable[i])
    }else if (chc == "r"){
      somas[i] <- hemibrain_choice(prompt = "Is this neuron missing its soma/obvious tract to soma? yes/no ")
      message("Soma status is: ", splittable[i])
    }else{
      i <- i + 1
    }
    reset3d(brain=brain)
  }
  list(selected = selected,
       notes = notes,
       truncated = truncated,
       splittable = splittable,
       somas = somas
  )
}

# hidden
splitcheck_phaseII <- function(selected,
                               someneuronlist,
                               manual,
                               selected_file,
                               update_regularly = TRUE,
                               brain = NULL,
                               ...){
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
    me = manually_assign_labels(someneuronlist[as.character(selected)])
    edited = unlist(sapply(me, function(x) x$tags$manually_edited))
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
                                brain = NULL,
                                ...){
  phases = match.arg(phases)
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
      mes.sp = subset(manual, bodyid %in% names(mes))
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
      message("Manual edits to be saved for: ",paste(new.select, sep =", "))
      save = hemibrain_choice(prompt = "Are you satisfied with your edits (y)? Or do you want to review them again (n)? ")
    }
    mes = mes[names(mes)%in%new.select]
  }
  mes
}

# hidden
check_undoneids <- function(undone.ids, bodyids, gs, phases, assignments, initials, manual){
  if(!is.null(bodyids)){
    undone.ids = intersect(bodyids, undone.ids)
  }
  if(assignments){
    gs.a = gs[gs$user==initials,]
    undone.ids.a = intersect(undone.ids,gs.a$bodyid)
    if(is.issue(undone.ids.a)){
      warning("There are no unchecked neurons to select assigned to ", initials)
    }else{
      message("Choosing only those neuron assigned to ", initials)
      undone.ids = undone.ids.a
    }
  }
  if(phases %in% c("II","III")){
    undone = subset(gs, split == "manual")
    if(phases == "II"){
      undone.ids = setdiff(unique(undone$bodyid), manual$bodyid)
    }
    if(!length(undone.ids)){
      stop("There are no neurons marked for manual splitting (in Phase I) that have not had a manual split saved.")
    }
  }
  if(!length(undone.ids)){
    stop("There are no unchecked neurons to select for Phase:", phases)
  }
  as.character(undone.ids)
}
