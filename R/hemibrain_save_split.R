#####################################################################################
################################ Save Manual Split ##################################
#####################################################################################

# hidden
## Set up Google Drive file as database
setup_splitcheck_sheet <-function(){
  roots = subset(hemibrain_all_splitpoints, point == "root")
  lhn.gs = as.data.frame(googlesheets4::read_sheet(ss = "1dH4d3-9aWLvoA8U3bz94JWMXOtW1krzLPnFc7tkgZnM"))
  rownames(lhn.gs) = lhn.gs$bodyId
  i = intersect(roots$bodyid,lhn.gs$bodyId)
  roots$soma = hemibrain_metrics[as.character(roots$bodyid),"soma"]
  roots[as.character(i),]$soma = unlist(nullToNA(lhn.gs[as.character(i),]$hasSoma))
  roots$cut = hemibrain_metrics[as.character(roots$bodyid),"cropped"]
  roots[as.character(i),]$cut = unlist(nullToNA(lhn.gs[as.character(i),]$cropped2))
  roots$split = "automatic_synapses"
  roots$checked = FALSE
  roots$user = "flyconnectome"
  roots$time = Sys.time()
  roots$note = ""
  googlesheets4::write_sheet(roots[,],
                             ss = "1YjkVjokXL4p4Q6BR-rGGGKWecXU370D1YMc1mgUYr8E",
                             sheet = "roots")
  googlesheets4::write_sheet(roots[0,],
                             ss = "1YjkVjokXL4p4Q6BR-rGGGKWecXU370D1YMc1mgUYr8E",
                             sheet = "manual")
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
#' @param db a \code{neuronlist} object to act as a database of neuron skeletons. If you want to provide all
#' hemibrain neurons that could possibly be split, then you can use \code{\link{hemibrain_download_neurons}} to download
#' all hemibrain neurons as a \code{nat::neuronlistfh} file. All neurons can be read
#' using \code{db = hemibrain_read_neurons(savedir=TRUE)}. This can then be given as the argument to \code{db}.
#' @param check_thresh when checking skeleton splits and correcting them, we do not want to look at skeletons that have already been checked by others many times.
#' If a skeleton has already been checked by this many users or more, it will not be given to the user to edit.
#' @param batch_size when reviewing splits for hemibrain neurons, we examine neurons in batchs of this size. We then go through two phases,
#' choosing neurons that are incorrect in some way (wrongly split, cropped) and then editing them manually, before savign them to a Google Sheet.
#' You may not have access to this Google Sheet, in which case you will not be able to contribute manual splits for neurons to this package. Contact us
#' if you would like access.
#' @param update_regularly whether or not to read from the Google Sheet database regulalry, to avoid duplicating manual splits. Should be set to \code{TRUE} unless
#' reading the Google Sheet introduces a prohibitive time delay.
#' @param motivate whether or not to plot 2D motivationals from \href{(https://inspirobot.me/)}{InspiroBot} to keep to entertained while splitting neurons.
#' @param clean whether or not to set synapse-less branches to \code{Label = 0}.
#' @param ... additional arguments passed to plotting functions, such as \code{plot3d_split} and \code{nlscan_split}
#' @return The neuron/neuronlist object with axon/dendrite info assigned in SWC format to neuron$d
#' @examples
#' \donttest{
#' \dontrun{
#' # Plot the split to check it, correcting any errors, and saving your changes
#' hemibrain_adjust_saved_split()
#' }}
#' @export
#' @rdname hemibrain_adjust_saved_split
hemibrain_adjust_saved_split <- function(bodyids = NULL,
                          db = NULL,
                          check_thresh = 1,
                          batch_size = 5,
                          brain = NULL,
                          update_regularly = TRUE,
                          motivate = TRUE,
                          clean = TRUE,
                          Wait = TRUE,
                          sleep = 0.1,
                          selected_file = "1YjkVjokXL4p4Q6BR-rGGGKWecXU370D1YMc1mgUYr8E",
                          selected_col = "#fadadd",
                          ...){
  ### Get google sheet, assess what is to be done
  if(motivate){plot_inspirobot()}
  initials = readline(prompt = "What are your initials? ")
  say_hello(greet = initials)
  gs = googlesheets4::read_sheet(ss = selected_file, sheet = "roots")
  manual = googlesheets4::read_sheet(ss = selected_file, sheet = "manual")
  undone = subset(gs, checked < check_thresh)
  undone.ids = unique(undone$bodyid)
  undone.ids = sample(undone.ids, length(undone.ids), replace = FALSE)
  done = subset(gs, checked >= check_thresh)
  message(length(unique(done$bodyid)), "/",length(unique(gs$bodyid))," have been checked by at least ",check_thresh," user")
  users = gs[gs$checked,"user"]
  print(table(users))
  message("Manual splits by users: ")
  print(table(manual[!duplicated(manual$bodyid),]$user))
  ### Get chosen bodyIDs
  if(!is.null(bodyids)){
    message("Reading chosen neurons from the hemibrain project ...")
    bodyids = intersect(bodyids, undone.ids)
    batches = split(bodyids, ceiling(seq_along(bodyids)/batch_size))
  }else{
    batches = split(undone.ids, ceiling(seq_along(undone.ids)/batch_size))
  }
  if(!length(batches)){
    stop("None of the given bodyIDs have been checked by less than ",check_thresh," user")
  }else{
    for(batch in batches){
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
      notes = as.character(gs[match(names(someneuronlist), gs$bodyid),"note"][[1]])
      satisfied = FALSE
      while(!satisfied){
        ### Select neurons to edit and make notes on neurons
        message("Phase I : select neurons that need to be edited manually")
        phaseI <- splitcheck_phaseI(someneuronlist = someneuronlist,
                                    notes = notes,
                                    brain = brain,
                                    selected_col = selected_col)
        selected = phaseI["selected"][[1]]
        notes = phaseI["notes"][[1]]
        ### Manually fix selected skeletons
        message("Phase II : manually edit ", length(selected), " neurons")
        mes <- splitcheck_phaseII(someneuronlist = someneuronlist,
                                  selected = selected,
                                  brain = brain,
                                  manual = manual,
                                  update_regularly = update_regularly)
        ### Review the changes
        message("Phase II : review ", length(mes), " neurons")
        mes <- splitcheck_phaseIII(mes)
        ### Continue?
        satisfied = hemibrain_choice(prompt = "Continue to saving changes (y), or would you like to run the edit process
                                     for neurons in this batch again (keeping curent changes) (n)? ")
        someneuronlist[names(mes)] = mes
      }
      ### Save manual split to Google Sheet
      if(motivate){plot_inspirobot()}
      if(length(mes)){
        mes.sp = hemibrain_splitpoints(x = mes)
        mes.sp$split = "manual"
        mes.sp$checked = TRUE
        mes.sp$user = initials
        mes.sp$cut = unlist(sapply(mes[as.character(mes.sp$bodyid)], function(x) x$tags$cropped))
        mes.sp$soma = unlist(sapply(mes[as.character(mes.sp$bodyid)], function(x) x$tags$soma))
        mes.sp$time = Sys.time()
        mes.sp$note = ""
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
      rows = match(names(someneuronlist), undone$bodyid)
      splits = as.character(gs$split[rows])
      cuts = as.character(gs$cut[rows])
      somas = as.character(gs$soma[rows])
      if(length(mes)){
        which.mes = which(names(someneuronlist)%in%selected)
        splits[which.mes] = "manual"
        cuts[which.mes] = unlist(sapply(mes, function(x) x$tags$cropped))
        somas[which.mes] = unlist(sapply(mes, function(x) x$tags$soma))
      }
      update = data.frame(
        soma = somas,
        cut = cuts,
        split = splits,
        checked = gs[rows,]$checked+1,
        user = rep(initials, length(rows)),
        time = Sys.time(),
        note = notes)
      range = paste0("G",min(rows)+1,":M",max(rows)+1)
      update[is.na(update)] = ""
      gsheet_manipulation(FUN = googlesheets4::sheets_edit,
                          ss = selected_file,
                          range = range,
                          data = update[,intersect(colnames(gs),colnames(update))],
                          sheet = "roots",
                          col_names = FALSE)
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
    g = FUN(...)
    # g = tryCatch(FUN(...),
    #              error = function(e) NULL)
    if(!is.null(g)){
      success = TRUE
    }else{
      message("Google sheet read failures, re-trying ...")
    }
  }
}

# hidden
splitcheck_phaseI <- function(someneuronlist,
                              notes,
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
      rgl::clear3d()
      if(length(selected)){
        rgl::plot3d(someneuronlist[selected], col = hemibrain_bright_colour_ramp(length(selected)))
      }
      end = hemibrain_choice("Done selecting neurons to edit (shown)? yes/no ")
      rgl::clear3d()
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
    pl <- plot3d_split(someneuronlist[i])
    more_rgl_ids <- list()
    chc <- must_be(prompt = "Return to continue, b to go back, s to select, t to stop (with selection), c to cancel (no selection) and e to make a note: ",
                   answers = c("","b","s","t","c","e"))
    if (chc == "c" || chc == "t") {
      sapply(pl, rgl::rgl.pop, type = "shape")
      sapply(more_rgl_ids, rgl::rgl.pop, type = "shape")
      if (is.null(chc) || chc == "c")
        chc = must_be(prompt = "No selection made. Continue (c) or go back and make selection (b)? ", answers = c("c","b"))
      break
    }
    if (chc == "s") {
      if (n %in% selected) {
        message("Deselected: ", n)
        selected <- setdiff(selected, n)
      }
      else selected <- union(selected, n)
    }
    if (chc == "b")
      i <- i - 1
    if (chc == "e")
      notes[i] <- readline(prompt = "Add your note here: ")
    else i <- i + 1
    sapply(pl, rgl::rgl.pop, type = "shape")
    sapply(more_rgl_ids, rgl::rgl.pop, type = "shape")
    reset3d(brain=brain)
  }
  list(selected = selected,
       notes = notes
  )
}

# hidden
splitcheck_phaseII <- function(selected,
                               someneuronlist,
                               manual,
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
splitcheck_phaseIII <- function(mes, ...){
  if(length(mes)){
    save = hemibrain_choice(prompt = "Are you satisfied with your edits (y)? Or do you want to review them (n)? ")
    new.select = names(mes)
    while(!save){
      message("Select (s) which neurons to save")
      mes.sp = hemibrain_splitpoints(x = mes)
      mes2 = hemibrain_use_splitpoints(x = mes, df = mes.sp)
      new.select = nlscan_split(mes2)
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
