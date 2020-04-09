################################################################################
################################ Manual Split ##################################
################################################################################

#' Manually assign the dendrite and axon to a neuron
#'
#' @description 3D interactive functions to edit neurons in an rgl viewer,
#'   including re-assigned cable as being axon, dendrite, primary neurite, etc.,
#'   choosing soma locations and pruning bits of cable as a user chooses.
#'
#' @param x a neuron/neuronlist object
#' @param brain which brain to plot while splitting neuron.
#' @param is.cropped logical. If TRUE, the user is asked whether the neuron being considered is complete or not.
#' @param Label the type of cable to be assigned/edited. Use \code{standardise_labels} to see how numbers convert to cable types.
#' @param lock a vector of number representing cable types to lock (so that they cannot be edited interactively be a user).
#' If \code{NULL} (default) no cable type is treated as locked. Locked cable will appear in green.
#' @param ... additional arguments passed to plotting functions, namely \code{plot3d_split}
#' @return The neuron/neuronlist object with axon/dendrite info assigned in SWC format to \code{neuron$d}
#' @examples
#' \donttest{
#' # Choose some neurons
#'exemplars = c("202916528", "1279775082",  "203253072",
#'              "326530038",  "203253253", "5813079341")
#'
# Get neurons
#'neurons = neuprintr::neuprint_read_neurons(exemplars)
#'
# Now use a pre-saved axon-dendrite split
#'neurons.flow = hemibrain_flow_centrality(neurons)
#'
#' \dontrun{
#' # Plot the split to check it, correcting any errors
#' nat::nopen3d()
#' neurons.manual = manually_assign_labels(neurons.flow)
#' }}
#' @export
#' @rdname manually_assign_labels
manually_assign_labels <-function(x, brain = NULL, is.cropped = FALSE, ...) UseMethod("manually_assign_labels")

#' @export
#' @rdname manually_assign_labels
manually_assign_labels.neuron <- function(x, brain = NULL, is.cropped = FALSE, ...){
  bodyid = ifelse(is.null(x$bodyid),x$skid,x$bodyid)
  WithConnectors = TRUE
  reset3d(brain=brain);plot3d_split(x, WithConnectors = WithConnectors)
  x$tags$manually_edited = FALSE
  x$tags$cut = FALSE
  x$tags$soma = "automatic"
  x.safe = x
  happy = hemibrain_engage(Label = NULL)
  if(happy){
    happy = hemibrain_engage(prompt = "Are you sure this neuron has no problems? yes/no ", Label = NULL)
  }
  ### Cropped?
  if(!happy & is.cropped){
    cropped = hemibrain_iscropped()
    x$tags$cut = cropped
    crop = TRUE
    happy = hemibrain_engage(Label = NULL, prompt = "Do you want to edit further? yes/no ")
  }
  ### Manually edit
  locks = c()
  while(!happy%in%c("y","yes",TRUE)){
    reset3d(brain=brain)
    plot3d_split(x, WithConnectors = WithConnectors)
    if(!is.issue(locks)){
      message("locks: ", paste(locks,collapse = ", "))
    }
    happy = hemibrain_edit_cable()
    if(happy==0){
      break
    }else{
      if(happy == 1){
        ### Soma
        message("Please check the soma of your neuron")
        x = hemibrain_correctsoma(x, ...)
      }else if(happy == 8){
        ### Invert?
        if(sum(c(2,3)%in%unique(x$d$Label))>0){
          i = TRUE
          while(i){
            x = hemibrain_invert(x=x)
            reset3d(brain=brain)
            plot3d_split(x)
            i = hemibrain_choice("Should we invert the dendrite and axon assignment again? yes/no ")
          }
          x$tags$manually_edited = i
        }
      }else if(happy == 9){
        x = cycle_branches(x=x,brain=brain)
        x$tags$manually_edited = TRUE
      }else if(happy==10){
        x = x.safe
        x$tags$manually_edited = FALSE
        message("Neuron has been reverted to its original state")
      }else if(happy == 11){
        WithConnectors = !WithConnectors
        message("Plot synapses: ", WithConnectors)
      }else if(happy == 12){
        x = hemibrain_makenote(x, ...)
      }else if (grepl("L",happy)){
        if(happy %in% locks){
          locks = setdiff(locks,happy)
        }else{
          locks = unique(c(locks,happy))
        }
      }else{
        x = hemibrain_select_cable(x=x, Label = happy, lock = locks)
        x$tags$manually_edited = TRUE
      }
    }
  }
  x = internal_assignments(x)
  x$bodyid = bodyid
  x = hemibrain_neuron_class(x)
  x
}

#' @export
#' @rdname manually_assign_labels
manually_assign_labels.neuronlist<-function(x, brain = NULL, is.cropped = FALSE, ...){
  nat::nlapply(x, manually_assign_labels.neuron, brain = brain, ...)
}

#' @export
#' @rdname manually_assign_labels
hemibrain_prune_online <-function(x, brain = NULL, Label = NULL, lock = NULL, ...) UseMethod("hemibrain_prune_online")

#' @export
#' @rdname manually_assign_labels
hemibrain_prune_online.neuron <- function (x, brain = NULL, Label = NULL, lock = NULL, ...) {
  if(lockdown(Label = Label, lock = lock, verbose = TRUE)){
    return(NULL)
  }
  continue = "no"
  ids = integer()
  while (!continue %in% c("y", "yes")) {
    pids = plot3d_split(x, ...)
    selected = hemibrain_select_points(x$d, clear_plot_on_exit = TRUE, Label = Label, lock = lock)
    reset3d(brain=brain)
    if(nrow(selected)){
      message("Selected cable in red")
      v = match(data.frame(t(selected)), data.frame(t(nat::xyzmatrix(x))))
      x = add_Label(x = x, PointNo = v, Label = Label, erase = TRUE)
      frag = nat::prune_vertices(x, verticestoprune = v, invert = TRUE, ...)
      ids = rgl::plot3d(nat::as.neuronlist(frag), col = hemibrain_bright_colours["cerise"], ...)
    }else{
      message("Note: no points selected")
      frag = NULL
    }
    ids = c(ids, rgl::plot3d(nat::as.neuronlist(x), col = "grey30", ...))
    continue = must_be("Finished with this selection? yes/no ", answers = c("y","yes","n","no"))
    rgl::pop3d(id = unlist(ids))
  }
  frag = carryover_tags(x = x, y = frag)
  frag
}

#' @export
#' @rdname manually_assign_labels
hemibrain_prune_online.neuronlist <- function (x, brain = NULL, Label = NULL, lock = NULL, ...){
  nat::nlapply(x, hemibrain_prune_online.neuron, brain = brain, Label = Label, lock = lock, ...)
}

# Hidden
hemibrain_select_points <- function(x, clear_plot_on_exit = FALSE, Label = NULL, lock = NULL, ...) {
  if(nat::is.neuron(x)){
    x = x$d
  }
  if(is.null(Label)){
    selected.points <- data.frame()
  }else{
    selected.points <- nat::xyzmatrix(x[x$Label==Label,])
  }
  if(lockdown(Label=Label, lock = lock)){
    return(selected.points)
  }
  if(!is.null(lock)){
    points.locked <- nat::xyzmatrix(x[x$Label%in%lock,])
    x <- x[!x$Label%in%lock,]
  }else{
    points.locked <- NULL
  }
  points <- nat::xyzmatrix(x)
  ids = rgl::points3d(selected.points, col = hemibrain_bright_colours["cerise"])
  ids = c(ids, rgl::points3d(points.locked, col = hemibrain_bright_colours["green"]))
  ids = c(ids, rgl::points3d(points, col = "grey30"))
  progress = "progressive"
  while (progress != "c") {
    progress = must_be(prompt = "Add (a) or remove (r) points, or continue (c)?  ",
                       answers = c("a","r","c"))
    if (progress == "a") {
      keeps = rgl::select3d()
      keep.points <- keeps(points)
      keep.points = subset(points, keep.points)
      selected.points = rbind(selected.points, keep.points)
    }
    if (progress == "r") {
      remove.points <- rgl::select3d()
      removed.points <- remove.points(selected.points)
      selected.points = subset(selected.points, !removed.points)
    }
    rgl::pop3d(id = ids)
    ids = integer()
    if (length(selected.points) > 0) {
      ids = rgl::points3d(selected.points, col = hemibrain_bright_colours["cerise"], ...)
    }
    ids = union(ids, rgl::points3d(points.locked, col = hemibrain_bright_colours["green"], ... ))
    ids = union(ids, rgl::points3d(points, col = "grey30", ... ))
  }
  if (clear_plot_on_exit)
    rgl::pop3d(id = ids)
  return(selected.points)
}

#' @export
#' @rdname manually_assign_labels
hemibrain_correctsoma <- function(x, ...) {
  plot3d_split(x, ...)
  eps.xyz = x$d[c(nat::endpoints(x), nat::branchpoints(x)), ]
  progress = "no"
  is.there.a.soma = hemibrain_choice(prompt = "Is there even a soma for this neuron? yes/no ")
  if(is.there.a.soma){
    if(is.null( x$tags$soma)){
    x$tags$soma = "automatic"
    }
  }else{
    x$tags$soma = "missing"
    message("Soma recorded as missing. Please still choose a sensible 'root' point ")
  }
  while(progress %in% c("n", "no")) {
      message("Rotate brain and then hit [enter] to continue")
      line <- readline()
      message("Select new root from highlighted end/branch points")
      selected.point <- rgl::select3d()
      selected.point <- selected.point(nat::xyzmatrix(eps.xyz))
      selected.point <- eps.xyz$PointNo[selected.point]
      if (length(selected.point) > 1 | length(selected.point) == 0) {
        message("Multiple end points selected, try again ...")
      } else if (length(selected.point)) {
        corrected = nat::as.neuron(nat::as.ngraph(x), origin = selected.point)
        corrected = hemibrain_carryover_labels(y = corrected, x = x)
        reset3d(...)
        plot3d_split(corrected, ...)
      }
      x$tags$soma = "manual"
      x$tags$manually_edited = TRUE
      progress = must_be(prompt = "Good enough? yes/no  ", answers = c("y","yes","n","no"))
  }
  x = add_Label(x = x, PointNo = selected.point, Label = 1, erase = TRUE)
  corrected = carryover_tags(x = x, y = corrected)
  corrected
}

# Hidden
hemibrain_carryover_labels <- function(x, y){
  for(label in unique(x$d$Label)){
    y = add_Label(x = y, PointNo = x$d[x$d$Label==label,]$PointNo, Label = label, erase = FALSE)
  }
  y
}

# hidden
reset3d <- function(brain = NULL){
  rgl::clear3d()
  if(!is.null(brain)){
    rgl::plot3d(brain, col = "grey70", alpha = 0.1)
  }
}

# hidden
hemibrain_engage <- function(Label = NULL, prompt = "Does this neuron require editing? yes/no ", ...){
  if(is.null(Label)){
    !hemibrain_choice(prompt = prompt)
  }else{
    comp = standard_compartments(Label)
    !hemibrain_choice(paste("Is the",comp,"looking good? yes/no ",sep = " "))
  }
}

# hidden
hemibrain_edit_cable <- function(){
  options = c(0:4,7:12)
  cables = c(0:7)
  locks = c(paste0("L",cables),paste0("l",cables))
  choice = must_be(prompt = "What cable do you wish to edit?:
  soma (1), axon (2), dendrite (3), linker (4), primary neurite (7),
  invert the neuron (8) or cycle through branches (9), revert all edits (10),
  toggle synapses (11), make a note (12) do nothing further (0).
  L+number to lock cable, ",
          answers  = c(options,locks))
  choice = gsub("l","L",choice)
}

# hidden
hemibrain_choice <- function(prompt){
  choice = must_be(prompt = prompt, answers = c("y","yes","n","no"))
  if(choice %in% c("y","yes")){
    TRUE
  }else{
    FALSE
  }
}

# hidden
hemibrain_select_cable <- function(x, Label, lock, ...){
  locks = gsub("L","",lock)
  if(lockdown(Label = Label, lock = locks)){
    return(x)
  }else if(hemibrain_engage(Label = Label)){
    comp = standard_compartments(Label)
    message(paste0("Please select: ",comp))
    frag = hemibrain_prune_online(x, Label = Label, lock = locks, ...)
    if(is.null(frag)){
      pnos = "DELETE"
    }else{
      pnos = frag$d$PointNo
    }
    x = add_Label(x = x, PointNo = pnos, Label = Label, erase = TRUE)
  }
  x
}

# hidden
hemibrain_invert <- function(x){
  y = x
  y = add_Label(x = y, PointNo = x$d$PointNo[x$d$Label==2], Label = 3)
  y = add_Label(x = y, PointNo = x$d$PointNo[x$d$Label==3], Label = 2)
  y
}

# hidden
## re-calculate axon/dendrite startpoints from label locations in the neuron.
## Useful if the labels on a neuron have been modified since using flow_centrality function.
internal_assignments <- function(x){
  ## Get cable points
  n = nat::as.ngraph(x)
  root = nat::rootpoints(n)
  dendrites = subset(rownames(x$d), x$d$Label == 3)
  dendrites = setdiff(dendrites, root)
  axon = subset(rownames(x$d), x$d$Label == 2)
  axon = setdiff(axon, root)
  p.d = subset(rownames(x$d), x$d$Label == 4)
  p.n = subset(rownames(x$d), x$d$Label == 7)
  nulls = subset(rownames(x$d), !x$d$Label %in% c(2,3))
  ## Get possible cable change points
  if(length(nulls)){
    s.n = change_points(x = x, v = nulls)
  }else{
    s.n = NULL
  }
  if(length(dendrites)){
    s.d = change_points(x = x, v = dendrites)
  }else{
    s.d = NULL
  }
  if(length(axon)){
    s.a = change_points(x = x, v = axon)
  }else{
    s.a = NULL
  }
  ## Get cable start points
  d.starts = tryCatch(sapply(s.d,function(s) igraph::neighbors(n, v=s, mode = c("in"))), error = function(e) NA)
  d.starts = tryCatch(sapply(d.starts,function(s) !s%in%dendrites), error = function(e) NA)
  a.starts = tryCatch(sapply(s.a,function(s) igraph::neighbors(n, v=s, mode = c("in"))), error = function(e) NA)
  a.starts = tryCatch(sapply(a.starts,function(s) !s%in%axon), error = function(e) NA)
  axon.starts = s.a[unlist(a.starts)]
  dendrites.starts = s.d[unlist(d.starts)]
  ## Assign primaries
  if(length(p.d)){
    possible = c(p.d[1],p.d[length(p.d)])
  }else{
    possible = ""
  }
  pd.dists = tryCatch(igraph::distances(n, v = p.d, to = as.numeric(p.d), mode = c("all")),
                     error = function(e) NA)
  linkers = suppress(tryCatch(rownames(which(pd.dists == max(pd.dists), arr.ind = TRUE)),error = function(e) NA))
  ### Primary (furthest from root) cable starts
  d.dists = tryCatch(igraph::distances(n, v = root, to = as.numeric(dendrites.starts), mode = c("all")),
                  error = function(e) NA)
  dendrite.primary = tryCatch(dendrites.starts[which.max(d.dists)], error = function(e) NA)
  a.dists = tryCatch(igraph::distances(n, v = root, to = as.numeric(axon.starts), mode = c("all")),
                     error = function(e) NA)
  axon.primary = tryCatch(axon.starts[which.max(a.dists)], error = function(e) NA)
  ## Get primary branch point
  dists = tryCatch(igraph::distances(n, v = root, to = as.numeric(p.n), mode = c("all")),
                   error = function(e) NA)
  primary.branch.point = tryCatch(p.n[which.max(dists)], error = function(e) NA)
  ## Assign
  x$primary.branch.point = purify(primary.branch.point)
  x$axon.start = purify(axon.starts)
  x$dendrite.start = purify(dendrites.starts)
  x$axon.primary = purify(axon.primary)
  x$dendrite.primary = purify(dendrite.primary)
  x$linker = purify(linkers)
  x
}

# hidden
#' @importFrom rgl plot3d
cycle_branches <- function(x, brain = NULL){
  message("Removing assigned primary neurite and linker cable to find sub-branches ...")
  sbt = break_into_subtrees(x, prune = TRUE)
  i <- 1
  message("Cycling through sub-branches (in red). Bear in mind that some branches can be very small.")
  while(length(sbt)) {
    if (i > length(sbt) || i < 1){
      end = hemibrain_choice("Done selecting branches to edit? yes/no ")
      if(end){
        return(x)
      }else{
        i <- 1
      }
    }
    reset3d(brain=brain)
    plot3d_split(x, soma = FALSE)
    cat("Current branch:", i, " "," (", i, "/", length(sbt),")","\n")
    plot3d(sbt[i], col = hemibrain_bright_colours["cerise"], lwd = 5)
    chc <- must_be(prompt = "What cable is this?:
    uncertain/erroneous (0), axon (2), dendrite (3),
    linker (4), primary neurite (7).
    Back (b), next (enter), cancel (c) ",
          answers = c(0,2:4,7,"","b","c"))
    if (chc == "c") {
      if (is.null(chc) || chc == "")
        break
    }
    if (chc %in%  c(0,2:4,7) ){
      message("Updating ", standard_compartments(chc))
      x = add_Label(x = x, PointNo = sbt[i][[1]]$orig.PointNo, Label = as.numeric(chc), erase = FALSE)
    }
    if (chc == "b"){
      i <- i - 1
    }else{
      i <- i + 1
    }
  }
  x
}

#hidden
lockdown <- function(Label = NULL, lock = NULL, verbose = FALSE){
  if(!is.issue(lock)&!is.issue(Label)){
    lock = setdiff(lock,Label)
    if(length(lock)==0){
      if(verbose){message("Label and lock cannot be the same value")}
      TRUE
    }else{
      if(verbose){message("Cable type locked: ", paste(standard_compartments(lock),collapse=", "))}
      FALSE
    }
  }else{
    FALSE
  }
}

# hidden
hemibrain_iscropped <- function(){
  crop = TRUE
  while(crop){
    cropped = crop = hemibrain_choice("Is this neuron cropped so badly it is mising a big chunk of axon/dendrite? yes/no ")
    if(crop){
      crop2 = hemibrain_choice("Can you still split the neuron into axon and/or dendrite? yes/no ")
      cropped = ifelse(crop2,"MB",FALSE)
    }
    message("The cropped status if this neuron will be recorded as: ", cropped, "  (options are TRUE/FALSE/slightly(MB))")
    crop = !hemibrain_choice("Happy with this answer? yes/no ")
  }
  cropped
}

