################################################################################
################################ Manual Split ##################################
################################################################################

#' Manually assign the dendrite and axon to a neuron
#'
#' @description 3D interactive functions to edit neurons in an rgl viewer, including re-assigned cable as
#' being axon, dendrite, rimary neurite, etc., choosing soma locations and pruning bits of
#' cable as a user chooses.
#'
#' @param x a neuron/neuronlist object
#' @param brain which brain to plot while splitting neuron.
#' @param Label the type of cable to be assigned/edited. Use \code{standardise_labels} to see how numbers convert to cable types.
#' @param ... additional arguments passed to plotting functions, namely \code{plot3d_split}
#' @return The neuron/neuronlist object with axon/dendrite info assigned in SWC format to neuron$d
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
manually_assign_labels <-function(x, brain = NULL, ...) UseMethod("manually_assign_labels")

#' @export
#' @rdname manually_assign_labels
manually_assign_labels.neuron <- function(x, brain = NULL, ...){
  bodyid = ifelse(is.null(x$bodyid),x$skid,x$bodyid)
  reset3d(brain=brain);plot3d_split(x)
  x$tags$manually_edited = FALSE
  x$tags$cropped = FALSE
  x$tags$soma = "automatic"
  happy = hemibrain_engage(Label = NULL)
  if(happy){
    happy = hemibrain_engage(prompt = "Are you sure this neuron has no problems? yes/no", Label = NULL)
  }
  ### Cropped?
  if(!happy){
    crop = TRUE
    while(crop){
      crop = hemibrain_choice("Is this neuron cropped so badly it is mising a big chunk of axon/dendrite? yes/no ")
      x$tags$cropped = crop
      if(crop){
        crop2 = must_be("Can you still split the neuron into axon and/or dendrite? yes/no ")
        x$tags$cropped = ifelse(crop2,"MB",FALSE)
      }
      crop = !hemibrain_choice("Happy with this answer? yes/no ")
    }
    happy = hemibrain_engage(Label = NULL, prompt = "Do you want to edit further? yes/no ")
  }
  ### Manually edit
  while(!happy%in%c("y","yes",TRUE)){
    reset3d(brain=brain);plot3d_split(x)
    happy = hemibrain_edit_cable()
    if(happy==0){
      break
    }else{
      if(happy == 1){
        ### Soma
        message("Please check the soma of your neuron")
        x = hemibrain_correctsoma(x, ...)
        x$tags$manually_edited = TRUE
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
          happy = !hemibrain_engage(Label = NULL, prompt = "Do you want to edit further? yes/no ")
        }
      }else{
        x = hemibrain_select_cable(x=x, Label = happy)
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
manually_assign_labels.neuronlist<-function(x, brain = NULL, ...){
  nat::nlapply(x, manually_assign_labels.neuron, brain = brain, ...)
}

#' @export
#' @rdname manually_assign_labels
hemibrain_prune_online <-function(x, ...) UseMethod("hemibrain_prune_online")

#' @export
#' @rdname manually_assign_labels
hemibrain_prune_online.neuron <- function (x, brain = NULL, Label = NULL, ...) {
  continue = "no"
  ids = integer()
  while (!continue %in% c("y", "yes")) {
    pids = plot3d_split(x, ...)
    selected = hemibrain_select_points(x$d, clear_plot_on_exit = TRUE, Label = Label)
    reset3d(brain=brain)
    if(nrow(selected)){
      message("Selected cable in cerise")
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
  frag
}

#' @export
#' @rdname manually_assign_labels
hemibrain_prune_online.neuronlist <- function (x, Label = NULL, ...){
  nat::nlapply(x, hemibrain_prune_online.neuron, Label = Label, ...)
}

# Hidden
hemibrain_select_points <- function(points, clear_plot_on_exit = FALSE, Label = NULL, ...) {
  if(is.null(Label)){
    selected.points <- data.frame()
  }else{
    selected.points <- nat::xyzmatrix(points[points$Label==Label,])
  }
  points <- nat::xyzmatrix(points)
  ids = rgl::points3d(selected.points, col = hemibrain_bright_colours["cerise"])
  ids = c(ids, rgl::points3d(points, col = "grey30"))
  progress = must_be(prompt = "Unselected points in grey, selected in cerise. Add (a) or remove (r) points, or continue (c)?  ",
                     answers = c("a","r","c"))
  while (progress != "c") {
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
    ids = union(ids, rgl::points3d(points, col = "grey30", ... ))
    progress = must_be(prompt = "Add (a) or remove (r) points, or continue (c)?  ",
                        answers = c("a","r","c"))
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
      progress = must_be(prompt = "Good enough? yes/no  ", answers = c("y","yes","n","no"))
  }
  x = add_Label(x = x, PointNo = selected.point, Label = 1, erase = TRUE)
  corrected$tags$soma = "manual"
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
hemibrain_engage <- function(Label = NULL, prompt = "Does this neuron require editing, or is cropped? yes/no ", ...){
  if(is.null(Label)){
    !hemibrain_choice(prompt = prompt)
  }else{
    comp = standard_compartments(Label)
    !hemibrain_choice(paste("Is the",comp,"looking good? yes/no ",sep = " "))
  }
}

# hidden
hemibrain_edit_cable <- function(){
  as.numeric(must_be(prompt = "What cable do you wish to edit?:
                     soma (1), axon (2), dendrite (3),
                     linker (4), primary neurite (7),
                     invert the neuron (8) or do nothing further (0) ",
          answers  = c(0:4,7:8)))
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
hemibrain_select_cable <- function(x, Label, ...){
  if(hemibrain_engage(Label = Label)){
    comp = standard_compartments(Label)
    message(paste0("Please select: ",comp))
    frag = hemibrain_prune_online(x, Label = Label, ...)
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
  dendrites = subset(rownames(x$d), x$d$Label == 3)
  axon = subset(rownames(x$d), x$d$Label == 2)
  p.d = subset(rownames(x$d), x$d$Label == 4)
  p.n = subset(rownames(x$d), x$d$Label == 7)
  nulls = subset(rownames(x$d), !x$d$Label %in% c(2,3))
  ## Get possible cable change points
  s.d = change_points(x = x, v = dendrites)
  s.a = change_points(x = x, v = axon)
  s.n = change_points(x = x, v = nulls)
  ## Get cable start points
  d.starts = sapply(unique(c(s.n,s.a)),function(s) igraph::neighbors(n, v=s, mode = c("out")))
  d.starts = as.character(unique(unlist(d.starts)))
  a.starts = sapply(unique(c(s.n,s.d)),function(s) igraph::neighbors(n, v=s, mode = c("out")))
  a.starts = as.character(unique(unlist(a.starts)))
  axon.starts = intersect(axon, a.starts)
  dendrites.starts = intersect(dendrites, d.starts)
  ## Assign primaries
  possible = c(p.d[1],p.d[length(p.d)])
  starts = lapply(possible,function(s) igraph::neighbors(n, v=s, mode = c("all")))
  starts.p = sapply(starts, function(start) sum(start%in%axon))
  dendrite.primary = possible[which.min(starts.p)] ## not necessarily strictly correct ...
  axon.primary = possible[which.max(starts.p)] ## not necessarily strictly correct ...
  ## Get primary branch point
  primary.branch.point = p.n[length(p.n)]
  ## Assign
  x$primary.branch.point = nullToNA(primary.branch.point)
  x$axon.start = nullToNA(axon.starts)
  x$dendrite.start = nullToNA(dendrites.starts)
  x$axon.primary = nullToNA(axon.primary)
  x$dendrite.primary = nullToNA(dendrite.primary)
  x
}
