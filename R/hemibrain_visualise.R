#############################################################################
################################ visualise ##################################
#############################################################################

#' Plot neurons split up by flow centrality
#'
#' @param someneuronlist a neuronlist or neuron object that has been modified by flow.centrality
#' @param col colours of sections. Defaults to orange or axons, green for primary dendrite, blue for dendrites and pink for nodes with no flow.
#' @param splitnode if TRUE, a magenta sphere is placed at the location of the axon-dendrite split. Possible a putative action potential initiation site?
#' @param WithConnectors whether ot plot the anatomical location of pre (red) and post (cyan) synapses.
#' @param soma whether to plot a soma, and what the radius should be
#' @param WithNodes whether to plot branch points
#' @param lwd Line width (default 1)
#' @param radius For connectors and axon-dendrite split node (default 1)
#' @param highflow whether to plot the nodes of highest (with in one standard deviation less than maximum) flow centrality (pink points)
#' @param Verbose logical indicating that info about each selected neuron should be printed (default TRUE)
#' @param Wait logical indicating that there should be a pause between each displayed neuron
#' @param sleep time to pause between each displayed neuron when Wait=TRUE
#' @param extrafun an optional function called when each neuron is plotted, with two arguments: the current neuron name and the current selected neurons
#' @param selected_file an optional path to a yaml file that already contains a selection
#' @param selected_col the color in which selected neurons (such as those specified in selected_file) should be plotted
#' @param yaml a logical indicating that selections should be saved to disk in (human-readable) yaml rather than (machine-readable) rda format
#' @param ... additional arguments passed to methods.
#'
#' @return Plots coloured neuron(s)
#' @export
#' @seealso \code{\link{flow_centrality}}
#' @importFrom stats sd
plot3d_split = function(someneuronlist,
                        col = c("#1BB6AF", "#EF7C12", "#C70E7B", "#8FDA04", "#4D4D4D", "#FC6882"),
                        splitnode = FALSE, WithConnectors = TRUE, WithNodes = F, soma = 100, highflow = F, lwd = 1, radius = 1, ...){
  someneuronlist = nat::as.neuronlist(someneuronlist)
  for (n in 1:length(someneuronlist)){
    neuron = someneuronlist[[n]]
    if(is.null(neuron$d$flow.cent)){
      stop("No flow centrality calculated, dropping neuron")
    }
    dendrites.v = subset(rownames(neuron$d), neuron$d$Label == 3)
    axon.v = subset(rownames(neuron$d), neuron$d$Label == 2)
    #nulls.v = subset(rownames(neuron$d), neuron$d$Label == 0)
    p.d.v = subset(rownames(neuron$d), neuron$d$Label == 4)
    p.n.v = subset(rownames(neuron$d), neuron$d$Label == 7)
    dendrites = nat::prune_vertices(neuron, verticestoprune = as.integer(c(axon.v, p.d.v, p.n.v)))
    axon = nat::prune_vertices(neuron, verticestoprune = as.integer(c(dendrites.v, p.d.v, p.n.v)))
    #nulls = nat::prune_vertices(neuron, verticestoprune = as.integer(c(axon.v, dendrites.v, p.d.v, p.n.v)))
    p.d = nat::prune_vertices(neuron, verticestoprune = as.integer(c(axon.v, dendrites.v, p.n.v)))
    p.n = nat::prune_vertices(neuron, verticestoprune = as.integer(c(axon.v, dendrites.v, p.d.v)))
    rgl::plot3d(dendrites, col = col[1], WithNodes = WithNodes, lwd = lwd,...)
    rgl::plot3d(axon, col = col[2], WithNodes = WithNodes, soma = FALSE, lwd = lwd,...)
    rgl::plot3d(p.n, col = col[3], WithNodes = WithNodes, soma = soma, lwd = lwd,...)
    rgl::plot3d(p.d, col = col[4], WithNodes = WithNodes, soma = FALSE, lwd = lwd,...)
    #rgl::plot3d(nulls, col = col[5], WithNodes = WithNodes, soma = FALSE, lwd = lwd)
    #rgl::plot3d(neuron, col = col[3], WithNodes = WithNodes, soma = soma)
    if (WithConnectors){
      rgl::spheres3d(subset(nat::xyzmatrix(neuron$d),neuron$d$post>0), col = 'cyan', radius = radius,...)
      rgl::spheres3d(subset(nat::xyzmatrix(neuron$d),neuron$d$pre>0), col = 'red', radius = radius,...)
    }
    if (highflow == T){
      highest = max(neuron$d[,"flow.cent"])
      s.d = sd(neuron$d[,"flow.cent"], na.rm = T)
      high = subset(neuron$d, neuron$d[,"flow.cent"] > (highest - s.d))
      rgl::points3d(nat::xyzmatrix(high), col = col[6],...)
    }
    if(splitnode==T){
      ais = which(apply(neuron$d, 1, function(x) x["flow.cent"] == max(neuron$d[,"flow.cent"])))
      rgl::spheres3d(nat::xyzmatrix(neuron$d[ais,]),radius=radius,col="magenta",...)
    }
  }
}


#' @export
#' @rdname plot3d_split
nlscan_split <- function (someneuronlist, col = c("#1BB6AF", "#EF7C12", "#C70E7B", "#8FDA04", "#4D4D4D", "#FC6882"),
                          WithConnectors = T, WithNodes = F, soma = 100, highflow = F, Verbose = T, Wait = T,
                       sleep = 0.1, extrafun = NULL, selected_file = NULL, selected_col = "black",
                       yaml = TRUE, ...)
{
  if(!requireNamespace('yaml', quietly = TRUE))
    stop("Suggested package yaml is required to use this function!")
  if (nat::is.neuronlist(someneuronlist)) {
    db = someneuronlist
    neurons = as.data.frame(db)$name
  }
  frames <- length(neurons)
  selected <- character()
  i <- 1
  if (!is.null(selected_file) && file.exists(selected_file)) {
    selected <- yaml::yaml.load_file(selected_file)
    if (!all(names(selected) %in% neurons))
      stop("Mismatch between selection file and neurons.")
  }
  savetodisk <- function(selected, selected_file) {
    if (is.null(selected_file))
      selected_file <- file.choose(new = TRUE)
    if (yaml) {
      if (!grepl("\\.yaml$", selected_file))
        selected_file <- paste(selected_file, sep = "",
                               ".yaml")
      message("Saving selection to disk as ", selected_file,
              ".")
      writeLines(yaml::as.yaml(selected), con = selected_file)
    }
    else {
      if (!grepl("\\.rda$", selected_file))
        selected_file <- paste(selected_file, sep = "",
                               ".rda")
      save(selected, file = selected_file)
      message("Saving selection to disk as ", selected_file)
    }
    selected_file
  }
  chc <- NULL
  while (TRUE) {
    if (i > length(neurons) || i < 1)
      break
    n <- neurons[i]
    cat("Current neuron:", n, "(", i, "/", length(neurons),
        ")\n")
    pl <- plot3d_split(someneuronlist[i], col = col, WithConnectors = WithConnectors, WithNodes = WithNodes, soma = soma, highflow = highflow)
    message("segregation index: ", someneuronlist[[i]]$segregation.index)
    more_rgl_ids <- list()
    if (!is.null(extrafun))
      more_rgl_ids <- extrafun(n, selected = selected)
    if (Wait) {
      chc <- readline("Return to continue, b to go back, s to select, d [save to disk], t to stop, c to cancel (without returning a selection): ")
      if (chc == "c" || chc == "t") {
        sapply(pl, rgl::rgl.pop, type = "shape")
        sapply(more_rgl_ids, rgl::rgl.pop, type = "shape")
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
      else if (chc == "d")
        savetodisk(selected, selected_file)
      else i <- i + 1
    }
    else {
      Sys.sleep(sleep)
      i <- i + 1
    }
    sapply(pl, rgl::rgl.pop, type = "shape")
    sapply(more_rgl_ids, rgl::rgl.pop, type = "shape")
    rgl::clear3d()
  }
  if (is.null(chc) || chc == "c")
    return(NULL)
  if (!is.null(selected_file))
    savetodisk(selected, selected_file)
  selected
}
