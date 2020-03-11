#############################################################################
################################ visualise ##################################
#############################################################################

#' Plot neurons split up by flow centrality
#'
#' @param someneuronlist a neuronlist or neuron object that has been modified by flow.centrality
#' @param col colours of sections. Defaults to orange or axons, green for primary dendrite, blue for dendrites and pink for nodes with no flow.
#' @param splitnode if TRUE, a magenta sphere is placed at the location of the axon-dendrite split. Possible a putative action potential initiation site?
#' @param WithConnectors whether to plot the anatomical location of pre (red) and post (cyan) synapses.
#' @param soma whether to plot a soma, and what the radius should be. If \code{NULL}, an appropriate value is guessed.
#' @param WithNodes whether to plot branch points
#' @param lwd Line width (default 1)
#' @param radius For connectors and axon-dendrite split node (default 1). If \code{NULL}, an appropriate value is guessed.
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
#' @importFrom nat xyzmatrix
plot3d_split = function(someneuronlist,
                        col = c("#1BB6AF", "#EF7C12", "#C70E7B", "#8FDA04", "#4D4D4D", "#FC6882"),
                        splitnode = FALSE,
                        WithConnectors = TRUE,
                        WithNodes = F,
                        soma = NULL,
                        highflow = FALSE,
                        lwd = 1,
                        radius = NULL,
                        ...){
  someneuronlist = nat::as.neuronlist(someneuronlist)
  temps = nat.templatebrains::all_templatebrains()
  temps.microns = c(temps[temps$W<2000,"name"],"JRCFIB2018F")
  reg = nat.templatebrains::regtemplate(someneuronlist)
  if(is.null(soma)){
    if(is.null(reg)|!reg%in%temps.microns){
      soma = 500
    }else{
      soma = 4
    }
  }
  if(is.null(radius)){
    if(is.null(reg)|!reg%in%temps.microns){
      radius = 100
    }else{
      radius = 8
    }
  }

  for (n in 1:length(someneuronlist)){
    neuron = someneuronlist[[n]]
    dendrites.v = subset(rownames(neuron$d), neuron$d$Label == 3)
    axon.v = subset(rownames(neuron$d), neuron$d$Label == 2)
    p.d.v = subset(rownames(neuron$d), neuron$d$Label == 4)
    p.n.v = subset(rownames(neuron$d), neuron$d$Label == 7)
    dendrites = tryCatch( nat::prune_vertices(neuron, verticestoprune = as.integer(c(axon.v, p.d.v, p.n.v))),
                          error = function(e) NULL)
    axon = tryCatch( nat::prune_vertices(neuron, verticestoprune = as.integer(c(dendrites.v, p.d.v, p.n.v))),
                     error = function(e) NULL)
    p.d = tryCatch( nat::prune_vertices(neuron, verticestoprune = as.integer(c(axon.v, dendrites.v, p.n.v))),
                    error = function(e) NULL)
    p.n = tryCatch( nat::prune_vertices(neuron, verticestoprune = as.integer(c(axon.v, dendrites.v, p.d.v))),
                    error = function(e) NULL)
    tryCatch(rgl::plot3d(dendrites, col = col[1], WithNodes = WithNodes, lwd = lwd, add = TRUE, ...),
             error = function(e) NULL)
    tryCatch(rgl::plot3d(axon, col = col[2], WithNodes = WithNodes, soma = FALSE, lwd = lwd, add = TRUE, ...),
             error = function(e) NULL)
    tryCatch(rgl::plot3d(p.n, col = col[3], WithNodes = WithNodes, soma = soma, lwd = lwd, add = TRUE, ...),
    error = function(e) NULL)
    tryCatch(rgl::plot3d(p.d, col = col[4], WithNodes = WithNodes, soma = FALSE, lwd = lwd, add = TRUE,...),
    error = function(e) NULL)
    tryCatch(rgl::plot3d(neuron, col = col[3], WithNodes = WithNodes, soma = soma, add = TRUE,),
    error = function(e) NULL)
    if (WithConnectors){
      conns=neuron$connectors
      tryCatch(rgl::spheres3d(xyzmatrix(conns[conns$prepost==1,,drop=FALSE]),
                              col = "#132157", radius = radius/2, add = TRUE, ...),
               error = function(e) NULL)
      tryCatch(rgl::spheres3d(xyzmatrix(conns[conns$prepost==0,,drop=FALSE]),
                              col = "#EE4244", radius = radius, add = TRUE, ...),
               error = function(e) NULL)
    }
    if (highflow == T){
      highest = max(neuron$d[,"flow.cent"])
      s.d = sd(neuron$d[,"flow.cent"], na.rm = T)
      high = subset(neuron$d, neuron$d[,"flow.cent"] > (highest - s.d))
      rgl::points3d(nat::xyzmatrix(high), col = col[6],...)
    }
    if(splitnode==T){
      ais = which(apply(neuron$d, 1, function(x) x["flow.cent"] == max(neuron$d[,"flow.cent"])))
      rgl::spheres3d(nat::xyzmatrix(neuron$d[ais,]),radius=radius,col="magenta", add = TRUE, ...)
    }
  }
}


#' @export
#' @rdname plot3d_split
nlscan_split <- function (someneuronlist, col = c("#1BB6AF", "#EF7C12", "#C70E7B", "#8FDA04", "#4D4D4D", "#FC6882"),
                          WithConnectors = TRUE, WithNodes = FALSE, soma = NULL, radius = NULL, highflow = FALSE, Verbose = TRUE, Wait = TRUE,
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
    pl <- plot3d_split(someneuronlist[i], col = col, WithConnectors = WithConnectors, WithNodes = WithNodes, soma = soma, highflow = highflow, radius = radius, ...)
    message("segregation index: ", someneuronlist[[i]]$AD.segregation.index)
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
