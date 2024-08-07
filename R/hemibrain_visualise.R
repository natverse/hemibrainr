#############################################################################
################################ visualise ##################################
#############################################################################

#' Plot neurons split up by flow centrality
#'
#' @param someneuronlist a neuronlist or neuron object that has been modified by flow.centrality.
#' @param already_selected a vector of IDs for neurons that are already considered selected. If \code{NULL} then no neurons are considered selected to start with.
#' @param col colours of sections. Defaults to orange or axons, green for primary dendrite, blue for dendrites and pink for nodes with no flow.
#' @param splitnode if TRUE, a magenta sphere is placed at the location of the axon-dendrite split. Possible a putative action potential initiation site?
#' @param WithConnectors whether to plot the anatomical location of pre (red) and post (cyan) synapses.
#' @param prepost whether to plot input (POST-) or output (PRE-) synapses, or both (default).
#' @param soma whether to plot a soma, and what the radius should be. If \code{NULL}, an appropriate value is guessed.
#' @param soma.alpha numeric, alpha transparency value for plotting the soma. Passed to \code{rgl::spheres3d}
#' @param WithNodes whether to plot branch points.
#' @param lwd Line width (default 1).
#' @param radius For connectors and axon-dendrite split node (default 1). If \code{NULL}, an appropriate value is guessed.
#' @param brain a template brain to plot. \code{FALSE} results in no brain plotted.
#' @param highflow whether to plot the nodes of highest (with in one standard deviation less than maximum) flow centrality (pink points).
#' @param transmitters logical. If \code{TRUE}, and transmitter identities are given in each neuron's meta data at neuron$connectors, then synapses.
#' are plotted in their transmitter colours. Transmitter identity is given by the colour of a halo around each synapse. The colours correspond to entries in \code{\link{paper_colours}}.
#' @param transmitter.alpha numeric, alpha transparency value for a coloured halo indicative of transmitter identity around synapses. Passed to \code{rgl::spheres3d}.
#' @param volume a \code{mesh3d} or \code{hxsurf} object. Only somas outside this volume will be plotted.
#' @param check.template check which template space \code{someneuronlist} is in, in order to set default plotting settings.
#' @param invert logical, if \code{TRUE} only somas outside \code{volume} will be plotted, if \code{FALSE}, only those inside.
#' @param Verbose logical indicating that info about each selected neuron should be printed (default TRUE)
#' @param Wait logical indicating that there should be a pause between each displayed neuron
#' @param sleep time to pause between each displayed neuron when Wait=TRUE
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
plot3d_split <- function(someneuronlist,
                        prepost = c("BOTH","PRE","POST"),
                        col = c("#1BB6AF", "#EF7C12", "#C70E7B", "#8FDA04", "#4D4D4D", "#FC6882"),
                        splitnode = FALSE,
                        WithConnectors = TRUE,
                        WithNodes = F,
                        soma = NULL,
                        soma.alpha = 1,
                        transmitters = FALSE,
                        transmitter.alpha = 0.3,
                        highflow = FALSE,
                        lwd = 1,
                        radius = NULL,
                        already_selected = NULL,
                        check.template = FALSE,
                        ...){
  someneuronlist = nat::as.neuronlist(someneuronlist)
  prepost = match.arg(prepost)
  if(check.template){
    temps = nat.templatebrains::all_templatebrains()
    temps.microns = c(temps[temps$W<2000,"name"],"JRCFIB2018F")
    reg = nat.templatebrains::regtemplate(someneuronlist)
  }else{
    reg = NULL
  }
  if(is.null(soma)){
    if(is.null(reg)){
      soma = 500
    }else if(!reg%in%temps.microns){
      soma = 500
    }else{
      soma = 4
    }
  }
  if(is.null(radius)){
    if(is.null(reg)){
      radius = 100
    }else if(!reg%in%temps.microns){
      radius = 100
    }else{
      radius = 1
    }
  }
  for (n in 1:length(someneuronlist)){
    neuron = someneuronlist[[n]]
    dendrites.v = subset(rownames(neuron$d), neuron$d$Label == 3)
    axon.v = subset(rownames(neuron$d), neuron$d$Label == 2)
    p.d.v = subset(rownames(neuron$d), neuron$d$Label == 4)
    p.n.v = subset(rownames(neuron$d), neuron$d$Label == 7)
    null.v = subset(rownames(neuron$d), neuron$d$Label == 0 | is.na(neuron$d$Label))
    dendrites = tryCatch( nat::prune_vertices(neuron, verticestoprune = as.integer(c(axon.v, p.d.v, p.n.v, null.v))),
                          error = function(e) NULL)
    axon = tryCatch( nat::prune_vertices(neuron, verticestoprune = as.integer(c(dendrites.v, p.d.v, p.n.v, null.v))),
                     error = function(e) NULL)
    p.d = tryCatch( nat::prune_vertices(neuron, verticestoprune = as.integer(c(axon.v, dendrites.v, p.n.v, null.v))),
                    error = function(e) NULL)
    p.n = tryCatch( nat::prune_vertices(neuron, verticestoprune = as.integer(c(axon.v, dendrites.v, p.d.v, null.v))),
                    error = function(e) NULL)
    tryCatch(rgl::plot3d(dendrites, col = col[1], WithNodes = WithNodes, lwd = lwd, add = TRUE),
             error = function(e) NULL)
    tryCatch(rgl::plot3d(axon, col = col[2], WithNodes = WithNodes, soma = FALSE, lwd = lwd, add = TRUE),
             error = function(e) NULL)
    tryCatch(rgl::plot3d(p.n, col = col[3], WithNodes = WithNodes, soma = FALSE, lwd = lwd, add = TRUE),
    error = function(e) NULL)
    tryCatch(rgl::plot3d(p.d, col = col[4], WithNodes = WithNodes, soma = FALSE, lwd = lwd, add = TRUE),
    error = function(e) NULL)
    tryCatch(rgl::plot3d(neuron, col = "grey30", WithNodes = WithNodes, soma = FALSE, add = TRUE), error = function(e) NULL)
    if(soma){
      rgl::spheres3d(nat::xyzmatrix(neuron)[neuron$StartPoint,], radius = soma, col = col[3], alpha = soma.alpha)
    }
    if(WithConnectors){
      conns=neuron$connectors
      if("BOTH"%in%prepost){
        plot.prepost = c("PRE","POST")
      }else{
        plot.prepost = prepost
      }
      if(transmitters){
        input.synapses = conns[conns$prepost==1,,drop=FALSE]
        output.synapses = conns[conns$prepost==0,,drop=FALSE]
        cols1 = hemibrainr::paper_colours[input.synapses$top_nt]
        cols1[is.na(cols1)] = "black"
        cols2 = hemibrainr::paper_colours[output.synapses$top_nt]
        cols2[is.na(cols2)] = "black"
        if("PRE"%in%plot.prepost){
          tryCatch(rgl::spheres3d(nat::xyzmatrix(output.synapses),
                                  col = cols2, radius = radius, add = TRUE, alpha = transmitter.alpha), #"#EE4244"
                   error = function(e) NULL)
          tryCatch(rgl::points3d(nat::xyzmatrix(output.synapses),
                                 col = hemibrainr::paper_colours["pre"], radius = radius+(radius*0.1), add = TRUE), #"#EE4244"
                   error = function(e) NULL)
        }
        if("POST"%in%plot.prepost){
          tryCatch(rgl::spheres3d(nat::xyzmatrix(input.synapses),
                                  col = cols1, radius = radius/2, add = TRUE, alpha = transmitter.alpha),
                   error = function(e) NULL)
          tryCatch(rgl::points3d(nat::xyzmatrix(input.synapses),
                                 col = hemibrainr::paper_colours["post"], radius = (radius/2)+(radius*0.05), add = TRUE),
                   error = function(e) NULL)
        }
      }else{
        if("POST"%in%plot.prepost){
          tryCatch(rgl::spheres3d(xyzmatrix(conns[conns$prepost==1,,drop=FALSE]),
                                  col = hemibrainr::paper_colours["post"], radius = radius/2, add = TRUE),
                   error = function(e) NULL)
        }
        if("PRE"%in%plot.prepost){
          tryCatch(rgl::spheres3d(xyzmatrix(conns[conns$prepost==0,,drop=FALSE]),
                                  col = hemibrainr::paper_colours["pre"], radius = radius, add = TRUE), #"#EE4244"
                   error = function(e) NULL)
        }
      }
    }
    if(isTRUE(highflow)){
      highest = max(neuron$d[,"flow.cent"])
      s.d = sd(neuron$d[,"flow.cent"], na.rm = T)
      high = subset(neuron$d, neuron$d[,"flow.cent"] > (highest - s.d))
      rgl::points3d(nat::xyzmatrix(high), col = col[6],...)
    }
    if(isTRUE(splitnode)){
      ais = which(apply(neuron$d, 1, function(x) x["flow.cent"] == max(neuron$d[,"flow.cent"])))
      rgl::spheres3d(nat::xyzmatrix(neuron$d[ais,]),radius=radius,col="magenta", add = TRUE, ...)
    }
  }
}


#' @export
#' @rdname plot3d_split
nlscan_split <- function (someneuronlist,
                          already_selected = NULL,
                          prepost = c("BOTH","PRE","POST"),
                          col = c("#1BB6AF", "#EF7C12", "#C70E7B", "#8FDA04", "#4D4D4D", "#FC6882"),
                          splitnode = FALSE,
                          WithConnectors = TRUE,
                          WithNodes = F,
                          soma = NULL,
                          soma.alpha = 1,
                          transmitters = FALSE,
                          transmitter.alpha = 0.3,
                          highflow = FALSE,
                          lwd = 1,
                          radius = NULL,
                          check.template = FALSE,
                          Verbose = TRUE,
                          Wait = TRUE,
                          sleep = 0.1,
                          selected_file = NULL,
                          selected_col = "#fadadd",
                          yaml = TRUE,
                          brain = FALSE, ...){
  prepost = match.arg(prepost)
  if(!requireNamespace('yaml', quietly = TRUE))
    stop("Suggested package yaml is required to use this function!")
  if (nat::is.neuronlist(someneuronlist)) {
    db = someneuronlist
    neurons = names(someneuronlist)
    nams = as.data.frame(someneuronlist, stringsAsFactors = FALSE)$name
  }
  frames <- length(neurons)
  if(!is.issue(already_selected)){
    selected <- as.character(intersect(already_selected, names(someneuronlist)))
  }else{
    selected <- character()
  }
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
  if (!isFALSE(brain)){
    brain = brain
  }
  while (TRUE) {
    if (i > length(neurons) || i < 1){
      rgl::clear3d()
      rgl::bg3d(color = "white")
      if(length(selected)){
        if (!isFALSE(brain)){
          rgl::plot3d(brain, col = "grey70", alpha = 0.1)
        }
        rgl::plot3d(someneuronlist[selected], ..., col = hemibrain_bright_colour_ramp(length(selected)))
      }
      end = hemibrain_choice("Done selecting/scanning neurons? yes/no ")
      rgl::clear3d()
      if(end){
        break
      }else{
       i <- 1
      }
    }
    n <- neurons[i]
    cat("Current neuron:", n, " ", nams[i]," (", i, "/", length(neurons),")",ifelse(n %in% selected,"SELECTED","NOT SELECTED"), "\n")
    si = someneuronlist[[i]]$AD.segregation.index
    if(is.null(si)){
      si <- someneuronlist[i,"segregation.index"]
    }
    cat("segregation index: ", si)
    if(n %in% selected){
      rgl::bg3d(color = selected_col)
    }else{
      rgl::bg3d(color = "white")
    }
    if (!isFALSE(brain)){
      rgl::plot3d(brain, col = "grey70", alpha = 0.1)
    }
    pl <- plot3d_split(someneuronlist[i],
                       col = col,
                       WithConnectors = WithConnectors,
                       WithNodes = WithNodes,
                       highflow = highflow,
                       radius = radius,
                       prepost = prepost,
                       splitnode = splitnode,
                       soma = soma,
                       soma.alpha = soma.alpha,
                       transmitters = transmitters,
                       transmitter.alpha = transmitter.alpha,
                       lwd = lwd,
                       check.template = check.template,
                       ...)
    if (Wait) {
      chc <- readline("Return to continue, b to go back, s to (de)select, d [save to disk], t to stop, c to cancel (without returning a selection): ")
      if (chc == "c" || chc == "t") {
        sapply(pl, rgl::rgl.pop, type = "shape")
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
    rgl::clear3d()
    rgl::bg3d(color = "white")
  }
  if (is.null(chc) || chc == "c"){
    rgl::bg3d(color = "white")
    return(NULL)
  }
  if (!is.null(selected_file)){
    savetodisk(selected, selected_file)
  }
  selected
}

#' @export
#' @rdname plot3d_split
plot3d_somas <- function(someneuronlist,
                         col = hemibrain_bright_colors["green"],
                         radius = 300,
                         soma.alpha = 1,
                         volume = NULL,
                         invert = FALSE){
  somapos <- as.data.frame(catmaid::soma(someneuronlist), stringsAsFactors = FALSE)
  somapos <- somapos[!is.na(somapos$X),]
  if(!is.null(volume)){
    p = nat::pointsinside(somapos, surf = volume)
    if(invert){
      somapos = somapos[p,]
    }else{
      somapos = somapos[!p,]
    }
  }
  rgl::spheres3d(somapos, radius = radius, alpha = soma.alpha, col = col)
}


#' Plot chosen neurons and see all neurons in that type
#'
#' @description Plot your selected neuron (rightmost pane) and then see all of the neurons in the hemibrain project that belong to
#' the types you have plotted (left pane).
#'
#' @param someneuronlist a neuronlist or neuron object to plot. If \code{bodyids} is \code{NULL} then \code{bodyids = names(someneuronlist)}.
#' @param bodyids bodyids for the neurons you want to plot. If \code{someneuronlist} is null then \code{hemibrain_neurons} will be used
#' to load a neuron database.
#' @param print logical, whether or not to print the name of the cell_types being plotted and their Body IDs.
#' @param meta a \code{data.frame} with a \code{bodyid} and a \code{type} field, such as one pulled by \code{neuprint_get_meta}. If \code{NULL}
#' then meta will be acquired as \code{neuprint_get_meta(bodyids)}.
#' @param brain a \code{mesh3d} or \code{hxsurf} object, such as a template brain, to plot alongside neurons of interest.
#' @param ... vector lists of bodyids, each of which will be plotted in a different set of colours.
#'
#' @return Plots coloured neuron(s)
#' @examples
#' \donttest{
#' \dontrun{
#'
#' # Choose neurons
#' ids = c("605063393", "573708714", "485387561", "5813083717")
#'
#' # Plot
#' hemibrain_type_plot(ids)
#'
#' # Plot multi
#' open3d()
#' a = c("546511447", "573333903", "577546843", "573333284", "573683455",
#'  "609941228", "604368022")
#' b = c("573670049", "579576294", "727838223", "696432315")
#' hemibrain_multi3d(a, b)
#'
#' }}
#' @export
#' @seealso \code{\link{plot3d_split}}
#' @importFrom grDevices rainbow
hemibrain_type_plot <- function(bodyids = NULL,
                                someneuronlist = NULL,
                                brain = hemibrainr::hemibrain.surf,
                                print = TRUE,
                                meta = NULL){

  # Get bodyids
  if(is.null(bodyids)){
    if(is.null(someneuronlist)){
      stop("Either bodyids or someneuronlist must be given")
    }
    bodyids = names(someneuronlist)
  }

  # Get metadata
  if(is.null(meta)){
    meta = neuprintr::neuprint_get_meta(bodyids)
  }else{
    types = subset(meta, meta$bodyid %in% bodyids)$type
    meta = subset(meta, meta$type %in% unique(types))
    missing = setdiff(bodyids,meta$bodyid)
    if(length(missing)){
      meta = plyr::rbind.fill(meta,
                              neuprintr::neuprint_get_meta(missing))
    }
  }
  # NB this means that NA -> "NA"
  ito.cts = as.character(unique(meta$type))

  # Neurons
  bodyids = as.character(bodyids)
  all.bodyids = unique(c(bodyids,as.character(unique(meta$bodyid))))
  if(is.null(someneuronlist)){
    someneuronlist = hemibrain_neurons()
  }
  neurons = tryCatch(someneuronlist[all.bodyids], error = function(e) NULL)
  if(is.null(neurons)){
      neurons = neuprintr::neuprint_read_neurons(all.bodyids)
  }

  # Set up rgl
  rgl::clear3d()
  mat <- matrix(1:2, nrow = 1, 2)
  rgl::layout3d(mat, height = 3, sharedMouse = TRUE)

  # Plot
  rgl::rgl.viewpoint(userMatrix = structure(c(0.988000273704529, -0.0847902745008469,
                                         -0.129098162055016, 0, 0.0781147107481956, -0.446752369403839,
                                         0.891240775585175, 0, -0.133243337273598, -0.890630483627319,
                                         -0.434768080711365, 0, 0, 0, 0, 1), .Dim = c(4L, 4L)), zoom = 1)
  rgl::plot3d(hemibrainr::hemibrain.surf, col = "grey", alpha = 0.1, add = TRUE)
  cols = rainbow(length(ito.cts))
  for(i in 1:length(ito.cts)){
    ids = subset(meta,meta$type==ito.cts[i])$bodyid
    n = neurons[as.character(ids)]
    col = grDevices::colorRampPalette(colors = c(cols[i],"grey10"))
    col = col(length(n)+2)[1:length(n)]
    rgl::plot3d(n, lwd = 2, soma = TRUE, col = col)
    if(print){
      message("Type: ", ito.cts[i])
      print(dput(names(n)))
    }
  }

  # Plot the given IDs
  rgl::next3d()
  rgl::rgl.viewpoint(userMatrix = structure(c(0.988000273704529, -0.0847902745008469,
                                         -0.129098162055016, 0, 0.0781147107481956, -0.446752369403839,
                                         0.891240775585175, 0, -0.133243337273598, -0.890630483627319,
                                         -0.434768080711365, 0, 0, 0, 0, 1), .Dim = c(4L, 4L)), zoom = 1)
  rgl::plot3d(brain, col = "grey", alpha = 0.1, add = TRUE)
  rgl::plot3d(neurons[bodyids], col = hemibrain_bright_colour_ramp(length(bodyids)), lwd=2, soma = TRUE)

  # Body Ids
  if(print){
    message("right: ")
    print(dput(bodyids))
  }
}


#' @export
#' @rdname hemibrain_type_plot
hemibrain_multi3d <- function(..., someneuronlist = hemibrain_neurons()){
  m = as.list(match.call())
  count = length(m)-1
  cols = rainbow(count)
  for(i in 1:count){
    j = i+1
    n = as.character(get(as.character(m[[j]])))
    n = n[n%in%names(someneuronlist)]
    col = grDevices::colorRampPalette(colors = c(cols[i],"grey10"))
    col = col(length(n)+2)[1:length(n)]
    rgl::plot3d(someneuronlist[n], lwd = 2, col = col, soma = TRUE)
  }
}

