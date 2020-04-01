#######################################################################################
################################ axon-dendrite split ##################################
#######################################################################################

#' Determine a neuron's dendrite and axon by calculating flow centrality
#'
#' @description implementation of the algorithm for calculating flow
#'   centralities from Schneider-Mizell et al. (2016) (slightly modified).
#'   Can be used on \code{nat::neuronlist}
#'   objects read using \code{neuprintr::neuprint_read_neurons}. The example code
#'   below gives the recommended arguments when using hemibrain data.
#'
#' @param x a \code{nat::neuronlist} or \code{nat::neuron} object. It is assumed
#' that this neuron has been read in by \code{neuprintr::neuprint_read_neurons} or
#' possibly \code{catmaid::read.neurons.catmaid}.
#' @param mode type of flow centrality to calculate. There are three flavours:
#'   (1) centrifugal, which counts paths from proximal inputs to distal outputs;
#'   (2) centripetal, which counts paths from distal inputs to proximal outputs;
#'   and (3) the sum of both.
#' @param polypre whether to consider the number of presynapses as a multiple of
#'   the numbers of connections each makes
#' @param soma logical, whether or not the given neuron has a soma.
#' The soma should be its root. If it does, \code{\link{primary_neurite}}
#' will be called to make sure the primary neurite is labelled correctly.
#' @param primary.dendrite whether to try to assign nodes to a 'primary
#'   dendrite'. Defaults to considering nodes of 0.9*maximal flow centrality.
#'   Assigning to NULL will prevent generating this compartment.
#' @param primary.branchpoint the proportion of the maximum branchpoint score needed
#' for a point to be assigned as the primary branchpoint. Used only when \code{soma = TRUE}.
#' See the \code{first} argument in \code{\link{primary_neurite}} for details.
#' @param split the algorithm will assign two main neurite compartments,
#' which as per SWC format will be indicates as either axon (Label =2)
#' or dendrite (Label = 3) in the returned objects, at neuron$d$Label.
#' This assignment can be based which compartment contains the most
#' postsynapses ("postsynapses") or presynapses ("presynapses"),
#' or the geodesic distance of its first branch point from the primary branch
#' point (i.e. the first branch point from the soma) ("distance").
#' The compartment closest to this branchpoint is most usually the dendrite
#' (Bates & Schlegel 2020).
#' @param ... methods sent to \code{nat::nlapply}
#'
#' @details From Schneider-Mizell et al. (2016): "We use flow centrality for
#'   four purposes. First, to split an arbor into axon and dendrite at the
#'   maximum centrifugal SFC, which is a preliminary step for computing the
#'   segregation index, for expressing all kinds of connectivity edges (e.g.
#'   axo-axonic, dendro-dendritic) in the wiring diagram, or for rendering the
#'   arbor in 3d with differently colored regions. Second, to quantitatively
#'   estimate the cable distance between the axon terminals and dendritic arbor
#'   by measuring the amount of cable with the maximum centrifugal SFC value.
#'   Third, to measure the cable length of the main dendritic shafts using
#'   centripetal SFC, which applies only to insect neurons with at least one
#'   output syn- apse in their dendritic arbor. And fourth, to weigh the color
#'   of each skeleton node in a 3d view, providing a characteristic signature of
#'   the arbor that enables subjective evaluation of its identity."
#'
#' @references Schneider-Mizell, C. M., Gerhard, S., Longair, M., Kazimiers, T.,
#'   Li, F., Zwart, M. F., … Cardona, A. (2015). Quantitative neuroanatomy for
#'   connectomics in Drosophila. bioRxiv, 026617. http://doi.org/10.1101/026617
#'
#' @return the neuron or neuron list object inputted, with centripetal flow
#'   centrality information added to neuron$d and a segregation index score.
#'   The neuron$d$Label now gives the compartment, where axon is Label = 2,
#'   dendrite Label = 3, primary dendrite Label = 9 and
#'   primary neurite Label = 7.
#'   Soma is Label = 1.
#'
#' @examples
#' \donttest{
#'
#' # Choose neurons
#' ## These neurons are some 'tough' examples from the hemibrain:v1.0.1
#' ### They will split differently depending on the parameters you use.
#' tough = c("5813056323", "579912201", "5813015982", "973765182", "885788485",
#' "915451074", "5813032740", "1006854683", "5813013913", "5813020138",
#' "853726809", "916828438", "5813078494", "420956527", "486116439",
#' "573329873", "5813010494", "5813040095", "514396940", "665747387",
#' "793702856", "451644891", "482002701", "391631218", "390948259",
#' "390948580", "452677169", "511262901", "422311625", "451987038"
#' )
#' # for documentation purposes only run first 5 examples
#' tough=tough[1:5]
#'
#' # Get neurons
#' neurons = neuprintr::neuprint_read_neurons(tough)
#'
#' # Now make sure the neurons have a soma marked
#' ## Some hemibrain neurons do not, as the soma was chopped off
#' neurons.checked = hemibrain_skeleton_check(neurons, meshes = hemibrain.surf)
#'
#' # Split neuron
#' ## These are the recommended parameters for hemibrain neurons
#' neurons.flow = flow_centrality(neurons.checked,
#' polypre = TRUE,
#' mode = "centrifugal",
#' split = "distance")
#'
#' \dontrun{
#' # Plot the split to check it
#' nat::nopen3d()
#' nlscan_split(neurons.flow, WithConnectors = TRUE)
#'
#' }}
#' @export
#' @seealso \code{\link{primary_neurite}}, \code{\link{hemibrain_skeleton_check}}, \code{\link{hemibrain_flow_centrality}}, \code{\link{hemibrain_use_splitpoints}}
#' \code{\link{hemibrain_splitpoints}},
flow_centrality <-function(x,
                           mode = c("centrifugal", "centripetal", "sum"),
                           polypre = TRUE,
                           soma = TRUE,
                           primary.dendrite = 0.9,
                           primary.branchpoint = 0.25,
                           split = c("distance","synapses"),
                           ...) UseMethod("flow_centrality")

#' @export
flow_centrality.neuron <- function(x,
                                   mode = c("centrifugal", "centripetal", "sum"),
                                   polypre = TRUE,
                                   soma = TRUE,
                                   primary.dendrite = 0.9,
                                   primary.branchpoint = 0.25,
                                   split = c("distance", "synapses"),
                                   ...){
  split = match.arg(split)
  mode = match.arg(mode)
  x$d$Label = 0
  el = x$d[x$d$Parent != -1, c("Parent", "PointNo")]
  n = nat::ngraph(data.matrix(el[, 2:1]), x$d$PointNo, directed = TRUE, xyz = nat::xyzmatrix(x$d), diam = x$d$W)
  igraph::V(n)$name = igraph::V(n)
  leaves = which(igraph::degree(n, v = igraph::V(n), mode = "in") == 0, useNames = T)
  root = which(igraph::degree(n, v = igraph::V(n), mode = "out") == 0, useNames = T)
  segs = x$SegList
  nodes = x$d
  nodes[, c("post","pre","up.syns.in","up.syns.out","flow.cent")] = 0
  nodes[,"Label"] = 3
  nodes = nodes[unlist(c(root, lapply(segs, function(x) x[-1]))),]
  nodes = nodes[order(as.numeric(rownames(nodes))), ]
  syns.in = x$connectors[x$connectors$prepost == 1, ][, "treenode_id"]
  if (polypre) {
    syns.out = x$connectors[x$connectors$prepost == 0,][, "treenode_id"]
  }else {
    syns.out = x$connectors[x$connectors$prepost ==0 & !duplicated(x$connectors$connector_id), ][, "treenode_id"]
  }
  point.no.in = rownames(nodes)[match(syns.in, nodes[, "PointNo"])]
  nodes.in = rep(1, length(point.no.in))
  names(nodes.in) = point.no.in
  nodes.in = tapply(nodes.in, point.no.in, sum)
  point.no.out = rownames(nodes)[match(syns.out, nodes[, "PointNo"])]
  nodes.out = rep(1, length(point.no.out))
  names(nodes.out) = point.no.out
  nodes.out = tapply(nodes.out, point.no.out, sum)
  nodes[names(nodes.in), "post"] <- nodes.in
  nodes[names(nodes.out), "pre"] <- nodes.out
  ins = c(0, lapply(segs, function(x) if (!is.null(x)) {
    rev(cumsum(rev(unlist(lapply(x, function(x) ifelse(x %in% names(nodes.in), nodes[as.character(x), "post"],  0))))))
  }))
  outs = c(0, lapply(segs, function(x) if (!is.null(x)) {
    rev(cumsum(rev(unlist(lapply(x, function(x) ifelse(x %in%
                                                         names(nodes.out), nodes[as.character(x), "pre"],
                                                       0))))))
  }))
  nodes[, "up.syns.in"] = nodes[, "up.syns.in"] + c(0, unlist(lapply(ins, function(x) x[-1])))
  nodes[, "up.syns.out"] = nodes[, "up.syns.out"] + c(0, unlist(lapply(outs,  function(x) x[-1])))
  in.bps = unlist(lapply(ins, function(x) x[1]))[-1]
  out.bps = unlist(lapply(outs, function(x) x[1]))[-1]
  names(in.bps) = names(out.bps) = bps = c(root, unlist(lapply(segs, function(x) x[1]))[-1])
  in.bps.child = tapply(in.bps, names(in.bps), function(x) ifelse(names(x)[1] %in%
                                                                    names(nodes.in), sum(x) - (nodes.in[names(x)[1]] * length(x)),
                                                                  sum(x)))
  out.bps.child = tapply(out.bps, names(out.bps), function(x) ifelse(names(x)[1] %in%
                                                                       names(nodes.out), sum(x) - (nodes.out[names(x)[1]] *
                                                                                                     length(x)), sum(x)))
  nodes[names(in.bps.child), "up.syns.in"] = in.bps.child
  nodes[names(out.bps.child), "up.syns.out"] = out.bps.child
  bps = as.numeric(names(in.bps.child))
  for (i in 1:length(bps)) {
    bp = bps[i]
    new.in = new.out = c(rep(0, nrow(nodes)))
    names(new.in) = names(new.out) = rownames(nodes)
    vertices = unlist(igraph::shortest_paths(n, bp, to = root)$vpath)[-1]
    new.in[as.character(vertices)] = new.in[as.character(vertices)] + in.bps.child[i]
    new.out[as.character(vertices)] = new.out[as.character(vertices)] + out.bps.child[i]
    nodes[, "up.syns.in"] = nodes[, "up.syns.in"] + new.in
    nodes[, "up.syns.out"] = nodes[, "up.syns.out"] + new.out
  }
  in.total = sum(nodes[,"post"])
  out.total = sum(nodes[,"pre"])
  if(mode == "centrifugal") {
    nodes[, "flow.cent"] = (in.total - nodes[, "up.syns.in"]) * nodes[, "up.syns.out"]
  }else if (mode == "centripetal") {
    nodes[, "flow.cent"] = (out.total - nodes[, "up.syns.out"]) * nodes[, "up.syns.in"]
  }else {
    nodes[, "flow.cent"] = ((in.total - nodes[, "up.syns.in"]) * nodes[, "up.syns.out"]) + ((out.total - nodes[, "up.syns.out"]) * nodes[, "up.syns.in"])
  }
  for (bp in bps) {  # bending flow
    down = unlist(igraph::ego(n, 1, nodes = bp, mode = "in", mindist = 0))[-1]
    bf = c()
    for (u in down) {
      this.seg.posts = nodes[u, ]$up.syns.in
      other.segs.pre = sum(nodes[down[!down == u], ]$up.syns.out)
      bf = c(bf, sum(this.seg.posts * other.segs.pre))
    }
    nodes[bp, "flow.cent"] = nodes[bp, "flow.cent"] + sum(bf)
  }
  high = max(nodes[!rownames(nodes)%in%bps, "flow.cent"]) # some bps are too high!
  if (!is.null(primary.dendrite)) {
    highs = subset(rownames(nodes), nodes[, "flow.cent"] >= (primary.dendrite * high))
  }else {
    primary.dendrite = 0.9
    highs = subset(rownames(nodes), nodes[, "flow.cent"] >= (primary.dendrite * high))
  }
  highs = as.numeric(unique(unlist(igraph::shortest_paths(n,highs, to = highs, mode = "all")))) # fill in any missed points
  if(soma){
    primary.branch.point = primary_branchpoint(x, primary_neurite = TRUE, first = primary.branchpoint)
  }else{
    primary.branch.point = primary_branchpoint(x, primary_neurite = FALSE)
  }
  p.n = suppressWarnings(unique(unlist(igraph::shortest_paths(n, as.numeric(root), to = as.numeric(primary.branch.point), mode = "all")$vpath)))
  primary.branch.points = p.n[p.n%in%nat::branchpoints(nodes)]
  primary.branch.point.downstream = suppressWarnings(unique(unlist(igraph::shortest_paths(n, as.numeric(primary.branch.point), to = as.numeric(leaves), mode = "in")$vpath)))
  p.n = p.n[1:ifelse(sum(p.n%in%highs)>1,min(which(p.n%in%highs)),length(p.n))]
  select.highs = highs[!highs%in%c(root)]
  ais = select.highs[nodes[select.highs,"flow.cent"] >= high ] # Point of highest flow
  if (length(ais) > 0) {
    runstosoma = unlist(lapply(ais, function(x) length(unlist(igraph::shortest_paths(n,x, to = root)$vpath))))
    ais = as.numeric(ais[match(min(runstosoma), runstosoma)])
  }
  ### Put the main branches into two groups, which will become axon and dendrite
  if(ais %in% c(primary.branch.points)){
    down = unlist(igraph::ego(n, 1, nodes = ais, mode = "in", mindist = 0))[-1]
    downstreams = list()
    order = nrow(x$d)
    for(i in 1:length(down)){
      downstream = unlist(igraph::ego(n, order = order, nodes = down[i], mode = c("in"), mindist = 0))
      downstreams[[i]] = downstream
    }
    downstream.lengths = sapply(downstreams,length)
    two.main.branches = downstreams[(1:length(downstreams))[order(downstream.lengths,decreasing = TRUE)][1:2]]
    upstream = two.main.branches[[1]]
    downstream = two.main.branches[[2]]
  }else{
    downstream = suppressWarnings(unique(unlist(igraph::shortest_paths(n, ais, to = leaves, mode = "in")$vpath)))
    downstream.unclassed = downstream[!downstream %in% c(ais)]
    remove = rownames(nodes)[!rownames(nodes) %in% downstream.unclassed]
    downstream.g = igraph::delete_vertices(n, v = as.character(remove))
    main = igraph::components(downstream.g)
    club = main$membership
    not.in.club = rownames(nodes)[!rownames(nodes) %in% names(main$membership)]
    going.into.club = rep(max(main$membership)+1,length(not.in.club))
    names(going.into.club) = not.in.club
    club = c(club,going.into.club)
    directions = c()
    for(component in unique(club)){
      branch = names(club[club == component])
      branch.out = sum(nodes[branch,"pre"])
      branch.in = sum(nodes[branch,"post"])
      flow.into.branch = branch.out*(in.total-branch.in)
      flow.outof.branch = branch.in*(out.total-branch.out)
      direction = flow.into.branch>flow.outof.branch
      directions = c(directions, direction)
    }
    if(length(which(directions))==0){
      upstream = rownames(nodes)[!rownames(nodes) %in% downstream]
    }else{
      hits = which(directions)
      downstream = names(club[club %in% hits])
      upstream = names(club[!club %in% hits])
    }
  }
  ### Get main axon/dendrite start points
  downstream.unclassed = downstream[!downstream %in% c(p.n, highs, root, primary.branch.point)]
  remove = rownames(nodes)[!rownames(nodes) %in% downstream.unclassed]
  downstream.g = igraph::delete_vertices(n, v = as.character(remove))
  main1 = igraph::components(downstream.g)
  main1 = names(main1$membership[main1$membership %in% which.max(table(main1$membership))])
  nodes.downstream = nodes[as.character(main1), ]
  downstream.tract.parent = unique(nodes.downstream$Parent[!nodes.downstream$Parent %in% nodes.downstream$PointNo])
  downstream.tract.parent = rownames(nodes)[match(downstream.tract.parent,nodes$PointNo)]
  if (sum(downstream.tract.parent %in% c(p.n,highs,primary.branch.point)) > 0) {
    bps.all = intersect(nat::branchpoints(nodes),main1)
    bps.downstream = bps.all[bps.all %in% downstream.unclassed]
    runstoprimarybranchpoint = unlist(lapply(bps.downstream,
                                             function(x) length(unlist(suppressWarnings(igraph::shortest_paths(n,to = as.numeric(downstream.tract.parent), from = x)$vpath)))))
    downstream.tract.parent = bps.downstream[which.min(runstoprimarybranchpoint)]
  }
  upstream.unclassed = upstream[!upstream %in% c(p.n, highs, root, primary.branch.point)]
  remove = rownames(nodes)[!rownames(nodes) %in% intersect(upstream.unclassed,primary.branch.point.downstream)]
  upstream.g = igraph::delete_vertices(n, v = as.character(remove))
  main2 = igraph::components(upstream.g)
  main2 = names(main2$membership[main2$membership %in% which.max(table(main2$membership))])
  nodes.upstream = nodes[as.character(main2), ]
  upstream.tract.parent = unique(nodes.upstream$Parent[!nodes.upstream$Parent %in%nodes.upstream$PointNo])
  upstream.tract.parent = ifelse(upstream.tract.parent==-1,root,upstream.tract.parent)
  upstream.tract.parent = rownames(nodes)[match(upstream.tract.parent,nodes$PointNo)]
  if (sum(upstream.tract.parent %in% c(p.n,highs,primary.branch.point)) > 0) {
    bps.all = intersect(nat::branchpoints(nodes),main2)
    bps.upstream = bps.all[bps.all %in% upstream.unclassed]
    runstoprimarybranchpoint = unlist(lapply(bps.upstream,function(x) length(unlist(suppressWarnings(igraph::shortest_paths(n,to = as.numeric(upstream.tract.parent), from = x)$vpath)))))
    upstream.tract.parent = bps.upstream[which.min(runstoprimarybranchpoint)]
  }
  ### Assign axon versus dendrite
  if (grepl("synapses", split)) {
    ## flow in and out of 'upstream' complex
    upstream.out = sum(nodes[upstream.unclassed,"pre"])
    upstream.in = sum(nodes[upstream.unclassed,"post"])
    flow.outof.upstream = upstream.out*(in.total-upstream.in)
    flow.into.upstream = upstream.in*(out.total-upstream.out)
    upstream.flow = flow.outof.upstream-flow.into.upstream
    ## flow in and out of 'downstream' complex
    downstream.out = sum(nodes[downstream.unclassed,"pre"])
    downstream.in = sum(nodes[downstream.unclassed,"post"])
    flow.outof.downstream = downstream.out*(in.total-downstream.in)
    flow.into.downstream = downstream.in*(out.total-downstream.out)
    downstream.flow = flow.outof.downstream-flow.into.downstream
    ## Compare net flow
    choice = downstream.flow > upstream.flow
    if (downstream.flow == upstream.flow) {
      split = "distance"
      warning("synapse numbers are the same,
              splitting based on branch point distances to primary branchpoint")
    }
    else if (choice) {
      nodes[as.character(downstream.unclassed), "Label"] = 2
      axon.nodes = downstream.unclassed
      dendrite.nodes = upstream.unclassed
    }
    else {
      nodes[as.character(upstream.unclassed), "Label"] = 2
      axon.nodes = upstream.unclassed
      dendrite.nodes = downstream.unclassed
    }
  }
  if (split == "distance") {
    dist.upstream.to.primary.branchpoint = suppressWarnings(length(unlist(igraph::shortest_paths(n,to = as.numeric(primary.branch.point), from = as.numeric(upstream.tract.parent))$vpath)))
    dist.downstream.to.primary.branchpoint = suppressWarnings(length(unlist(igraph::shortest_paths(n,to = as.numeric(primary.branch.point), from = as.numeric(downstream.tract.parent))$vpath)))
    if (dist.upstream.to.primary.branchpoint < dist.downstream.to.primary.branchpoint) {
      nodes[as.character(downstream.unclassed), "Label"] = 2
      axon.nodes = downstream.unclassed
      dendrite.nodes = upstream.unclassed
    }
    else if (dist.upstream.to.primary.branchpoint > dist.downstream.to.primary.branchpoint) {
      nodes[as.character(upstream.unclassed), "Label"] = 2
      axon.nodes = upstream.unclassed
      dendrite.nodes = downstream.unclassed
    }
    else {
      warning("branch point distances are the same, splitting based on postsynapses")
      choice = sum(nodes[as.character(downstream.unclassed),
                         "post"]) < sum(nodes[as.character(upstream.unclassed),
                                              "post"])
      if (choice) {
        nodes[as.character(downstream.unclassed), "Label"] = 2
        axon.nodes = downstream.unclassed
        dendrite.nodes = upstream.unclassed
      }
      else {
        nodes[as.character(upstream.unclassed), "Label"] = 2
        axon.nodes = upstream.unclassed
        dendrite.nodes = downstream.unclassed
      }
    }
  }
  p.d = as.numeric(unique(unlist(igraph::shortest_paths(n, primary.branch.point, to = highs, mode = "all"))))
  nodes[p.d, "Label"] = 4
  nodes[p.n, "Label"] = 7
  if(length(p.d)==nrow(nodes)){
    nodes[p.d, "Label"] = 3
    p.d = NULL
  }
  if(soma){
    if(is.null(p.d)|length(dendrite.nodes)==nrow(nodes)|length(axon.nodes)==nrow(nodes)){
      p.n = primary_neurite(x, neuron = FALSE)
      primary.branch.point = p.n[length(p.n)]
      nodes[p.n, "Label"] = 7
    }
  }
  dendrites = subset(nodes, nodes$Label == 3)
  dendrites.post = sum(subset(dendrites$post, dendrites$post >0))
  dendrites.pre = sum(subset(dendrites$pre, dendrites$pre >0))
  dendrites.both = dendrites.post + dendrites.pre
  dendrites.pi = dendrites.post/dendrites.both
  dendrites.si = -(dendrites.pi * log(dendrites.pi) + (1 -dendrites.pi) * log(1 - dendrites.pi))
  if (is.nan(dendrites.si)) {
    dendrites.si = 0
  }
  axon = subset(nodes, nodes$Label == 2)
  axon.post = sum(subset(axon$post, axon$post > 0))
  axon.pre = sum(subset(axon$pre, axon$pre > 0))
  axon.both = axon.post + axon.pre
  axon.pi = axon.post/axon.both
  axon.si = -(axon.pi * log(axon.pi) + (1 - axon.pi) * log(1 -axon.pi))
  if (is.nan(axon.si)) {
    axon.si = 0
  }
  entropy.score = (1/(dendrites.both + axon.both)) * ((axon.si *axon.both) + (dendrites.si * dendrites.both))
  both.comps = (dendrites.post + axon.post)/(dendrites.both +axon.both)
  control.score = -(both.comps * log(both.comps) + (1 - both.comps) *log(1 - both.comps))
  segregation.index = 1 - (entropy.score/control.score)
  if (is.na(segregation.index)) {
    segregation.index = 0
  }
  secondary.branch.points = as.numeric(c(downstream.tract.parent, upstream.tract.parent))
  starts = sapply(unique(c(p.n,p.d)),function(s) igraph::neighbors(n, v=s, mode = c("in")))
  starts = unique(unlist(starts))
  axon.starts = as.numeric(rownames(axon[rownames(axon)%in%starts,]))
  dendrite.starts = as.numeric(rownames(dendrites[rownames(dendrites)%in%starts,]))
  x$d = nodes
  x$AD.segregation.index = segregation.index
  x$primary.branch.point = as.numeric(primary.branch.point)
  x$axon.start = nullToNA(axon.starts)
  x$dendrite.start = nullToNA(dendrite.starts)
  nsbp <- nodes[secondary.branch.points,]
  x$axon.primary = nullToNA(nsbp$PointNo[!nsbp$PointNo%in%dendrite.nodes])
  x$dendrite.primary = nullToNA(nsbp$PointNo[nsbp$PointNo%in%dendrite.nodes])
  x$secondary.branch.points = secondary.branch.points
  x$max.flow.centrality = as.numeric(ais)
  x$split = TRUE
  x = hemibrain_neuron_class(x)
  x
}

#' @export
flow_centrality.neuronlist <- function(x,
                                       mode = c("centrifugal", "centripetal", "sum"),
                                       polypre = TRUE,
                                       soma = TRUE,
                                       primary.dendrite = 0.9,
                                       primary.branchpoint = 0.25,
                                       split = c("distance","synapses"),
                                       ...){
  split = match.arg(split)
  mode = match.arg(mode)
  neurons = nat::nlapply(x, flow_centrality.neuron, mode = mode,
                         polypre = polypre, soma = soma,
                         primary.dendrite = primary.dendrite,
                         split = split,
                         ...)
  neurons
}


#' Get the positions of key points in a 'split' neuron
#'
#' @description Get the positions of key points, i.e. primary and secondary branchpoints and
#' a neuron's root, in a 'split' neuron.
#'
#' @inheritParams flow_centrality
#'
#' @return a \code{data.frame}
#' @seealso \code{\link{flow_centrality}}, \code{\link{hemibrain_use_splitpoints}}
#' @export
hemibrain_splitpoints <- function(x){
  x = nat::as.neuronlist(x)
  splits = data.frame()
  for(i in 1:length(x)){
    bi = names(x)[i]
    n = x[[i]]
    if(is.null(bi)){
      bi = nullToNA(n$bodyid)
    }
    points = c(root = nullToNA(nat::rootpoints(n)),
              primary.branch.point = nullToNA(n$primary.branch.point),
              axon.primary = nullToNA(n$axon.primary),
              dendrite.primary = nullToNA(n$dendrite.primary),
              axon.start = nullToNA(n$axon.start),
              dendrite.start = nullToNA(n$dendrite.start),
              linker = nullToNA(n$linker))
    points.m = reshape2::melt(points)
    points.m$point = rownames(points.m)
    rownames(points.m) = NULL
    xyz = nat::xyzmatrix(n$d[points,])
    pos = suppressWarnings( do.call(cbind, list(bodyid = bi, points.m, xyz)) )
    splits = rbind(splits,pos)
  }
  colnames(splits) = c("bodyid", "position", "point", "X", "Y", "Z")
  splits
}


#' Split a neuron using stored 'splitpoints'
#'
#' @description Split a neuron using a \code{data.frame} of stored split points.
#'
#' @inheritParams flow_centrality
#' @param df a \code{data.frame} of splitpoints from running \code{\link{flow_centrality}},
#' as produced by \code{\link{hemibrain_splitpoints}}.
#' @param knn logical, whether or not to find corresponding points
#' between \code{df} and \code{x$d} using a nearest neighbour search in 3D
#' space (TRUE) or just use the given point IDs in \code{df} (i.e. if neurons
#' have not been resampled or their skeletons otherwise modified).
#'
#' @inherit flow_centrality return
#'
#' @return a \code{neuronlist}
#' @seealso \code{\link{flow_centrality}}
#' @examples
#' \donttest{
#'
#' # Choose neurons
#' ## These neurons are some 'tough' examples from the hemibrain:v1.0.1
#' ### They will split differently depending on the parameters you use.
#' tough = c("5813056323", "579912201", "5813015982", "973765182", "885788485",
#' "915451074", "5813032740", "1006854683", "5813013913", "5813020138",
#' "853726809", "916828438", "5813078494", "420956527", "486116439",
#' "573329873", "5813010494", "5813040095", "514396940", "665747387",
#' "793702856", "451644891", "482002701", "391631218", "390948259",
#' "390948580", "452677169", "511262901", "422311625", "451987038"
#' )
#' # for documentation purposes only run first 5 examples
#' tough=tough[1:5]
#'
#' # Get neurons
#' neurons = neuprintr::neuprint_read_neurons(tough)
#'
#' # Now make sure the neurons have a soma marked
#' ## Some hemibrain neurons do not, as the soma was chopped off
#' neurons.checked = hemibrain_skeleton_check(neurons, meshes = hemibrain.surf)
#'
#' # Split neuron
#' ## These are the recommended parameters for hemibrain neurons
#' neurons.flow = flow_centrality(neurons.checked, polypre = TRUE,
#' mode = "centrifugal",
#' split = "distance")
#'
#' # Save the results
#' splitpoints = hemibrain_splitpoints(neurons.flow)
#'
#' # Re-use the results
#' neurons.flow.2 = hemibrain_use_splitpoints(neurons, splitpoints, knn = FALSE)
#'
#' \dontrun{
#' # Plot the split to check it
#' nat::nopen3d()
#' nlscan_split(neurons.flow2, WithConnectors = TRUE)
#' }}
#' @export
#' @seealso \code{\link{hemibrain_splitpoints}}, \code{\link{flow_centrality}}, \code{\link{hemibrain_precomputed_splitpoints}}
hemibrain_use_splitpoints <-function(x,
                                     df,
                                     knn = FALSE,
                                     ...) UseMethod("hemibrain_use_splitpoints")

#' @export
hemibrain_use_splitpoints.neuron <-function(x, df, knn = FALSE, ...){

  # Get splitpoints
  df = df[df$bodyid == x$bodyid,]
  df = df[,c("bodyid", "position", "point", "X", "Y", "Z")]
  df[df==""] = NA

  if(!nrow(df)){
    warning("bodyid ", x$bodyid," not in splitpoints, returning neuron unmodified")
    y = x
    y$split = FALSE
  }else{
    # Find nearest points on skeleton
    if(knn){
      near = nabor::knn(query = nat::xyzmatrix(df), data = nat::xyzmatrix(x), k  = 1)
      df = cbind(df, position = near$nn.idx)
    }

    # Find point indexes
    root = as.numeric(df[grepl("root",df$point),"position"])
    primary.branch.point = as.numeric(df[grepl("primary.branch.point",df$point),"position"])
    axon.start = as.numeric(df[grepl("axon.start",df$point),"position"])
    dendrite.start = as.numeric(df[grepl("dendrite.start",df$point),"position"])
    axon.primary = as.numeric(df[grepl("axon.primary",df$point),"position"])
    dendrite.primary = as.numeric(df[grepl("dendrite.primary",df$point),"position"])
    linkers = as.numeric(df[grepl("linker",df$point),"position"])

    # Work around errors
    if(is.na(dendrite.primary)){
      dendrite.primary = dendrite.start[1]
    }
    if(is.na(axon.primary)){
      axon.primary = axon.start[1]
    }

      # Assign root and mark soma
      y = nat::as.neuron(nat::as.ngraph(x$d), origin = root)
      y$connectors = x$connectors
      n = nat::as.ngraph(y$d)
      y$tags$soma = root
      y$d$Label = 0
      y$connectors$Label = 0

      # Assign other cable
      if(!is.issue(root)&!is.issue(primary.branch.point)){
        pnt = suppressWarnings(unique(unlist(igraph::shortest_paths(n, from = as.numeric(root), to = as.numeric(primary.branch.point), mode = "out")$vpath)))
      }else{
        pnt = NULL
      }
      if(!is.issue(linkers)){
        pd = suppressWarnings(unique(unlist(igraph::shortest_paths(n, from = as.numeric(linkers[1]), to = as.numeric(linkers[length(linkers)]), mode = "all")$vpath)))
      }else if(!is.issue(axon.primary)&!is.issue(dendrite.primary)){
        pd = suppressWarnings(unique(unlist(igraph::shortest_paths(n, from = as.numeric(dendrite.primary), to = as.numeric(axon.primary), mode = "all")$vpath)))
      }else{
        pd = NULL
      }
      leaves = as.numeric(nat::endpoints(y))
      if(!is.issue(axon.primary)){
        axon = suppressWarnings(
          unique(unlist(sapply(as.numeric(axon.start),function(s) unlist(igraph::shortest_paths(n, from=s, to = leaves, mode = c("out"))$vpath))))
        )
      }else{
        axon = NULL
      }
      if(!is.issue(dendrite.primary)){
        dendrite = suppressWarnings(
          unique(unlist(sapply(dendrite.start,function(s) unlist(igraph::shortest_paths(n, from=s, to = leaves, mode = c("out"))$vpath))))
        )
      }else{
        dendrite = NULL
      }

      # Assign labels
      if(!is.null(axon)){
        y = add_Label(x = y, PointNo = axon, Label = 2, erase = TRUE)
      }
      if(!is.null(dendrite)){
        y = add_Label(x = y, PointNo = dendrite, Label = 3, erase = TRUE)
      }
      if(!is.null(pd)){
        y = add_Label(x = y, PointNo = pd, Label = 4, erase = TRUE)
      }
      if(!is.null(pnt)){
        y = add_Label(x = y, PointNo = pnt, Label = 7, erase = TRUE)
      }
      if(!is.null(root)){
        y = add_Label(x = y, PointNo = root, Label = 1, erase = TRUE)
      }
      #y$d[dendrite,]$Label = 3

      # Calculate segregation score

      # Add in branch points
      y$primary.branch.point = primary.branch.point
      y$axon.start = axon.start
      y$dendrite.start = dendrite.start
      y$secondary.branch.points = c(axon.start,dendrite.start)

      # Assign bodyid
      y$bodyid = y$d$bodyid = x$bodyid
      y$split = TRUE
  }
  # Return split skeleton
  y = hemibrain_neuron_class(y)
  y
}

#' @export
hemibrain_use_splitpoints.neuronlist <-function(x, df, knn = FALSE, ...){
  x = add_field_seq(x,x[,"bodyid"],field="bodyid")
  neurons = nat::nlapply(x, hemibrain_use_splitpoints.neuron,
                         df = df,
                         knn = knn,
                         ...)
  neurons
}

#' Split neurons into axon and dendrites using pre-computed split points
#'
#' @description Split a neuron into its putative axon and dendrite,
#'  using a \code{data.frame} of precomputed split points. We have already pre-computed
#'  splits for all neurons and stored them in this package as \code{\link{hemibrain_precomputed_splitpoints}}.
#'  This is helpful because  \code{\link{flow_centrality}} can take a long time to run on a large
#'  number of neurons.
#'
#' @inheritParams flow_centrality
#' @param splitpoints a \code{data.frame} of splitpoints from running \code{\link{flow_centrality}},
#' as produced by \code{\link{hemibrain_splitpoints}}. If a custom set of splitpoints is not given, precomputed splitpoints
#' are used, \code{\link{hemibrain_precomputed_splitpoints}}. This defaults to \code{hemibrain_splitpoints},
#' however to see the available precomputations (which have used \code{\link{flow_centrality}}, in this case `polypre = TRUE`, `mode = "centrifugal"`
#' and `split = "distance"`) with different parameters)
#' please see \code{\link{hemibrain_precomputed_splitpoints}}.
#' @inheritParams hemibrain_use_splitpoints
#' @param knn logical, whether or not to find corresponding points
#' between \code{splitpoints} and \code{x$d} using a nearest neighbour search in 3D
#' space (TRUE) or just use the given point IDs in \code{splitpoints} (i.e. if neurons
#' have not been resampled or their skeletons otherwise modified).
#'
#' @inherit flow_centrality return details references
#'
#' @return a \code{neuronlist}
#' @seealso \code{\link{flow_centrality}}
#' @examples
#' \donttest{
#'
#' # Choose some neurons
#' exemplars = c("202916528", "1279775082",  "203253072",
#' "326530038",  "203253253", "5813079341")
#'
#' # Get neurons
#' neurons = neuprintr::neuprint_read_neurons(exemplars)
#'
#' # Now use a pre-saved axon-dendrite split
#' neurons.flow = hemibrain_flow_centrality(neurons)
#'
#' \dontrun{
#' # Plot the split to check it
#' nat::nopen3d()
#' nlscan_split(neurons.flow, WithConnectors = TRUE)
#' }}
#' @export
#' @seealso \code{\link{hemibrain_splitpoints}}, \code{\link{flow_centrality}}, \code{\link{hemibrain_use_splitpoints}}, \code{\link{hemibrain_precomputed_splitpoints}}
hemibrain_flow_centrality <-function(x,
                                     splitpoints = hemibrainr::hemibrain_all_splitpoints,
                                     knn = FALSE,
                                     ...) UseMethod("hemibrain_flow_centrality")

#' @export
hemibrain_flow_centrality.neuron <- function(x, splitpoints = hemibrainr::hemibrain_all_splitpoints, knn = FALSE, ...){
  bi = x$bodyid
  if(is.null(x$bodyid)|is.na(x$bodyid)){
    stop("No bodyid given at x$bodyid")
  }
  df = dplyr ::filter(splitpoints, .data$bodyid == bi)
  if(!nrow(df)){
    warning("bodyid ", bi," not in splitpoints, returning neuron unmodified")
    y = x
  }else{
    y = hemibrain_use_splitpoints(x, df, knn = knn, ...)
  }
  y
}

#' @export
hemibrain_flow_centrality.neuronlist <- function(x, splitpoints = hemibrainr::hemibrain_all_splitpoints, knn = FALSE, ...){
  cropped = subset(x, cropped)
  if(length(cropped)){
    warning(length(cropped), " neurons cropped, split likely to be inaccurate for: ", paste(names(cropped),collapse=", "))
  }
  untraced = subset(x, status!="Traced")
  if(length(untraced)){
    warning(length(untraced), " neurons do not have 'traced' status, split likely to be inaccurate for: ", paste(names(untraced),collapse=", "))
  }
  nosoma = subset(x, !soma)
  if(length(untraced)){
    warning(length(nosoma), " neurons have no soma tagged, split could be inaccurate for: ", paste(names(nosoma),collapse=", "))
  }
  missed = setdiff(x[,"bodyid"], splitpoints$bodyid)
  if(length(missed)>0){
    warning(length(missed), " neurons are not represented in the given splitpoints: ", paste(names(missed),collapse=", "))
  }
  y = hemibrain_use_splitpoints(x, splitpoints, knn = knn, ...)
  y
}

#' Manually add a Label annotation to a neuron
#'
#' @description Assign all points on a skeleton with a certain Label.
#'
#' @inheritParams flow_centrality
#' @param Label the Label to be added. See \code{\link{flow_centrality}}.
#' @param PointNo the points in the neuron for which the \code{Label}
#' will be added.If \code{NULL} all points are used.
#' @param erase if \code{TRUE}, all instance of \code{Label} not
#' in \code{PointNo} are set to \code{0}.
#'
#' @return a \code{neuron} or \code{neuronlist}
#' @seealso \code{\link{flow_centrality}}, \code{\link{add_field}}
#' @export
add_Label <-function(x, PointNo = NULL, Label = 2, erase = FALSE, ...) UseMethod("add_Label")

#' @export
add_Label.neuron <- function(x, PointNo = NULL, Label = 2, erase = FALSE, ...){
  Label = as.numeric(Label[1])
  if(!is.null(PointNo)){
    if(erase){
      erasure = x$d$PointNo[x$d$Label==Label]
      x$d$Label[match(erasure, x$d$PointNo)] = 0
    }
    x$d$Label[x$d$PointNo%in%PointNo] = Label
    if(!is.null(x$connectors)){
      if(erase){
        x$connectors$Label[x$connectors$treenode_id%in%erasure] = 0
      }
      x$connectors$Label[x$connectors$treenode_id%in%PointNo] = Label
    }
  }else{
    x$d$Label = Label
    if(!is.null(x$connectors)){
      x$connectors$Label = Label
    }
  }
  x
}

#' @export
add_Label.neuronlist <- function(x, PointNo = NULL, Label = 2, erase = FALSE, ...){
  nat::nlapply(x,
               add_Label.neuron,
               PointNo = PointNo,
               Label = Label,
               erase = erase,
               ...)
}

#' Add a field to neuron objects (inc. in a neuronlist)
#'
#' @description At an item to the the list object that comprises a neuron
#'
#' @inheritParams flow_centrality
#' @param entry item to add.
#' @param entries a vector/list of entries, the same length as \code{x}.
#' Each neuron in \code{x} will have the entry in entries as the same index, added to it
#' as \code{field}.
#' @param field name of new field, or field to overwrite.
#'
#' @return a \code{neuron} or \code{neuronlist}
#' @seealso \code{\link{add_Label}}
add_field <-function(x, entry, field = "bodyid", ...) UseMethod("add_field")
#' @export
add_field.neuron <- function(x, entry, field = "bodyid", ...){
  x[[field]] = entry
  x
}
#' @export
add_field.neuronlist <- function(x, entry, field = "bodyid", ...){
  nat::nlapply(x, add_Label.neuron, entry, field = "bodyid", ...)
}
#' @export
#' @rdname add_field
add_field_seq <- function(x, entries, field = "bodyid", ...){
  x = nat::as.neuronlist(x)
  if(length(entries)!=length(x)){
    stop("The length of the entries to add must be the same as the length of the neuronlist, x")
  }
  nl = nat::neuronlist()
  for(i in 1:length(x)){
    y = x[[i]]
    entry = entries[i]
    y = add_field(y, entry = entry, field = field, ...)
    nl = c(nl, nat::as.neuronlist(y))
  }
  names(nl) = names(x)
  nl[,] = x[,]
  nl
}

# # test set of tricky neurons:
# friends = c("327499164", "328861282", "487144598","480590566","574688051",
#             "514375643", "5813087438", "421641859", "604709727","328533761",
#             "329225149", "329897255","330268940", "5813068669","360255138",
#             "360284300", "511271574","579912201","5813021291", "5813075020",
#             "517506265","5813020988",
#             "1203070528","5812982779","925799233",
#             "1141976953", "1096919756", "364061776", "486837442", "693500652",
#             "846952407", "860266123", "664511977", "611646132", "632069728",
#                "1295855722","455159380","456174330","613398716",
#                "488559952","5813096699","1234481054","487899562",
#                "549955973","1233773217","1204124270","1233428010",
#                "581332573","1109957857","1266871604","1047525093",
#                "1202397143","582355376", "730562988","5813056323",
#                "579912201", "5813015982", "973765182", "885788485",
#                "915451074", "5813032740", "1006854683", "5813013913", "5813020138",
#                "853726809", "916828438", "5813078494", "420956527", "486116439",
#                "573329873", "5813010494", "5813040095", "514396940", "665747387",
#                "793702856", "451644891", "482002701", "391631218", "390948259",
#                "390948580", "452677169", "511262901", "422311625", "451987038",
#                 "612738462","487925037"
# )
# neurons = neuprint_read_neurons(friends[23:length(friends)])
# hemibrain.rois = hemibrain_roi_meshes()
# neurons.checked = hemibrain_skeleton_check(neurons, meshes = hemibrain.rois)
# neurons.flow = flow_centrality(neurons.checked, polypre = TRUE,
#                                mode = "centrifugal",
#                                split = "synapses",
#                                primary.branchpoint = 0.25)
# nlscan_split(neurons.flow, WithConnectors = TRUE)
# # 20, 65
#

# al.local.neurons = c("1702323386", "2068966051", "2069311379", "1702305987", "5812996027",
#  "1702336197", "1793744512", "1976565858", "2007578510", "2101339904",
#  "5813003258", "2069647778", "1947192569", "1883788812", "1916485259",
#  "1887177026", "2101348562", "2132375072", "2256863785", "5813002313",
#  "5813054716", "5813018847", "5813055448", "1763037543", "2101391269",
#  "1794037618", "5813018729", "2013333009")
# neurons = neuprint_read_neurons(al.local.neurons)
# hemibrain.rois = hemibrain_roi_meshes()
# neurons.checked = hemibrain_skeleton_check(neurons, meshes = hemibrain.rois)
# neurons.flow = flow_centrality(neurons.checked, polypre = TRUE,
#                                mode = "centrifugal",
#                                split = "synapses",
#                                primary.branchpoint = 0.25)
# syns = nat::nlapply(neurons.flow, extract_synapses, unitary = FALSE)
# syns = do.call(rbind,syns)
# syns = subset(syns, partner %in% al.local.neurons)
# nlscan_split(neurons.flow, WithConnectors = TRUE)


# "5901222683"

