#######################################################################################
################################ axon-dendrite split ##################################
#######################################################################################

#' Determine a neuron's dendrite and axon by calculating flow centrality
#'
#' @description implementation of the algorithm for calculating flow
#'   centralities from Schneider-Mizell et al. (2016).
#'   Can be used on \code{nat::neuronlist}
#'   objects read using \code{\link{neuprint_read_neurons}}. The example code
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
#' @param bending.flow we may need to add the 'bending flow' to all
#' the branchpoints if looking at centripetal flow centrality.
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
#'
#' # Get neurons
#' neurons = neuprint_read_neurons(tough)
#'
#' # Now make sure the neurons have a soma marked
#' ## Some hemibrain neurons do not, as the soma was chopped off
#' neurons.checked = hemibrain_skeleton_check(neurons, meshes = hemibrain.surf)
#'
#' # Split neuron
#' ## These are the recommended parameters for hemibrain neurons
#' neurons.flow = flow_centrality(neurons.checked, polypre = FALSE,
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
                           polypre = FALSE,
                           soma = TRUE,
                           primary.dendrite = 0.9,
                           bending.flow = FALSE,
                           split = c("distance","postsynapses","presynapses"),
                           ...) UseMethod("flow_centrality")

#' @export
flow_centrality.neuron <- function(x,
                                   mode = c("centrifugal", "centripetal", "sum"),
                                   polypre = FALSE,
                                   soma = TRUE,
                                   primary.dendrite = 0.9,
                                   bending.flow = FALSE,
                                   split = c("distance", "postsynapses","presynapses"),
                                   ...){
  split = match.arg(split)
  mode = match.arg(mode)
  x$d$Label = 0
  el = x$d[x$d$Parent != -1, c("Parent", "PointNo")]
  n = nat::ngraph(data.matrix(el[, 2:1]), x$d$PointNo, directed = TRUE, xyz = nat::xyzmatrix(x$d), diam = x$d$W)
  leaves = which(igraph::degree(n, v = igraph::V(n), mode = "in") == 0, useNames = T)
  root = which(igraph::degree(n, v = igraph::V(n), mode = "out") == 0, useNames = T)
  segs = x$SegList
  nodes = x$d
  nodes[, c("post","pre","up.syns.in","up.syns.out","flow.cent")] = 0
  nodes[,"Label"] = 3
  nodes = nodes[unlist(c(root, lapply(segs, function(x) x[-1]))),]
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
  in.total = nodes[1, "up.syns.in"] = length(point.no.in)
  out.total = nodes[1, "up.syns.out"] = length(point.no.out)
  if (mode == "centrifugal") {
    nodes[, "flow.cent"] = (in.total - nodes[, "up.syns.in"]) * nodes[, "up.syns.out"]
  }else if (mode == "centripetal") {
    nodes[, "flow.cent"] = (out.total - nodes[, "up.syns.out"]) * nodes[, "up.syns.in"]
  }else {
    nodes[, "flow.cent"] = ((in.total - nodes[, "up.syns.in"]) * nodes[, "up.syns.out"]) + ((out.total - nodes[, "up.syns.out"]) * nodes[, "up.syns.in"])
  }
  nodes = nodes[order(as.numeric(rownames(nodes))), ]
  if (bending.flow) {
    for (bp in bps) {
      down = unlist(igraph::ego(n, 1, nodes = bp, mode = "in", mindist = 0))[-1]
      bending.flow = centrifugal.bending.flow = c()
      for (u in down) {
        this.seg.posts = nodes[u, ]$down.syns.in
        other.segs.pre = nodes[down[!down == u], ]$down.syns.out
        bending.flow = c(bending.flow, sum(this.seg.posts *
                                             other.segs.pre))
      }
      nodes[bp, "flow.cent"] = nodes[bp, "flow.cent"] +  bending.flow
    }
  }
  ais = which(apply(nodes, 1, function(x) x["flow.cent"] == max(nodes[, "flow.cent"])))
  if (length(ais) > 0) {
    runstosoma = unlist(lapply(ais, function(x) length(unlist(igraph::shortest_paths(n,x, to = root)$vpath))))
    ais = as.numeric(ais[match(min(runstosoma), runstosoma)])
  }
  downstream = suppressWarnings(unique(unlist(igraph::shortest_paths(n, ais, to = leaves, mode = "in")$vpath)))
  upstream = rownames(nodes)[!rownames(nodes) %in% downstream]
  if (bending.flow) {
    if (nodes[ais, ]$up.syns.in == 0 & nodes[ais, ]$up.syns.out == 0) {
      down = unlist(igraph::ego(n, 1, nodes = bp, mode = "in",mindist = 0))[-1]
      ais = down[1]
      downstream = suppressWarnings(unique(unlist(igraph::shortest_paths(n,ais, to = leaves, mode = "in")$vpath)))
      upstream = rownames(nodes)[!rownames(nodes) %in%downstream]
    }
  }
  igraph::V(n)$name = igraph::V(n)
  if(soma){
    p.n = primary_neurite(x, neuron = FALSE)
  }else{
    p.n = unlist(x$SegList[which.max(sapply(x$SegList, function(x) sum(nodes[x,"flow.cent"]==0)))])
  }
  nodes[p.n, "Label"] = 7
  if (!is.null(primary.dendrite)) {
    highs = subset(rownames(nodes), nodes[, "flow.cent"] >=primary.dendrite * max(nodes[, "flow.cent"]))
    nodes[as.character(highs), "Label"] = 4
  }else {
    primary.dendrite = 0.9
    highs = subset(rownames(nodes), nodes[, "flow.cent"] >=primary.dendrite * max(nodes[, "flow.cent"]))
  }
  p.n = p.n[1:ifelse(sum(p.n%in%highs)>1,min(which(p.n%in%highs)),length(p.n))]
  primary.branch.points = p.n[p.n%in%nat::branchpoints(nodes)]
  primary.branch.point = primary.branch.points[length(primary.branch.points)]
  primary.branch.point.downstream = suppressWarnings(unique(unlist(igraph::shortest_paths(n, as.numeric(primary.branch.points), to = as.numeric(leaves), mode = "in")$vpath)))
  downstream.unclassed = downstream[!downstream %in% c(p.n,highs, root, leaves, primary.branch.point)]
  remove = rownames(nodes)[!rownames(nodes) %in% intersect(downstream.unclassed,primary.branch.point.downstream)]
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
  upstream.unclassed = upstream[!upstream %in% c(p.n, highs, root, leaves, primary.branch.point)]
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
  if (grepl("synapses", split)) {
    synapse.choice = gsub("synapses", "", split)
    sdown = sum(nodes[as.character(downstream.unclassed),
                      synapse.choice], na.rm = TRUE)
    sup = sum(nodes[as.character(upstream.unclassed),
                    synapse.choice], na.rm = TRUE)
    choice = sdown > sup
    if (sdown == sup) {
      split == "distance"
      warning("synapse numbers are the same, splitting based on branch point distances to primary branchpoint")
    }
    else if ( (choice & synapse.choice=="pre") | (!choice & synapse.choice=="post") ) {
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
    dist.upstream.to.primary.branchpoint = length(unlist(igraph::shortest_paths(n,to = as.numeric(primary.branch.point), from = as.numeric(upstream.tract.parent))$vpath))
    dist.downstream.to.primary.branchpoint = length(unlist(igraph::shortest_paths(n,to = as.numeric(primary.branch.point), from = as.numeric(downstream.tract.parent))$vpath))
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
  starts = sapply(c(p.n,highs),function(s) igraph::neighbors(n, v=s, mode = c("in")))
  starts = unique(unlist(starts))
  axon.starts = as.numeric(axon.nodes[axon.nodes%in%starts])
  dendrite.starts = as.numeric(dendrite.nodes[dendrite.nodes%in%starts])
  x$d = nodes
  x$AD.segregation.index = segregation.index
  x$primary.branch.point = as.numeric(primary.branch.point)
  x$axon.start = nullToNA(axon.starts)
  x$dendrite.start = nullToNA(dendrite.starts)
  nsbp <- nodes[secondary.branch.points,]
  x$axon.primary = nullToNA(rownames(subset(nsbp, nsbp$Label==2)))
  x$dendrite.primary = nullToNA(rownames(subset(nsbp, nsbp$Label==3)))
  x$secondary.branch.points = secondary.branch.points
  x$max.flow.centrality = as.numeric(ais)
  x
}

#' @export
flow_centrality.neuronlist <- function(x,
                                       mode = c("centrifugal", "centripetal", "sum"),
                                       polypre = FALSE,
                                       soma = TRUE,
                                       primary.dendrite = 0.9,
                                       bending.flow = FALSE,
                                       split = c("distance","postsynapses","presynapses"),
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
              dendrite.start = nullToNA(n$dendrite.start))
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
#'
#' # Get neurons
#' neurons = neuprint_read_neurons(tough)
#'
#' # Now make sure the neurons have a soma marked
#' ## Some hemibrain neurons do not, as the soma was chopped off
#' neurons.checked = hemibrain_skeleton_check(neurons, meshes = hemibrain.surf)
#'
#' # Split neuron
#' ## These are the recommended parameters for hemibrain neurons
#' neurons.flow = flow_centrality(neurons.checked, polypre = FALSE,
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

  # Find nearest points on skeleton
  if(knn){
    near = nabor::knn(query = nat::xyzmatrix(df), data = nat::xyzmatrix(x), k  = 1)
    df = cbind(df, position = near$nn.idx)
  }

  # Find point indexes
  root = as.numeric(df[df$point=="root","position"])
  primary.branch.point = as.numeric(df[df$point=="primary.branch.point","position"])
  axon.start = as.numeric(df[grepl("axon.start",df$point),"position"])
  dendrite.start = as.numeric(df[grepl("dendrite.start",df$point),"position"])
  axon.primary = as.numeric(df[grepl("axon.primary",df$point),"position"])
  dendrite.primary = as.numeric(df[grepl("dendrite.primary",df$point),"position"])

  # Assign root and mark soma
  y = nat::as.neuron(nat::as.ngraph(x$d), origin = root)
  y$connectors = x$connectors
  n = nat::as.ngraph(y$d)
  y$tags$soma = root
  y$d$Label = 3 # dendrite by default

  # Assign other cable
  pnt = suppressWarnings(unique(unlist(igraph::shortest_paths(n, from = as.numeric(root), to = as.numeric(primary.branch.point), mode = "out")$vpath)))
  pd = suppressWarnings(unique(unlist(igraph::shortest_paths(n, from = as.numeric(dendrite.primary), to = as.numeric(axon.primary), mode = "all")$vpath)))
  leaves = as.numeric(nat::endpoints(y))
  axon = suppressWarnings(
    unique(unlist(sapply(as.numeric(axon.start),function(s) unlist(igraph::shortest_paths(n, from=s, to = leaves, mode = c("out"))$vpath))))
  )
  # dendrite = suppressWarnings(
  #   unique(unlist(sapply(dendrite.start,function(s) unlist(igraph::shortest_paths(n, from=s, to = leaves, mode = c("out"))$vpath))))
  # )

  # Assign labels
  y$d[axon,]$Label = 2
  #y$d[dendrite,]$Label = 3
  y$d[pd,]$Label = 4
  y$d[pnt,]$Label = 7
  y$d[root,]$Label = 1

  # Calculate segregation score

  # Add in branch points
  y$primary.branch.point = primary.branch.point
  y$axon.start = axon.start
  y$dendrite.start = dendrite.start
  y$secondary.branch.points = c(axon.start,dendrite.start)

  # Assign bodyid
  y$bodyid = y$d$bodyid = x$bodyid

  # Return split skeleton
  y

}

#' @export
hemibrain_use_splitpoints.neuronlist <-function(x, df, knn = FALSE, ...){
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
#' are used, \code{\link{hemibrain_precomputed_splitpoints}}. This defaults to \code{hemibrain_splitpoints_polypre_centrifugal_distance},
#' however to see the available precomputations (which have used \code{\link{flow_centrality}}, in this case `polypre = TRUE`, `mode = "centrifgual"`
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
#' neurons = neuprint_read_neurons(exemplars)
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
                                     splitpoints = hemibrainr::hemibrain_splitpoints_polypre_centrifugal_distance,
                                     knn = FALSE,
                                     ...) UseMethod("hemibrain_flow_centrality")

#' @export
hemibrain_flow_centrality.neuron <- function(x, splitpoints = hemibrainr::hemibrain_splitpoints_polypre_centrifugal_distance, knn = FALSE, ...){
  bi = x$bodyid
  df = filter(splitpoints, .data$bodyid == bi)
  y = hemibrain_use_splitpoints(x, df, knn = knn, ...)
  y$bodyid = bi
  y
}

#' @export
hemibrain_flow_centrality.neuronlist <- function(x, splitpoints = hemibrainr::hemibrain_splitpoints_polypre_centrifugal_distance, knn = FALSE, ...){
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
  y = hemibrain_use_splitpoints(x, splitpoints, knn = knn, ...)
  y
}



#' Manually add a Label annotation to a neuron
#'
#' @description Assign all points on a skeleton with a certain Label.
#'
#' @inheritParams flow_centrality
#' @param Label the Label to be added. See \code{\link{flow_centrality}}.
#'
#' @return a \code{neuron} or \code{neuronlist}
#' @seealso \code{\link{flow_centrality}}, \code{\link{hemibrain_use_splitpoints}}
#' @export
add_Label <-function(x, Label = 2, ...) UseMethod("add_Label")
#' @export
add_Label.neuron <- function(x, Label = 2, ...){
  x$d$Label = Label
  if(!is.null(x$connectors)){
    x$connectors$Label = Label
  }
  x
}
#' @export
add_Label.neuronlist <- function(x, Label = 2, ...){
  nat::nlapply(x, add_Label.neuron, Label = Label, ...)
}


#' Add a field to neuron objects (inc. in a neuronlist)
#'
#' @description At an item to the the list object that comprises a neuron
#'
#' @inheritParams flow_centrality
#' @param entry item to add.
#' @param entries a vector/list of entries, tha same length as \code{x}.
#' Each neuron in x will have the entry in entries as the same index, added to it
#' as \code{field}.
#' @param field name of new field, or field to overwrite.
#'
#' @return a \code{neuron} or \code{neuronlist}
#' @seealso
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
add_field_seq <- function(x, entries, field = "bodyid", ...){
  x = nat::as.neuronlist(x)
  if(length(entries)!=length(x)){
    stop("The length of the entries to add must be the same as the length of the neuronlist, x")
  }
  nl = nat::neuronlist()
  for(i in 1:length(x)){
    y = x[[1]]
    entry = entries[i]
    y = add_field(y, entry = entry, field = field, ...)
    nl = c(nl, nat::as.neuronlist(y))
  }
  names(nl) = names(x)
  nl[,] = x[,]
  nl
}







