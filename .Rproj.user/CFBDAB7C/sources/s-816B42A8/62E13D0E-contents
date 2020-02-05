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
#' @param x a \code{nat::neuronlist} or \code{nat::neuron} object
#' @param mode type of flow centrality to calculate. There are three flavors:
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
#'   dendrite Label = 3, primay dendrite Label = 9 and
#'   primary neurite Label = 7.
#'   Soma is Label = 1.
#'
#' @examples
#' \donttest{
#' \dontrun{
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
#' neurons.checked = hemibrain_skeleton_check(neurons, meshes = hemibrain.rois)
#'
#' # Split neuron
#' ## These are the recommended parameters for hemibrain neurons
#' neurons.flow = flow_centrality(neurons.checked, polypre = TRUE,
#' mode = "centrifugal",
#' split = "distance")
#'
#' # Plot the split to check it
#' nat::nopen3d()
#'
#' }}
#' @export
#' @seealso \code{\link{primary_neurite}}, \code{\link{hemibrain_skeleton_check}}
flow_centrality <-function(x,
                           mode = c("sum","centrifugal","centripetal"),
                           polypre = TRUE,
                           soma = TRUE,
                           primary.dendrite = 0.9,
                           bending.flow = FALSE,
                           split = c("distance","postsynapses","presynapses"),
                           ...) UseMethod("flow_centrality")

# hidden
flow_centrality.neuron <- function(x,
                                   mode = c("sum", "centrifugal", "centripetal"),
                                   polypre = TRUE,
                                   soma = TRUE,
                                   primary.dendrite = 0.9,
                                   bending.flow = FALSE,
                                   split = c("postsynapses","presynapses", "distance")){
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
    syns.out = unique(x$connectors[x$connectors$prepost ==0, ][, "treenode_id"])
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
    ais = ais[match(min(runstosoma), runstosoma)]
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
  downstream.unclassed = downstream[!downstream %in% c(p.n,highs)]
  remove = rownames(nodes)[!rownames(nodes) %in% downstream.unclassed]
  downstream.g = igraph::delete_vertices(n, v = as.character(remove))
  main1 = igraph::components(downstream.g)
  main1 = names(main1$membership[main1$membership %in% 1])
  nodes.downstream = nodes[as.character(main1), ]
  tract.parent = unique(nodes.downstream$Parent[!nodes.downstream$Parent %in% nodes.downstream$PointNo])
  downstream.tract.parent = match(tract.parent, nodes$PointNo)
  if (sum(match(tract.parent, nodes$PointNo) %in% p.n) > 0) {
    bps.all = rownames(nodes)[match(as.numeric(nat::branchpoints(nodes)),nodes$PointNo)]
    bps.downstream = bps.all[bps.all %in% downstream.unclassed]
    runstoprimarybranchpoint = unlist(lapply(bps.downstream,
                                             function(x) length(unlist(suppressWarnings(igraph::shortest_paths(n,to = downstream.tract.parent, from = x)$vpath)))))
    downstream.tract.parent = bps.downstream[which.min(runstoprimarybranchpoint)]
  }
  downstream.tract.parent = nodes[downstream.tract.parent,]
  upstream.unclassed = upstream[!upstream %in% c(p.n, highs)]
  remove = rownames(nodes)[!rownames(nodes) %in% upstream.unclassed]
  upstream.g = igraph::delete_vertices(n, v = as.character(remove))
  main1 = igraph::components(upstream.g)
  main1 = names(main1$membership[main1$membership %in% 1])
  nodes.upstream = nodes[as.character(main1), ]
  tract.parent = unique(nodes.upstream$Parent[!nodes.upstream$Parent %in%nodes.upstream$PointNo])
  upstream.tract.parent = ifelse(tract.parent==-1,root,match(tract.parent, nodes$PointNo))
  if (sum(match(tract.parent, nodes$PointNo) %in% p.n) > 0) {
    bps.all = rownames(nodes)[match(as.numeric(nat::branchpoints(nodes)),nodes$PointNo)]
    bps.upstream = bps.all[bps.all %in% upstream.unclassed]
    runstoprimarybranchpoint = unlist(lapply(bps.upstream,function(x) length(unlist(suppressWarnings(igraph::shortest_paths(n,to = upstream.tract.parent, from = x)$vpath)))))
    upstream.tract.parent = bps.upstream[which.min(runstoprimarybranchpoint)]
  }
  upstream.tract.parent = nodes[upstream.tract.parent, ]
  neurite.nodes = nodes[!rownames(nodes) %in% p.n, ]
  p.n.PointNo = nodes[p.n, "PointNo"]
  primary.branch.point = p.n[p.n.PointNo %in% neurite.nodes$Parent]
  if (grepl("synapses", split)) {
    synapse.choice = gsub("synapses", "", split)
    choice = mean(nodes[as.character(downstream.unclassed),
                        synapse.choice], na.rm = TRUE) < mean(nodes[as.character(upstream.unclassed),
                                                                    synapse.choice], na.rm = TRUE)
    if (choice) {
      nodes[as.character(downstream.unclassed), "Label"] = 2
    }
    else if (!choice) {
      nodes[as.character(upstream.unclassed), "Label"] = 2
    }
    else {
      split == "distance"
      warning("synapse numbers are the same, splitting based on branch point distances to primary branchpoint")
    }
  }
  if (split == "distance") {
    primary.branch.point.xyz = as.matrix(nat::xyzmatrix(nodes[primary.branch.point, ]))
    secondary.branch.points.xyz = nat::xyzmatrix(rbind(upstream.tract.parent,downstream.tract.parent))
    dist.upstream.to.primary.branchpoint = length(unlist(igraph::shortest_paths(n,to = primary.branch.point, from = rownames(upstream.tract.parent))$vpath))
    dist.downstream.to.primary.branchpoint = length(unlist(igraph::shortest_paths(n,to = primary.branch.point, from = rownames(downstream.tract.parent))$vpath))
    if (dist.upstream.to.primary.branchpoint < dist.downstream.to.primary.branchpoint) {
      nodes[as.character(downstream.unclassed), "Label"] = 2
    }
    else if (dist.upstream.to.primary.branchpoint > dist.downstream.to.primary.branchpoint) {
      nodes[as.character(upstream.unclassed), "Label"] = 2
    }
    else {
      warning("branch point distances are the same, splitting based on postsynapses")
      choice = sum(nodes[as.character(downstream.unclassed),
                         "post"]) < sum(nodes[as.character(upstream.unclassed),
                                              "post"])
      if (choice) {
        nodes[as.character(downstream.unclassed), "Label"] = 2
      }
      else {
        nodes[as.character(upstream.unclassed), "Label"] = 2
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
  x$d = nodes
  x$AD.segregation.index = segregation.index
  x$primary.branch.point = as.numeric(primary.branch.point)
  x$secondary.branch.points = as.numeric(c(downstream.tract.parent$PointNo, upstream.tract.parent$PointNo))
  x$max.flow.centrality = as.numeric(ais)
  x
}

# hidden
flow_centrality.neuronlist <- function(x,
                                       mode = c("sum","centrifugal","centripetal"),
                                       polypre = T,
                                       soma = TRUE,
                                       primary.dendrite = 0.9,
                                       bending.flow = FALSE,
                                       split = c("postsynapses","presynapses","distance"),
                                       ...){
  neurons = nat::nlapply(x, flow_centrality, mode = mode, polypre = polypre, soma = soma, primary.dendrite = primary.dendrite, OmitFailures = T, split = split, ...)
  neurons
}
