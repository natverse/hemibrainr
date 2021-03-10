#
#
#
#
# # Get flywire connectivty matrix
# flywire_elist <- function(x,
#                                          local = FALSE,
#                                          cloudvolume.url = NULL){
#   method = match.arg(method)
#
#   # Get Google drive folder
#   savedir = good_savedir(local = local)
#
#   # we need the hemibrain, not hemibrainr drive
#   # this is a temporary trick to get it
#   savedir = gsub("hemibrainr","hemibrain",savedir)
#
#   # Fetch
#   fw.syns = find_gsql(savedir,
#                       sql.db = "flywire_synapses.db",
#                       tab = NULL,
#                       folder = "fafbsynapses")
#
#   # Get rootids
#   rootids = x[,"flywire.id"]
#
#   # Assemble synapses
#   for(rootid in unique(rootids)){
#     message("Fetching supervoxel ids for id: ", rootid)
#     svids = get_flywire_svids(rootid)
#     syn.map = flywire_svid_synapses(svids = svids,
#                                     partners = "both",
#                                     method = method,
#                                     Verbose = TRUE)
#
#   }
#
#
#
# }
#
#


# Simplify synapses spatially
flywire_synapse_simplify <- function(x, method = c("cleft_scores","scores","mean")){
  method = match.arg(method)

  # for each neuron
  if("pre_id"%in%colnames(x)){
    if(length(unique(x$pre_id))>1){
      res = list()
      for(preid in unique(x$pre_id)){
        z = subset(x, x$pre_id==preid)
        y = flywire_synapse_simplify(x=z, method=method)
        res[[preid]]=y
      }
      return(do.call(plyr::rbind,res))
    }
  }

  # Get 3D points
  x = x[order(x$scores,decreasing=TRUE),]
  x = x[order(x$cleft_scores,decreasing=TRUE),]
  if("prepost"%in%colnames(x)){
    pre = subset(x, x$prepost==0)
    post = subset(x, x$prepost==1)
  }else{
    pre = x
    post = data.frame()
  }
  xyz = tryCatch(nat::xyzmatrix(pre), error = function(e){
    pre[,c("x","y","z")] = pre[,c("pre_x","pre_y","pre_z")]
    nat::xyzmatrix(pre)
  })

  # Filter synapses
  near = nabor::knn(query = xyz, data = xyz, k = nrow(xyz), radius = 5000)
  near$nn.dists[is.infinite(near$nn.dists)] = 5000
  hc = hclust(dist(near$nn.dists))
  ct = dendextend::cutree(tree=hc, h = 2500)
  pre$cluster = ct

  # For synapse in cluster, combine
  pre.new = data.frame()
  for(t in 1:max(ct)){
    np = subset(pre, pre$cluster==t)
    if(method=="mean"){
      nl = unlist(lapply(np,class))%in%c("numeric","integer")
      collap = as.data.frame(t(colMeans(np[,nl], na.rm = TRUE)))
      collap$offset = np$offset[1]
      poss.nts=c("gaba", "acetylcholine", "glutamate", "octopamine", "serotonin","dopamine")
      tops = colSums(collap[,poss.nts])
      collap$top.p = max(tops)
      collap$top.nt = names(which.max(tops))
    }else{
      np = np[order(np[[method]],decreasing=TRUE),]
      collap = np[1,]
    }
    if(nrow(post)){
      post$pre_id%in%np$pre_id = collap$pre_id
      post$pre_svid%in%np$pre_svid = collap$pre_svid
    }
    pre.new = rbind(pre.new,collap)
  }
  post[!(duplicated(post$post_svid)&duplicated(post$pre_svid)),]
  plyr::rbind.fill(pre.new,post)
}




