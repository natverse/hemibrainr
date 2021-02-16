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

