### Code to prepare datasets ###

### Split points
hemibrain_splitpoints_polypre_centrifugal_distance <- read.csv("data-raw/hemibrain_all_neurons_splitpoints_polypre_centrifugal_distance.csv")
hemibrain_splitpoints_pre_centrifugal_distance <- read.csv("data-raw/hemibrain_all_neurons_splitpoints_pre_centrifugal_distance.csv")
hemibrain_splitpoints_polypre_centrifugal_synapses <- read.csv("data-raw/hemibrain_all_neurons_splitpoints_polypre_centrifugal_synapses.csv")
usethis::use_data(hemibrain_splitpoints_pre_centrifugal_distance, overwrite = TRUE)
usethis::use_data(hemibrain_splitpoints_polypre_centrifugal_distance, overwrite = TRUE)
usethis::use_data(hemibrain_splitpoints_polypre_centrifugal_synapses, overwrite = TRUE)

### Metrics
hemibrain_metrics_polypre_centrifugal_distance <- read.csv("data-raw/hemibrain_all_neurons_metrics_polypre_centrifugal_distance.csv")
hemibrain_all_neurons_metrics_polypre_centrifugal_synapses <- read.csv("data-raw/hemibrain_all_neurons_metrics_polypre_centrifugal_synapses.csv")
rownames(hemibrain_metrics_polypre_centrifugal_distance) <- hemibrain_metrics_polypre_centrifugal_distance$bodyid
rownames(hemibrain_all_neurons_metrics_polypre_centrifugal_synapses) <- hemibrain_all_neurons_metrics_polypre_centrifugal_synapses$bodyid
usethis::use_data(hemibrain_metrics_polypre_centrifugal_distance, overwrite = TRUE)
usethis::use_data(hemibrain_all_neurons_metrics_polypre_centrifugal_synapses, overwrite = TRUE)

### Consensus
hemibrain_all_splitpoints <- hemibrain_splitpoints_polypre_centrifugal_distance
hemibrain_metrics <- hemibrain_metrics_polypre_centrifugal_distance
usethis::use_data(hemibrain_metrics, overwrite = TRUE)
usethis::use_data(hemibrain_all_splitpoints, overwrite = TRUE)

### Surface mesh
hemibrain.surf = readobj::read.obj("data-raw/hemibrain_raw.obj", convert.rgl = TRUE)
hemibrain.surf = nat::as.hxsurf(hemibrain.surf[[1]])
hemibrain_microns.surf = hemibrain.surf*(8/1000)
nat.templatebrains::regtemplate(hemibrain_microns.surf) = "JRCFIB2018F"
nat.templatebrains::regtemplate(hemibrain.surf) = "JRCFIB2018Fraw"
usethis::use_data(hemibrain_microns.surf, overwrite = TRUE)
usethis::use_data(hemibrain.surf, overwrite = TRUE)

### AL glomeruli meshes
al.objs = list.files("data-raw/AL/",pattern = ".obj",full.names = TRUE)
al.meshes = sapply(al.objs, readobj::read.obj, convert.rgl = TRUE)
gloms = gsub(".*raw\\.|\\_mesh.*","",names(al.meshes))
names(al.meshes) = gloms
Vertices = data.frame()
Regions = list()
count = 0
for(i in 1:length(al.meshes)){
  glom = names(al.meshes)[i]
  mesh = al.meshes[[glom]]
  h = nat::as.hxsurf(mesh)
  vertices = h$Vertices
  vertices$PointNo = vertices$PointNo+count
  regions = h$Regions$Interior + count
  count = max(vertices$PointNo)
  Vertices = rbind(Vertices, vertices)
  Regions[[glom]] = regions
}
hemibrain_al.surf = list(Vertices = Vertices,
     Regions = Regions,
     RegionList = gloms,
     RegionColourList = hemibrain_bright_colour_ramp(length(gloms)))
class(hemibrain_al.surf) = "hxsurf"
hemibrain_al_microns.surf = hemibrain_al.surf*(8/1000)
nat.templatebrains::regtemplate(hemibrain_al_microns.surf) = "JRCFIB2018F"
nat.templatebrains::regtemplate(hemibrain_al.surf) = "JRCFIB2018Fraw"
usethis::use_data(hemibrain_al_microns.surf, overwrite = TRUE)
usethis::use_data(hemibrain_al.surf, overwrite = TRUE)

# Just save the ORN and HRN bodyids
rn.ids = class2ids("RN", possible = TRUE)
orn.ids = class2ids("ORN")
hrn.ids = class2ids("HRN")
pn.ids = class2ids("PN")
upn.ids = class2ids("uPN")
mpn.ids = class2ids("mPN")
vppn.ids = class2ids("VPPN")
alln.ids = class2ids("ALLN")
dan.ids = class2ids("DAN")
mbon.ids = class2ids("MBON")
## Use them
usethis::use_data(rn.ids, overwrite = TRUE)
usethis::use_data(orn.ids, overwrite = TRUE)
usethis::use_data(hrn.ids, overwrite = TRUE)
usethis::use_data(pn.ids, overwrite = TRUE)
usethis::use_data(upn.ids, overwrite = TRUE)
usethis::use_data(mpn.ids, overwrite = TRUE)
usethis::use_data(vppn.ids, overwrite = TRUE)
usethis::use_data(alln.ids, overwrite = TRUE)
usethis::use_data(dan.ids, overwrite = TRUE)
usethis::use_data(dan.ids, overwrite = TRUE)
usethis::use_data(mbon.ids, overwrite = TRUE)


