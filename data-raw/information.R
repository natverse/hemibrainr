### Code to prepare neuron information ###
library(googledrive)

## Olfactory information
### Depends on having the repo: hemibrain_olf_data
odour_scenes = read.csv("/Users/GD/LMBD/Papers/hemibrain_olf_data/odour_scenes.csv")
odour.scenes = odour_scenes[,c("glomerulus", "key_ligand", "odour_scene", "valence")]
odour.scenes.agg.1 = odour.scenes.agg = aggregate(list(odour_scenes=odour.scenes$odour_scene),
                             list(glomerulus=odour.scenes$glomerulus),
                             function(x) paste(unique(x), collapse="/"))
odour.scenes.agg.2 = aggregate(list(key_ligand=odour.scenes$key_ligand),
                             list(glomerulus=odour.scenes$glomerulus),
                             function(x) paste(unique(x), collapse="/"))
odour.scenes.agg.3 = aggregate(list(valence=odour.scenes$valence),
                             list(glomerulus=odour.scenes$glomerulus),
                             function(x) paste(unique(x), collapse="/"))
hemibrain_glomeruli_summary = read.csv("/Users/GD/LMBD/Papers/hemibrain_olf_data/AL_gloms_RN_info.csv")
hemibrain_glomeruli_summary$receptor = hemibrain_glomeruli_summary$RN
hemibrain_glomeruli_summary = merge(hemibrain_glomeruli_summary, odour.scenes.agg.1, all.x = TRUE, all.y = TRUE)
hemibrain_glomeruli_summary = merge(hemibrain_glomeruli_summary, odour.scenes.agg.2, all.x = TRUE, all.y = TRUE)
hemibrain_glomeruli_summary = merge(hemibrain_glomeruli_summary, odour.scenes.agg.3, all.x = TRUE, all.y = TRUE)
hemibrain_glomeruli_summary$X = hemibrain_glomeruli_summary$RN = NULL

## Olfactory projection neurons
mpns.fib.info = read.csv("/Users/GD/LMBD/Papers/hemibrain_olf_data/FIB_mPNs.csv")
upns.fib.info = read.csv("/Users/GD/LMBD/Papers/hemibrain_olf_data/FIB_uPNs.csv")
upns.left.fib.info = read.csv("/Users/GD/LMBD/Papers/hemibrain_olf_data/FIB_LHS_uPNs.csv")
vpns.fib.info = read.csv("/Users/GD/LMBD/Papers/hemibrain_olf_data/FIB_VP_PNs.csv")
other.fib.info = read.csv("/Users/GD/LMBD/Papers/hemibrain_olf_data/FIB_other_PNs.csv")
pn.info = plyr::rbind.fill(upns.fib.info,mpns.fib.info,vpns.fib.info, other.fib.info, upns.left.fib.info)
pn.info = plyr::rbind.fill(pn.info, upns.left.fib.info)
pn.info = pn.info[!is.na(pn.info$bodyid),]
pn.info = pn.info[!duplicated(pn.info$bodyid),]
pn.info = subset(pn.info, !grepl("_L",pn.info$instance)) # Just take right-side
pn.info$glomerulus = gsub("\\+.*","",pn.info$glomerulus)
pn.info = merge(pn.info, odour.scenes.agg, all.x = TRUE, all.y = FALSE)
pn.meta = neuprint_get_meta(pn.info$bodyid)
rownames(pn.info) = pn.info$bodyid
pn.info$side = "right"
pn.info[rownames(pn.info)%in%upns.left.fib.info$bodyId] = "left"
pn.meta = hemibrain_get_meta(as.character(pn.info$bodyid))
pn.info = cbind(pn.info, pn.meta[as.character(pn.info$bodyid),setdiff(colnames(pn.meta),colnames(pn.info))])
pn.info$n.cell_type = sapply(pn.info$type, function(x) nrow(subset(pn.info,type==x)))
pn.info$class = pn.info$PN_type
pn.info$neurotransmitter = ifelse(pn.info$tract=="mlALT","GABA","acetylcholine")
pn.info$glomerulus = gsub(" ","",pn.info$glomerulus)
pn.info$glomerulus[pn.info$class=='mPN'] = "mPN"
pn.info=pn.info[,!colnames(pn.info)%in%c( "X", "PN_type", "fafb_type", "PN_type_revised", "best_fafb_match")]
pn.info$putative.classic.transmitter[is.na(pn.info$putative.classic.transmitter)] = "unknown"
pn.info$putative.classic.transmitter = hemibrainr::standardise(pn.info$putative.classic.transmitter)

# ALRN info
rn.info = read.csv("/Users/GD/LMBD/Papers/hemibrain_olf_data/FIB_RNs.csv")
rownames(rn.info) = rn.info$bodyid
rn.info$glomerulus = gsub(" ","",rn.info$glomerulus)
rn.meta = hemibrain_get_meta(as.character(rn.info$bodyid))
rn.info = cbind(rn.info, rn.meta[as.character(rn.info$bodyid),setdiff(colnames(rn.meta),colnames(rn.info))])
rn.info$class = ifelse(rn.info$modality=="olfactory","ORN","THRN")
rn.info$side = rn.info$soma_side
rn.info=rn.info[,!colnames(rn.info)%in%c( "X", "PN_type", "fafb_type", "PN_type_revised", "best_fafb_match", "soma_side")]
rn.info = rn.info[,!grepl("dend\\.|pd\\.|segregation",colnames(rn.info))]
rn.info$putative.classic.transmitter = "acetylcholine"

## MBON information
mbon.info = neuprint_search(".*MBON.*")
mb_ann = gsheet_manipulation(FUN = googlesheets4::read_sheet,
                                          ss = "1NROq0fw-WXGgu_mcWWVW9-l3YOdWHk-nESkI23fz-bo",
                                          sheet = "MBON",
                                          return = TRUE)#read.csv("data-raw/annotations/mb_annotations.csv")
mbon.info = cbind(mbon.info,mb_ann[match(mbon.info$bodyid,mb_ann$bodyid),])
mbon.info[is.na(mbon.info)] = "unknown"
# mbon.info = subset(mbon.info, compartments!=""&!is.na(compartments))
mbon.info$class = "MBON"
mbon.info$compartments = gsub("+\\).*","",mbon.info$name)
mbon.info$compartments = gsub("+.*\\(","",mbon.info$compartments)
mbon.info$compartments = gsub("_.*","",mbon.info$compartments)
rownames(mbon.info) = mbon.info$bodyid
mbon.meta = hemibrain_get_meta(as.character(mbon.info$bodyid))
mbon.info = cbind(mbon.info, mbon.meta[as.character(mbon.info$bodyid),setdiff(colnames(mbon.meta),colnames(mbon.info))])
mbon.info$naive.valence = "unknown"
mbon.info$naive.valence[mbon.info$valence=="app"] = "aversive"
mbon.info$naive.valence[mbon.info$valence=="aver"] = "appetitive"
mbon.info$neurotransmitter = hemibrainr::standardise(mbon.info$neurotransmitter)

### olfactory TON information
ton.info = lhns::hemibrain_tons
ton.info$class[ton.info$class=="TON"] = "TOON"
ton.info$class[ton.info$bodyid%in%cent.ids] = "LHCENT"
ton.meta = hemibrain_get_meta(as.character(ton.info$bodyid))
ton.meta$side = "right"
ton.meta$putative.classic.transmitter[ton.meta$putative.classic.transmitter=="acetlcholine"] = "acetylcholine"
ton.meta$putative.classic.transmitter[ton.meta$putative.classic.transmitter=="acetylcoline"] = "acetylcholine"
ton.meta$putative.classic.transmitter[ton.meta$putative.classic.transmitter=="NA"] = "unknown"
ton.meta$putative.classic.transmitter[is.na(ton.meta$putative.classic.transmitter)]= "unknown"
ton.info = cbind(ton.info, ton.meta[as.character(ton.info$bodyid),setdiff(colnames(ton.meta),colnames(ton.info))])
ton.info$putative.classic.transmitter = hemibrainr::standardise(ton.info$putative.classic.transmitter)

## Visual projection neuron information
lc.info = neuprint_search("LC.*",field="type")
lc.info = hemibrain_get_meta(lc.info$bodyid)
lc.info$class = "LC"

# AL local neuron informatipn
al.sheet = gsheet_manipulation(FUN = googlesheets4::read_sheet,
                                            ss = "124eTYqQ8evTGm_z75V8jNVmfBI763_s4h1EAPVMiSvI",
                                            sheet = "AL",
                                            return = TRUE)
# alln.summary = read.csv("/Users/GD/LMBD/Papers/hemibrain_olf_data/LN_type_description.csv")
alln.summary = read.csv("/Users/GD/LMBD/Papers/hemibrain_olf_data/ALLN_type_group_description.csv")
alln.info = hemibrain_get_meta(alln.ids)
alln.info = alln.info[!is.na(alln.info$type),]
alln.info$class = "ALLN"
alln.info$class = al.sheet[match(alln.info$bodyid,al.sheet$bodyid),"manual_morph_group"]
alln.info$class = paste0("ALLN_",alln.info$class)
alln.info$side = "right"
rownames(alln.info) = alln.info$bodyid
alln.info$anatomy.group = alln.summary$classification[match(alln.info$bodyid,alln.summary$bodyid)]
# alln.info$anatomy.group = gsub("_$","",alln.info$anatomy.group)
alln.info$group = gsub(" |  ","/",alln.summary$grouped[match(alln.info$bodyid,alln.summary$bodyid)])
alln.info$notes = paste0("Type: ",alln.summary$type.description[match(alln.info$bodyid,alln.summary$bodyid)],". Group: ", alln.summary$type.description[match(alln.info$bodyid,alln.summary$bodyid)])

# See what the status of our FAB matches is. We need presynapses for predictions.
# a = alln.info[,c("bodyid", "cellBodyFiber","type", "ito_lee_hemilineage", "FAFB.match", "FAFB.match.quality", "class")]
# a = subset(a, FAFB.match!='none' & !is.na(FAFB.match))
# a$catmaid_name = catmaid_get_neuronnames(a$FAFB.match)
# a$pre = sapply(a$FAFB.match,function(s) length(unique(catmaid_get_connector_table(s, direction = "outgoing")$connector_id)))
# catmaid_set_annotations_for_skeletons(a$FAFB.match, annotations = "matched_ALLNs")

# KC information
kc.info = neuprint_search("^KC.*",field="type")
kc.info = hemibrain_get_meta(kc.info$bodyid)
kc.info$class = "KC"
kc.info$putative.classic.transmitter = "acetylcholine"

# CENT information
cent.info = hemibrain_get_meta(unique(hemibrainr::cent.ids))
cent.info$class = "LHCENT"

# Dn information
dn.info = hemibrain_get_meta(unique(hemibrainr::dn.ids))
dn.info$class = "DN"
dn.info = dn.info[,!grepl("dend\\.|pd\\.|segregation|axon\\.",colnames(dn.info))]

# chatacter when needs
char_df <- function(df, cols = c("bodyid","match","FAFB.match","skid","flywire_id")){
  cols = intersect(cols,colnames(df))
  for(col in cols){
    df[[col]] = as.character(df[[col]])
  }
  df
}
rn.info = char_df(rn.info)
pn.info = char_df(pn.info)
mbon.info = char_df(mbon.info)
ton.info = char_df(ton.info)
lc.info = char_df(lc.info)
alln.info = char_df(alln.info)
kc.info = char_df(kc.info)
cent.info = char_df(cent.info)
dn.info = char_df(dn.info)
hemibrain_glomeruli_summary = char_df(hemibrain_glomeruli_summary)

# Save information
usethis::use_data(rn.info, overwrite = TRUE)
usethis::use_data(pn.info, overwrite = TRUE)
usethis::use_data(mbon.info, overwrite = TRUE)
usethis::use_data(ton.info, overwrite = TRUE)
usethis::use_data(lc.info, overwrite = TRUE)
usethis::use_data(alln.info, overwrite = TRUE)
usethis::use_data(kc.info, overwrite = TRUE)
usethis::use_data(cent.info, overwrite = TRUE)
usethis::use_data(dn.info, overwrite = TRUE)
usethis::use_data(hemibrain_glomeruli_summary, overwrite = TRUE)
usethis::use_data(odour_scenes, overwrite = TRUE)

### Supplementary data for Schlegel and Bates 2021:
supp.cols = c("bodyid", "pre", "post", "upstream", "downstream",
              "status", "name",  "voxels", "soma", "side",
              "connectivity.type", "cell_type", "group","anatomy.group", "class", "cellBodyFiber",
              # "ito_lee_lineage", "ito_lee_hemilineage", "hartenstein_lineage", "hartenstein_hemilineage",
              "putative.classic.transmitter", "putative.other.transmitter", "glomerulus", "presyn.glom",
              "FAFB.match", "FAFB.match.quality", "layer", "ct.layer",
              "axon_outputs", "dend_outputs",
              "axon_inputs", "dend_inputs", "total_length", "cable_length.um", "cable_length.glom.um","axon_length", "dend_length",
              "pd_length", "segregation_index", "notes")
toon.supp = ton.info[,colnames(ton.info)%in%supp.cols]
order = match(supp.cols,colnames(toon.supp))
order = order[!is.na(order)]
toon.supp = toon.supp[,order]
alpn.supp = pn.info[,colnames(pn.info)%in%supp.cols]
order = match(supp.cols,colnames(alpn.supp))
order = order[!is.na(order)]
alpn.supp = alpn.supp[,order]
alrn.supp = rn.info[,colnames(rn.info)%in%supp.cols]
order =match(supp.cols,colnames(alrn.supp))
order = order[!is.na(order)]
alrn.supp = alrn.supp[,order]
alln.supp = alln.info[,colnames(alln.info)%in%supp.cols]
order =match(supp.cols,colnames(alln.supp))
order = order[!is.na(order)]
alln.supp = alln.supp[,order]
dn.supp = dn.info[,colnames(dn.info)%in%supp.cols]
order =match(supp.cols,colnames(dn.supp))
order = order[!is.na(order)]
dn.supp = dn.supp[,order]
# clean data a little
wipe = c("FAFB.match", "FAFB.match.quality","putative.classic.transmitter","putative.other.transmitter")
for(w in wipe){
  alpn.supp[[w]] = NULL
  alrn.supp[[w]] = NULL
  alln.supp[[w]] = NULL
  dn.supp[[w]] = NULL
}
# Write supp data
write.csv(toon.supp, file = "/Users/GD/LMBD/Papers/hemibrain_olf_data/Schlegel2021_supp/hemibrain_TOON_meta.csv", row.names = FALSE)
write.csv(alpn.supp, file = "/Users/GD/LMBD/Papers/hemibrain_olf_data/Schlegel2021_supp/hemibrain_ALPN_meta.csv", row.names = FALSE)
write.csv(alrn.supp, file = "/Users/GD/LMBD/Papers/hemibrain_olf_data/Schlegel2021_supp/hemibrain_ALRN_meta.csv", row.names = FALSE)
write.csv(alln.supp, file = "/Users/GD/LMBD/Papers/hemibrain_olf_data/Schlegel2021_supp/hemibrain_ALLN_meta.csv", row.names = FALSE)
write.csv(dn.supp, file = "/Users/GD/LMBD/Papers/hemibrain_olf_data/Schlegel2021_supp/hemibrain_DN_meta.csv", row.names = FALSE)
write.csv(hemibrain_glomeruli_summary, file = "/Users/GD/LMBD/Papers/hemibrain_olf_data/Schlegel2021_supp/hemibrain_olfactory_information.csv", row.names = FALSE)
hemibrain.roots = hemibrain_somas[,c("bodyid","X","Y","Z")]
write.csv(hemibrain.roots, file = "/Users/GD/LMBD/Papers/hemibrain_olf_data/Schlegel2021_supp/hemibrain_root_points.csv", row.names = FALSE)
hemibrain.splitpoints = c("bodyid", "position", "point", "X", "Y", "Z")
write.csv(hemibrain.splitpoints, file = "/Users/GD/LMBD/Papers/hemibrain_olf_data/Schlegel2021_supp/hemibrain_compartment_startpoints.csv", row.names = FALSE)
