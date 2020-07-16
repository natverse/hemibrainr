###### Manipulate GoogleSheets for recording matches #######

## Code in the package lhns was first used to create these google sheets
### And make initial matches.

# The google sheet database:
# https://docs.google.com/spreadsheets/d/1OSlDtnR3B1LiB5cwI5x5Ql6LkZd8JOS5bBr-HTi0pOw/edit#gid=0
selected_file = "1OSlDtnR3B1LiB5cwI5x5Ql6LkZd8JOS5bBr-HTi0pOw"

# Add all hemibrain neurons to sheet
hemibrain_matching_rewrite()
fafb_matching_rewrite()
hemibrain_matching_transfers()

# Remove confounding annotations in segmentation instances
fafb_lineage_complete_wipe(server = "v14seg-Li-190411.0")

# Get hemibrain information into v14 CATMAID for matches
matches = hemibrain_matches()
matches = subset(matches, match.quality %in% c("good","medium","poor") & dataset == "FAFB")
skds = rownames(matches)
fafb_hemibrain_annotate(skds)

# Flag 5 neurons from each cell body fiber
## To match to at least tract level
gs = hemibrain_match_sheet(sheet = "hemibrain", selected_file = selected_file)
gs.undone = subset(gs, is.na(gs$FAFB.match))
gs.done = subset(gs, is.na(gs$FAFB.match))
cbfs = setdiff(hemibrain_hemilineages$cellBodyFiber, unique(gs.done$cellBodyFiber))
chosen = c()
for(cbf in cbfs){
  choose = subset(gs.undone,cellBodyFiber==cbf)$bodyid
  choose = choose[1:ifelse(length(choose)<5,length(choose),5)]
  chosen = c(chosen,choose)
}
gs[chosen,"User"] = "ASB4"
write_matches(gs=gs,
              ids = chosen,
              column = "User")

# Organise FAFB hemilineage annotations
a <- catmaid::catmaid_get_annotationlist()
volker.old <- subset(a$annotations, grepl("Volker_|volker_", name))
volker <- subset(a$annotations, grepl("Hartenstein", name))
itolee <- subset(a$annotations, grepl("Ito", name))
volker.lineages <-  subset(a$annotations, grepl("Hartenstein_Lineage: ", name))
volker.hemilineages <-  subset(a$annotations, grepl("Hartenstein_Hemilineage: ", name))
itolee.lineages <-  subset(a$annotations, grepl("ItoLee_Lineage: ", name))
itolee.hemilineages <-  subset(a$annotations, grepl("ItoLee_Hemilineage: ", name))
putative.transmitter <-  subset(a$annotations, grepl("Putative_Neurotransmitter: ", name))
meta = c("ItoLee_Lineage",
         "ItoLee_Hemilineage",
         "Hartenstein_Lineage",
         "Hartenstein_Hemilineage")
catmaid::catmaid_set_meta_annotations(meta_annotations = "ItoLee_Lineage", annotations = itolee.lineages$id)
catmaid::catmaid_set_meta_annotations(meta_annotations = "ItoLee_Hemilineage", annotations = itolee.hemilineages$id)
catmaid::catmaid_set_meta_annotations(meta_annotations = "Hartenstein_Lineage", annotations = volker.lineages$id)
catmaid::catmaid_set_meta_annotations(meta_annotations = "Hartenstein_Hemilineage", annotations = volker.hemilineages$id)
catmaid::catmaid_set_meta_annotations(meta_annotations = "Putative_Neurotransmitter", annotations = putative.transmitter$id)
catmaid::catmaid_set_meta_annotations(meta_annotations = "Lineage_annotations", annotations = c(itolee.lineages$id,
                                                                                       itolee.hemilineages$id,
                                                                                       volker.lineages$id,
                                                                                       volker.hemilineages$id))
for(ann in itolee.lineages$name){
  skids <- catmaid::catmaid_skids(paste0("annotation:", ann))
  catmaid::catmaid_set_annotations_for_skeletons(skids = skids, annotations = "Lineage_annotated")
}
wrongs  = c()
for(o in volker.old$name){
  w = catmaid::catmaid_skids(o)
  wrongs = c(wrongs, w)
}
catmaid::catmaid_remove_annotations_for_skeletons(unique(wrongs), volker.old$id, force = TRUE)

# Correct mis-assignment in FAFB. Example:
# ## Wrong FAFB Ito_Lee Hemilineage: CREa2
# ## Correct one: CREa2_ventral
# ## Hemibrain CBF: ADM05
# fafb_set_hemilineage(find = "ItoLee_Hemilineage: CREa2",
#                      putative = FALSE,
#                      ItoLee_Hemilineage = "CREa2_ventral")



