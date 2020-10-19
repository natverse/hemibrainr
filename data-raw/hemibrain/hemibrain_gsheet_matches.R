###### Manipulate GoogleSheets for recording matches #######

## Code in the package lhns was first used to create these Google sheets
### And make initial matches.

# The Google sheet database:
# https://docs.google.com/spreadsheets/d/1OSlDtnR3B1LiB5cwI5x5Ql6LkZd8JOS5bBr-HTi0pOw/edit#gid=0
selected_file = "1OSlDtnR3B1LiB5cwI5x5Ql6LkZd8JOS5bBr-HTi0pOw"

# Unify annotations
hemibrain_hemilineages = read.csv("data-raw/annotations/hemibrain_hemilineages_cbf.csv")
fafb.meta = elmr::fafb_get_meta("Lineage_annotated")
hemibrain_hemilineages[hemibrain_hemilineages$ItoLee_Hemilineage %in% fafb.meta$ItoLee_Hemilineage,"FAFB"] = TRUE
hemibrain_hemilineages[!hemibrain_hemilineages$ItoLee_Hemilineage %in% fafb.meta$ItoLee_Hemilineage,"FAFB"] = FALSE
hemibrain_hemilineages = write.csv(x=hemibrain_hemilineages, file = "data-raw/annotations/hemibrain_hemilineages_cbf.csv",row.names = FALSE)
hls = unique(hemibrain_hemilineages$ItoLee_Hemilineage)
while(length(hls)){
  for(hl in hls){
    print(hl)
    ann = paste0("ItoLee_Hemilineage: ", hl)
    t = tryCatch(fafb_set_hemilineage(find = ann,
                         putative = FALSE,
                         ItoLee_Hemilineage = hl),
             error = function(e) "error")
    if(!is.null(t)){
      message(t)
      next
    }
    hls = setdiff(hls,hl)
  }
}

# Add all hemibrain neurons to sheet
hemibrain_matching_rewrite()
fafb_matching_rewrite()
hemibrain_matching_transfers()
flywire_matching_rewrite()

# Remove confounding annotations in segmentation instances
fafb_lineage_complete_wipe(server = "v14seg-Li-190411.0")

# Get hemibrain information into v14 CATMAID for matches
matches = hemibrain_matches()
matches = subset(matches, quality %in% c("good","medium") & dataset == "FAFB")
skds = as.character(rownames(matches))
fafb_hemibrain_annotate(skds, flywire = TRUE)

# Flag 5 neurons from each cell body fiber
## To match to at least tract level
gs = hemibrain_match_sheet(sheet = "hemibrain", selected_file = selected_file)
gs.undone = subset(gs, !gs$FAFB.match.quality %in% c("good","medium"))
gs.done = subset(gs, gs$FAFB.match.quality %in% c("good","medium"))
cbfs = setdiff(hemibrainr::hemibrain_hemilineages$cellBodyFiber, unique(gs.done$cellBodyFiber))
cbfs = sort(cbfs[!grepl("miss",cbfs)])
chosen = c()
for(cbf in cbfs){
  choose = subset(gs.undone,cellBodyFiber==cbf)$bodyid
  if(!length(choose)){
    next
  }
  choose = choose[sample(length(choose),ifelse(length(choose)>5,5,length(choose)))]
  chosen = unique(c(chosen,choose))
}
gs[chosen%in%gs$bodyid,"User"] = "ASB7"
write_matches(gs=gs,
              ids = intersect(chosen,gs$bodyid),
              column = "User")

# Predict left-side seed planes in FAFB
seeds = gsheet_manipulation(FUN = googlesheets4::read_sheet,
                         ss = "1HI8RZJhV8aWC6hw5T0qh2__9D8V0zKdtxx2kkDsuI8Y",
                         sheet = "seeds",
                         return = TRUE)
seeds$left_guess_plane = sapply(seeds$right_seed_plane,elmr::mirror_fafb_url)
googlesheets4::write_sheet(seeds,
                           ss = "1HI8RZJhV8aWC6hw5T0qh2__9D8V0zKdtxx2kkDsuI8Y",
                           sheet = "seeds")

# Organise FAFB hemilineage annotations
a <- catmaid::catmaid_get_annotationlist()
volker.old <- subset(a$annotations, grepl("Volker_|volker_", name))
volker <- subset(a$annotations, grepl("Hartenstein_", name))
itolee <- subset(a$annotations, grepl("ItoLee_", name))
volker.lineages <-  subset(a$annotations, grepl("Hartenstein_Lineage: ", name))
volker.hemilineages <-  subset(a$annotations, grepl("Hartenstein_Hemilineage: ", name))
itolee.lineages <-  subset(a$annotations, grepl("ItoLee_Lineage: ", name))
itolee.hemilineages <-  subset(a$annotations, grepl("ItoLee_Hemilineage: ", name))
putative.transmitter <-  subset(a$annotations, grepl("known transmitter:|putative transmitter:", name))
catmaid::catmaid_set_meta_annotations(meta_annotations = "ItoLee_Lineage", annotations = itolee.lineages$id)
catmaid::catmaid_set_meta_annotations(meta_annotations = "ItoLee_Hemilineage", annotations = itolee.hemilineages$id)
catmaid::catmaid_set_meta_annotations(meta_annotations = "Hartenstein_Lineage", annotations = volker.lineages$id)
catmaid::catmaid_set_meta_annotations(meta_annotations = "Hartenstein_Hemilineage", annotations = volker.hemilineages$id)
catmaid::catmaid_set_meta_annotations(meta_annotations = "transmitter", annotations = putative.transmitter$id)
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

# Set transmitters for known matches in FAFB
lm = lm_matches()
lm.select = subset(lm, id %in% c(names(subset(lhns::lh.mcfo,classic.transmitter%in%c("acetylcholine","GABA","glutamate"))),
                          names(subset(lhns::lh.splits.dps,classic.transmitter%in%c("acetylcholine","GABA","glutamate"))))
                          & quality %in% c("good","medium"))
hemi = unique(lm.select$match)
for(h in hemi){
  match = lm[h,"match"]
  t = lhns::lh.mcfo[match,"classic.transmitter"]
  if(is.issue(t)){
    t = lhns::lh.splits.dps[match,"classic.transmitter"]
  }
  if(is.issue(t)|t%in%c("none"," ","","NA","unknown")){
    next
  }
  tryCatch(fafb_set_transmitter(find = h, dataset = "hemibrain", putative = FALSE, transmitter = t),
           error = function(e) message("Error for: ",h," ", t))
}

# Correct mis-assignment in FAFB. Example:
# ## Wrong FAFB Ito_Lee Hemilineage: CREa2
# ## Correct one: CREa2_ventral
# ## Hemibrain CBF: ADM05
# fafb_set_hemilineage(find = "ItoLee_Hemilineage: LHd2",
#                      putative = FALSE,
#                      ItoLee_Hemilineage = "LHd2")



