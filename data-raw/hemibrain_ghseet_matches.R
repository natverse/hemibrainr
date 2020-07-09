###### Manipulate GoogleSheets for recording matches #######

## Code in the package lhns was first used to create these google sheets
### And make initial matches.

# The google sheet database:
# https://docs.google.com/spreadsheets/d/1OSlDtnR3B1LiB5cwI5x5Ql6LkZd8JOS5bBr-HTi0pOw/edit#gid=0
selected_file = "1OSlDtnR3B1LiB5cwI5x5Ql6LkZd8JOS5bBr-HTi0pOw"

# Add all hemibrain neurons to sheet
all.ids = hemibrain_neuron_bodyids()
hemibrain_matching_add(ids = all.ids, sheet = "hemibrain")

# Transfer matches that have already been made
hemibrain_matching_transfers()

# Add uPN information
upn.hemi.info = subset(pn.info, class == "uPN")
upn.faf.info  = subset(lhns::pn.fafb, class == "uPN")
match.pns = merge(upn.hemi.info,upn.fafb.info, all = TRUE)
