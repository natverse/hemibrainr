## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----save.matches, eval = FALSE-----------------------------------------------
#  # Load package
#  library(hemibrainr)
#  # See matches
#  View(hemibrain_matched)

## ----fresh.matches, eval = FALSE----------------------------------------------
#  matches = hemibrain_matches() # You will be asked to log-in through a Google-enabled email address.
#  View(matches())

## ----ready.made, eval = FALSE-------------------------------------------------
#  made.matches = read.csv("my_matches.csv") # Must have the named columns: bodyid, skid, quality
#  hemibrain_matches(df = made.matches, direction = "both") # direction controls which tabs matches get written to

## ----add, eval = FALSE--------------------------------------------------------
#  # Add a mising FAFB projection neuron, so we can match it later:
#  hemibrain_matching_add(ids = "16", sheet = "fafb", User = "ASB")
#  ## the sheet argument specifies the worksheet or 'tab' on the Google sheet we want to add to

## ----interactive.pipeline, eval = FALSE---------------------------------------
#  # install package to bridge neurons between FAFB14 and hemibrain space
#  if (!requireNamespace("remotes")) install.packages("remotes")
#  remotes::install_github('natverse/nat.jrcbrains')
#  nat.jrcbrains::download_saalfeldlab_registrations()
#  
#  # Match hemibrain neurons!
#  hemibrain_matching() # Automatically, you can choose a User ID and you are given neurons that have this ID in the User column on the Google Sheet.
#  hemibrain_matching(ids=c("674108632","739256609")) # Otherwise you can select specific IDs
#  hemibrain_matching(ids=c("674108632","739256609"), overwrite = TRUE) # If a match has already been made you can overwrite it
#  # Otherwise neurons that have already been given a match will not be shown in the pipeline.
#  
#  # Match FAFB neurons!
#  fafb_matching()
#  fafb_matching(ids = "16") # Specify IDs
#  fafb_matching(ids = "16", overwrite = TRUE) # Overwrite

