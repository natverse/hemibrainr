#####################################################################################
################################ FIB-FAFB Matching ##################################
#####################################################################################

#' Match up neurons between the hemibrain and FAFB
#'
#' @description Match up neurons between the hemibrain and FAFB and save rhe result fo a \href{https://docs.google.com/spreadsheets/d/1OSlDtnR3B1LiB5cwI5x5Ql6LkZd8JOS5bBr-HTi0pOw/edit#gid=0}{Google Sheet} on the hemibrain Google Team Drive operated by
#' the FlyConnctome group at the University of Cambridge. Currently by default this function only considers lateral horn neurons. You must have access to the Team Drive in order to use this function. This function makes use of Google Filestream, which should be mounted
#' on your machine. Further, note that neurons are read from the FAFB CATMAID project, and you must have log-in details for this project recorded in your .renviron for this function to work.
#'
#' @param bodyids body IDs for hemibrain neurons present in the \href{https://docs.google.com/spreadsheets/d/1OSlDtnR3B1LiB5cwI5x5Ql6LkZd8JOS5bBr-HTi0pOw/edit#gid=0}{Google Sheet}, for which the user will attempt to make a match if one has not been made already.
#' @param hemibrain.fafb.nblast a FAFB (columns) - hemibrain (rows) normalised NBLAST matrix. By default this is read from the FlyConnectome Team Drive.
#' @param selected_file the Google Sheet databse to read and write from. For now, defaults to a \href{https://docs.google.com/spreadsheets/d/1OSlDtnR3B1LiB5cwI5x5Ql6LkZd8JOS5bBr-HTi0pOw/edit#gid=0}{Google Sheet for lateral horn neurons}. No other databases have been prepared.
#' @param batch_size the number of FAFB top matches to read from CATMAID in one go.
#'
#' @details Currently, the \href{https://docs.google.com/spreadsheets/d/1OSlDtnR3B1LiB5cwI5x5Ql6LkZd8JOS5bBr-HTi0pOw/edit#gid=0}{Google Sheet} is set up with  limited number of users, each of whom have been assigned a number of neurons to match up. In order
#' to add yourself as a user, simply open this Google Sheet in your browser and add your initials to neurons of your choosing on the rightmost column 'Users'. Once a match is recorded, the user selects a quality for that match. There can be
#' no match (n), a poor match (p) an okay match (o) or an exact match (e). As a rule of thumb, a poor match could be a neuron from a very similar same cell tye or a highly untraced neuron that may be the correct cell type. An okay match should be a neuron
#' that looks to be from the same morphological cell type but there may be some discrepancies in its arbour. An exact match is a neuron that corresponds well between FAFB and the hemibrain data.
#'
#' @examples
#' \donttest{
#' \dontrun{
#'
#' # install package to bridge neurons between FAFB14 and hemibrain space
#' if (!requireNamespace("remotes")) install.packages("remotes")
#' remotes::install_github('natverse/nat.jrcbrains')
#' nat.jrcbrains::download_saalfeldlab_registrations()
#'
#' # Load precomputed NBLAST from the flyconnectome Team Drive
#' load("/Volumes/GoogleDrive/Shared\ drives/flyconnectome/fafbpipeline/fafb (1).fib.twigs5.crossnblast.rda")
#'  # load("/Volumes/GoogleDrive/Shared\ drives/flyconnectome/fafbpipeline/fafb.fib.twigs5.crossnblast.rda") ## Or this one
#'
#' # Match!
#' hemibrain_FAFB_matching(hemibrain.fafb.nblast = fafb.fib.twigs5.crossnblast)
#' }}
#' @rdname hemibrain_FAFB_matching
#' @export
#' @seealso \code{\link{hemibrain_adjust_saved_split}}
hemibrain_FAFB_matching <- function(bodyids = NULL,
                         hemibrain.fafb.nblast = NULL,
                         selected_file = "1OSlDtnR3B1LiB5cwI5x5Ql6LkZd8JOS5bBr-HTi0pOw",
                         batch_size = 10){
  requireNamespace("nat.jrcbrains")
  requireNamespace("elmr")
  #
  # Motivate!
  nat::nopen3d()
  plot_inspirobot()
  unsaved = c()
  message("
          #######################Colours#######################
          black = hemibrain neuron,
          red = potential FAFB matches based on NBLAST score,
          green = a chosen neuron during scanning,
          blue = your selected match,
          #######################Colours#######################
          ")
  ## Get NBLAST
  if(is.null(hemibrain.fafb.nblast)){
    message("Loading FAFB-FIB NBLAST from flyconnectome Google Team Drive using Google Filestream: ")
    load("/Volumes/GoogleDrive/Shared\ drives/flyconnectome/fafbpipeline/fafb (1).fib.twigs5.crossnblast.rda")
    # load("/Volumes/GoogleDrive/Shared\ drives/flyconnectome/fafbpipeline/fafb.fib.twigs5.crossnblast.rda") ## Or this one
    hemibrain.fafb.nblast = t(fafb.fib.twigs5.crossnblast)
    rm("fafb.fib.twigs5.crossnblast")
  }
  # Read the Google Sheet
  gs = gsheet_manipulation(FUN = googlesheets4::read_sheet,
                           ss = selected_file,
                           sheet = "lhns",
                           return = TRUE)
  gs$bodyid = correct_id(gs$bodyid)
  rownames(gs) = gs$bodyid
  if(is.null(bodyids)){
    bodyids = gs$bodyid
  }else{
    bodyids = intersect(bodyids,gs$bodyid)
  }
  # Read neuron meta data
  meta = neuprintr::neuprint_get_meta(bodyids)
  meta = meta[order(meta$type),]
  meta = subset(meta, meta$bodyid%in%bodyids)
  db = tryCatch(hemibrain_neurons(), error = function(e) NULL)
  # How much is done?
  done = subset(gs, !is.na(gs$FAFB.match))
  message("Neurons matches: ", nrow(done))
  print(table(gs$FAFB.match.quality))
  # Choose user
  message("Users: ", paste(unique(gs$User),collapse = " "))
  initials = must_be("Choose a user : ", answers = unique(gs$User))
  rgl::bg3d("white")
  # Make matches!
  for(n in meta$bodyid){
    # Get bodyid
    n = as.character(n)
    # Remove neurons with matches
    donotdo = subset(gs, !is.na(gs$FAFB.match) | User != initials)
    if(n%in%donotdo$bodyid){
      next
    }
    # Plot brain
    clear3d()
    plot3d(elmr::FAFB14,alpha = 0.1, col ="grey")
    # Read hemibrain neuron
    if(is.null(db)){
      lhn  = neuprintr::neuprint_read_neurons(n)
    }else{
      lhn = tryCatch(db[as.character(n)], error = function(e) NULL)
    }
    if(is.null(lhn)){
      lhn  = neuprintr::neuprint_read_neurons(n)
    }
    # Transform hemibrain neuron to FAFB space
    lhn = scale_neurons.neuronlist(lhn, scaling = (8/1000))
    lhn = suppressWarnings(nat.templatebrains::xform_brain(lhn, reference = "FAFB14", sample = "JRCFIB2018F"))
    message("Hemibrain body ID: ", lhn[n,"bodyid"])
    message("Hemibrain-assigned cell type :",lhn[n,"type"])
    # Read top 10 FAFB matches
    message(sprintf("Reading the top %s FAFB hits",batch_size))
    r = sort(hemibrain.fafb.nblast[n,],decreasing = TRUE)
    plot3d(lhn[n], lwd = 2, soma = 500, col = "black")
    fafb = catmaid::read.neurons.catmaid(names(r)[1:batch_size])
    sel = c("go","for","it")
    k = 1
    j = batch_size
    # Cycle theough potential matches
    while(length(sel)>1){
      if(j>10){
        fafb = union(fafb, read.neurons.catmaid(names(r)[(k+1):j]))
      }
      sel = nat::nlscan(fafb[names(r)[1:j]], col = "red", soma = TRUE)
      if(length(sel)>1){
        message("You selected more than one neuron, please select again ... ")
      }
      prog = hemibrain_choice(sprintf("You selected %s neurons. Are you happy with that? ",length(sel)))
      if(!prog){
        sel = c("go","for","it")
        prog = hemibrain_choice("Do you want to read more neurons from CATMAID? ")
        if(prog){
          k = j
          j = j + batch_size
        }
      }
    }
    # Assign match and its quality
    gs[n,"FAFB.match"] = ifelse(length(sel)==0,'none',sel)
    if(length(sel)){
      rgl::plot3d(fafb[sel],col="blue",lwd=2,soma=TRUE)
      quality = must_be("What is the quality of this match? exact(e)/okay(o)/poor(p) ", answers = c("e","o","p"))
    }else{
      quality = "n"
    }
    gs[n,"FAFB.match.quality"] = quality
    gs = gs[!duplicated(gs$bodyid),]
    unsaved = c(unsaved, n)
    message(length(unsaved), " unsaved matches")
    print(knitr::kable(gs[unsaved,c("bodyid","type","FAFB.match","FAFB.match.quality")]))
    p = must_be("Continue (enter) or save (s)? ", answers = c("","s"))
    if(p=="s"){
      # Write!
      plot_inspirobot()
      write_matches(gs=gs,
                 bodyids = unsaved,
                 column = "FAFB.match")
      write_matches(gs=gs,
                 bodyids = unsaved,
                 column = "FAFB.match.quality")
      gs = gsheet_manipulation(FUN = googlesheets4::read_sheet,
                                            ss = selected_file,
                                            sheet = "lhns",
                                            return = TRUE)
      gs$bodyid = correct_id(gs$bodyid)
      rownames(gs) = gs$bodyid
      unsaved = c()
      rgl::bg3d("white")
    }
  }
}

# hidden
## Google sheet write function
write_matches <- function(gs,
                       bodyids = NULL,
                       selected_file = "1OSlDtnR3B1LiB5cwI5x5Ql6LkZd8JOS5bBr-HTi0pOw",
                       column = colnames(gs)){
  if((!identical(colnames(gs),column) & length(column) > 1)|(sum(column%in%colnames(gs))<1)){
    stop("Column must be one column of the google sheet, or all the columns")
  }
  rows = (1:nrow(gs))+1
  rownames(gs) = rows
  if(!is.null(bodyids)){
    gs = subset(gs, gs$bodyid %in% bodyids)
    message("Updating ", nrow(gs), " entries")
    rows = rownames(gs)
  }
  for(r in rows){
    if(length(column)==1){
      letter = LETTERS[match(column,colnames(gs))]
      range = paste0(letter,r)
    }else{
      range = paste0("A",r,":",LETTERS[ncol(gs)],r)
    }
    gsheet_manipulation(FUN = googlesheets4::range_write,
                                     ss = selected_file,
                                     range = range,
                                     data = as.data.frame(gs[as.character(r),column], stringsAsFactors = FALSE),
                                     sheet = "lhns",
                                     col_names = FALSE)
  }
}

# hidden
correct_id <-function(v){
  gsub(" ","",v)
}
