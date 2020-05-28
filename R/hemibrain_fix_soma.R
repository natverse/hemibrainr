

### Core Function

#' @export
#' @rdname hemibrain_adjust_saved_split
##### lets tear this function apart
hemibrain_adjust_saved_somas <- function(bodyids = NULL,
                                         brain = hemibrainr::hemibrain.surf,
                                         selected_file = "1YjkVjokXL4p4Q6BR-rGGGKWecXU370D1YMc1mgUYr8E",
                                         db = NULL,
                                         batch_size = 5){
  ### Get GoogleSheet Database
  gs = gsheet_manipulation(FUN = googlesheets4::read_sheet,
                           ss = selected_file,
                           sheet = "roots",
                           return = TRUE)
  # how would you like to go through neurons? one at a time, or based on cell body fibre
  # or Auto, if I ever add it...
  mode <- must_be(prompt = "Which mode do you want to use? single neurons (s) / Cell body fibre from neurons (n) / Indifical cbf (c)",answers = c("s","n","c"))
  # reset3d(brain=brain)
  if(mode=="s"){
    if(is.null(bodyids)){
      message("Please provide neuron bodyids as input, exiting pipeline")
      stop()
    }
    neurons = pipeline_read_neurons(batch = bodyids, db = NULL, clean = FALSE, motivate = FALSE)
    neurons = correct_singles(neurons, brain = brain)
    save_soma_to_gsheet(neurons = neurons, gs = gs, selected_file = selected_file)
    save_checked_soma_to_gsheet(neurons = neurons, gs = gs, selected_file = selected_file)
    message("Task updated! ")
  }
  if(mode=="n"){
    if(is.null(bodyids)){
      message("Please provide neuron bodyids as input, exiting pipeline")
      stop()
    }
    message("Collecting Cell fibre bodies from provided neurons")
    # get cell fibre body list from neurons
    cbfs = cbf_from_bodyid(bodyids)
    # how many are their, and list them
    message(c("The provided bodyids are from ", length(cbfs), " different cell fibre bodies"))
    if("unknown" %in% cbfs){
      message("This includes some neurons labeled with an 'unknown' cbf, these neurons will be checked individually")
    }
    # loop over each, printing which you are on
    for(c in cbfs){
      message(c("correcting cbf ", c))
      correct_cbf(cbf = c, gs = gs, selected_file = selected_file)
    }
  }
  if(mode=="c"){
    # input a cell fibre body
    cbfs = cbf_list()
    cbf = must_be(prompt = "Please input a Cell fibre body to correct",
                          answers = cbf_list())
    if( cbf == "no answer yet fool"){
      message("You must provide a cbf which exists in the hemibrain data. Use the cbf_list() function to check the list")
      stop()
    }
    # if yes, fix...
    message(c("Fixing Cell fibre body: ", cbf))
    correct_cbf(cbf = cbf, gs = gs, selected_file = selected_file)
  }
}

#### Correct single neurons

# hidden
#' @importFrom nat neuronlist nopen3d
#' @importFrom rgl clear3d points3d spheres3d
correct_singles <- function(neurons, brain = NULL){
  list = 1
  if (is.neuron(neurons)) {
    neurons = neuronlist(neurons)
    list = 0
  }
  nopen3d()
  correcting = TRUE
  while(correcting){
    for (n in neurons) {
      n.points =  nat::xyzmatrix(n)
      end_points = nat::endpoints(n)
      end_points = n.points[end_points,]
      clear3d()
      plot3d(brain, col = "grey70", alpha = 0.1)
      #
      points3d(n.points, lwd = 2, col = hemibrainr::hemibrain_bright_colors["yellow"], soma = FALSE)
      plot3d_somas(n)
      fix = hemibrain_choice(prompt = "Does the soma need fixing? Current possition in Green, if present yes/no ")
      if(isTRUE(fix)) {
        make.selection = TRUE
      }
      else {
        make.selection = FALSE
      }
      while(make.selection){
        angle = readline(prompt = "position view in order to best select soma (Enter to continue) ")
        message("Click and drag a box over an end point, in red, to select new soma location")
        clear3d()
        plot3d(brain, col = "grey70", alpha = 0.1)
        points3d(end_points, lwd = 2, col = hemibrainr::hemibrain_bright_colors["red"], soma = FALSE)
        points3d(n.points, lwd = 2, col = hemibrainr::hemibrain_bright_colors["yellow"], soma = FALSE)

        selection <- rgl::select3d()

        selected = selection(end_points)
        selected.point = end_points[selected,]
        if (length(selected.point) != 3) {
          selected.point = selected.point[1,]
        }
        clear3d()
        plot3d(brain, col = "grey70", alpha = 0.1)
        points3d(n.points, lwd = 2, col = hemibrainr::hemibrain_bright_colors["yellow"], soma = FALSE)
        spheres3d(selected.point, radius = 300, col = 'blue')
        plot3d_somas(n)
        make.selection = !hemibrain_choice(prompt = c("Happy with selection? (old soma in Green, new soma in Blue) yes/no "))

        # reroot neuron
        eps <- nat::endpoints(n)
        #
        ep.sel <- selection(nat::xyzmatrix(n)[eps,])
        ep.sel <- eps[ep.sel][1]
        soma.id <- n$d$PointNo[match(ep.sel, 1:nrow(n$d))]
        # create ourneuron as a graph, with the new origin point
        y <- nat::as.neuron(nat::as.ngraph(n$d), origin =26)
        # carryover labels and tags
        y <- hemibrain_carryover_labels(x=n,y=y)
        y <- hemibrain_carryover_tags(x=n,y=y)
        y$tags$soma = soma.id
        y$tags$soma.edit = TRUE
        y = hemibrain_neuron_class(y)
        neurons[[toString(n$bodyid)]] = y
      }
    }
    clear3d()
    plot3d(neurons, soma = TRUE, WithConnectors = FALSE)
    plot3d_somas(neurons)
    correcting = !hemibrain_choice(prompt = c("Final check, are you happy with the new soma possitions? yes/no "))
  }
  if (list == 0) {
    neurons = hemibrain_neuron_class(neurons)
  }
  neurons
}

### read/wrte functions

# Hey Nik, this function reads neurons for the pipline.
# hidden
pipeline_read_neurons <- function(batch,
                                  db = NULL,
                                  motivate = TRUE,
                                  clean = TRUE){
  readfail = FALSE
  if(!is.null(db)){
    message("Reading locally saved neurons ...")
    someneuronlist = tryCatch( db[as.character(batch)], error = function(e) NULL)
    if(is.null(someneuronlist)){
      message("Errors reading from given neuronlist, reading batch from neuPrint instead ...")
      readfail = TRUE
    }
  }
  if(is.null(db)|readfail){
    if(motivate){plot_inspirobot()}
    message("Reading and manipulating neurons from neuPrint ...")
    someneuronlist = hemibrain_read_neurons(x = as.character(batch),
                                            savedir = FALSE,
                                            remove.bad.synapses = FALSE,
                                            microns = FALSE,
                                            clean = FALSE,
                                            OmitFailures = TRUE)
  }
  if(clean){
    someneuronlist = hemibrain_clean_skeleton(someneuronlist, rval = "neuron")
  }
  someneuronlist = add_field_seq(someneuronlist,someneuronlist[,"bodyid"],field="bodyid")
  someneuronlist
}

# hidden
save_soma_to_gsheet = function(neurons = neurons,
                               gs = gs,
                               selected_file = selected_file){
  ### Save to Google Sheet
  update = prepare_update(someneuronlist=neurons,gs=gs,initials="flyconnectome",print=FALSE)
  update$position = update$soma
  update$point = "root"
  update = update[,!colnames(update)%in%c("checked", "user", "time","soma")]
  update = cbind(update, catmaid::soma(neurons))
  ### add neuron type here
  update$type = gs$type[as.numeric(row.names(update))]
  update$soma.edit = as.logical(update$soma.edit)
  update = subset(update, update$soma.edit)
  message("Your updates: ")
  print(knitr::kable(update))
  rows = as.numeric(rownames(update))
  for(r in rows){
    range = paste0("B",r,":H",r)
    up = update[as.character(r),intersect(colnames(gs),colnames(update))]
    if(sum(is.na(up))>1){
      message("Erroneous NAs generated for row ", r, ", dropping this update")
      print(up)
      next
    }
    gsheet_manipulation(FUN = googlesheets4::range_write,
                        ss = selected_file,
                        range = range,
                        data = up,
                        sheet = "roots",
                        col_names = FALSE)
  }
}

save_checked_soma_to_gsheet = function(neurons = neurons,
                                            gs = gs,
                                            selected_file = selected_file){
  update = prepare_update(someneuronlist=neurons,gs=gs,initials="flyconnectome",print=FALSE)
  rows = as.numeric(rownames(update))
  checked = as.data.frame(matrix(TRUE,nrow = length(neurons), ncol = 1))
  range = paste0("V", rows[1], ":V", rows[length(rows)])
  gsheet_manipulation(FUN = googlesheets4::range_write,
                      ss = selected_file,
                      range = range,
                      data = checked,
                      sheet = "roots",
                      col_names = FALSE)
}

#### Cell fibre body related functions

#' @importFrom rgl clear3d spheres3d legend3d
#' @importFrom RColorBrewer brewer.pal brewer.pal.info
correct_cbf = function(cbf = cbf, gs, selected_file){
  neurons = neurons_in_cbf(cbf)
  db = dbscan_neurons(neurons)

  ###
  # If only returns 0 in cluster...
  if (unique(db$cluster) == 0){
    message("All", length(db$cluster),  " somas have been labeled as noise")
    fix = hemibrain_choice(prompt = "Do you wish to correct each manually? yes|no")
    if(isTRUE(fix)){
      correct_singles(neurons)
    } else {
      message("You lazy person you... exiting pipeline")
      stop()
    }
  } else {
    # how many clusters and noise points are their?
    message("dbscan returns ", length(unique(db$cluster)), " cluster(s)")
    if (0 %in% db$cluster){
      message("This includes one cluster of ", sum(db$cluster == 0), " incorrect 'noise' somas")
    }
    ###
    # if only 1 cluster exists after removing 0's, is it correct?
    if (length(unique(db$cluster)) <= 2) {
      somas = soma_locations(neurons)
      # remove the 'noise' somas
      noise = somas[which(db$cluster == 0),]
      somas = somas[which(db$cluster == 1),]
      # plot and check if cluster is correct
      clear3d()
      plot3d(neurons,col = "grey70")
      spheres3d(somas[,c('X','Y','Z')],radius = 300,col = 'blue')
      spheres3d(noise[,c('X','Y','Z')],radius = 300, col = 'red')
      cluster_correct = hemibrain_choice(prompt = c("Has dbscan identified the correct soma cluster (in blue) yes|no "))
      ###
      # if yes, mark somas as correct and update sheet

      # if present correct 'noise' somas
      if(isTRUE(cluster_correct)){
        if (length(noise$PointNo) > 0){
          correct_noise = hemibrain_choice(prompt = c("Do you want to correct the incorrect somas? yes|no"))
          if(isTRUE(correct_noise)){
            neurons[which(db$cluster == 0)] = correct_singles(neurons[which(db$cluster == 0)])
          } else {
            message("fine, be lazy then")
          }
        }
      }
    } else {
      ### if multiple potential clusters were returned
      # get colour vector of distinct colours
      qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
      col = unlist(mapply(brewer.pal, qual_col_pals$maxcolors,
                          rownames(qual_col_pals)))
      # plot neurons with each cluster coloured, and select which is correct
      clear3d()
      plot3d(neurons,col = "grey70")
      count = 0
      for (c in unique(db$cluster)){
        count = count + 1
        cluster = neurons[which(db$cluster == c)]
        cluster = soma_locations(cluster)
        spheres3d(cluster[,c('X','Y','Z')],radius = 500,col = col[count])
      }
      legend3d("topright", legend = paste('Cluster', c(as.character(1:count))),
               pch = 10, col = col[1:count],
               cex=1, inset=c(0.02))
      clusters = as.numeric(must_be(prompt = "Which cluster is just the correctly identified somas? If none are correct, press n. Otheriwse, provide the number of the correct cluster. ",
                                    answers = c(as.character(1:count), "n")))
      if (clusters == "n"){
        # if no single cluster contains just somas
        message("Time to correct each soma individually, bad luck...")
        neurons = correct_singles(neurons)
      } else {
        # you have selected this cluster
        message(c(" You have chosen cluster ", as.character(clusters), " as the correct cluster of somas..."))
        message(c("sadly, now you have to fix the other ",as.character(sum(db$cluster != clusters))," somas individually."))
        #
        neurons[which(db$cluster == clusters)] = correct_singles(neurons[which(db$cluster == clusters)])
      }
    }
  }
  save_soma_to_gsheet(neurons = neurons, gs = gs, selected_file = selected_file)
  save_checked_soma_to_gsheet(neurons = neurons, gs = gs, selected_file = selected_file)
  message("Task updated! ")
}

# get a list of all annotated cell body fibres
#' @importFrom neuprintr neuprint_fetch_custom
cbf_list = function(){
  neu_cbf_list = neuprint_fetch_custom(
    cypher = "MATCH (n :`hemibrain_Neuron`) RETURN DISTINCT n.cellBodyFiber AS val")
  cbf_list = c()
  for (cbf in neu_cbf_list$data) {
    cbf_list = c(cbf_list,cbf[[1]][[1]])
  }
  cbf_list
}

# get neurons for a given cbf
neurons_in_cbf = function(cbf){
  ns = neuprintr::neuprint_search(search = cbf,field = "cellBodyFiber",all_segments = TRUE)
  neurons = pipeline_read_neurons(batch = ns$bodyid, db = NULL, clean = FALSE, motivate = FALSE)
  neurons
}

# get cbf for a given neuron/neuronlist
cbf_from_bodyid = function(bodyids){
  ### Get cell body fibre information
  meta = neuprint_get_meta(bodyids)
  cbfs = unique(meta$cellBodyFiber)
  if(NA %in% cbfs){
    cbfs = c(cbfs, "unknown")
  }
  cbfs
}

### Soma related miscelany
# get soma coords for neuron list
soma_locations = function(neurons){
  # get soma possitions, and identify missing
  somas = data.frame(matrix(NA,nrow = length(neurons),
                            ncol = length(colnames(neurons[[1]]$d))))
  colnames(somas) = colnames(neurons[[1]]$d)
  for (n in 1:length(neurons)){
    if (sum(is.na(neurons[[n]]$d[neurons[[n]]$soma,])) == 0) {
      if (!"bodyid" %in% names(neurons[[n]]$d[neurons[[n]]$soma,])) {
        add = neurons[[n]]$d[neurons[[n]]$soma,]
        add$bodyid = neurons[[n]]$bodyid
        somas[n,] = add
      } else {
        somas[n,] = neurons[[n]]$d[neurons[[n]]$soma,]
      }
    }
  }
  # remove some of the columns...
  somas = somas[,!colnames(somas)%in%c("Label", "W", "Parent")]
  somas
}

# Correct missing somas in neuron list
fix_missing_soma = function(neurons){
  somas = soma_locations(neurons)

  if (sum(is.na(somas$PointNo)) == 0) {
    message("All neurons have a soma! Yay!")
  } else {
    # message - x number of neurons have no soma, would you like to manually correct
    # them - all neurons need a soma for dbscan to work?
    message(" neurons without soma: ", sum(is.na(somas$PointNo)))
    check = hemibrain_choice(prompt = "Would you like to manualy place a soma? (yes|no)")
    if (check) {
      # if yes, single neuron pipeline for each
      # correct the neruons with no soma
      ind = which(is.na(somas[1]))
      for (i in ind)  {
        neurons[[i]] = correct_singles(neurons[[i]])
      }
    } else {
      message(" Fine, be lazy - we will use the StartPoint of the neurons with no soma")
      ind = which(is.na(somas[1]))
      for (i in ind)  {
        # use the neurons strt point
        neurons[[i]]$soma = neurons[[i]]$StartPoint
      }

    }
  }
  neurons
}

### dbscan clustering groups of somas

# run dbscan on a set of neurons - have this as a higher level function...
#' @importFrom dbscan dbscan
dbscan_neurons = function(neurons, eps = 1600, MinPts = 5){

  # check missing somas
  neurons = fix_missing_soma(neurons)

  # get soma points
  somas = soma_locations(neurons)

  # dbscan on soma points
  set.seed(123)
  db = dbscan(somas[,c('X','Y','Z')], eps = eps, MinPts = MinPts)
  # return dbscan
  db
}


#### Old functions from Alex not being used...

####
# Hey nik, this function does a lot of the leg work
# including saving to google sheet.
# hidden
cbf_check<-function(ids,
                    gs,
                    selected_file = "1YjkVjokXL4p4Q6BR-rGGGKWecXU370D1YMc1mgUYr8E",
                    cbf = "unknown batch",
                    db = NULL,
                    clean = FALSE,
                    motivate = FALSE,
                    brain = hemibrainr::hemibrain.surf){
  neurons = pipeline_read_neurons(batch = ids, db = db, clean = FALSE, motivate = motivate)
  continue = TRUE
  while(isTRUE(continue)){
    ### Plot
    reset3d(brain=brain)
    pnts = primary_neurite_cable(neurons, OmitFailures = TRUE)
    rgl::plot3d(pnts, lwd = 2, col = hemibrainr::hemibrain_bright_colors["purple"], soma = FALSE)
    rgl::plot3d(neurons, lwd = 2, col = hemibrainr::hemibrain_bright_colors["cerise"], soma = FALSE)
    plot3d_somas(neurons)
    ### Choose how to edit
    reset3d(brain=brain)
    neurons = correct_group(neurons, brain = brain)
    continue <- !hemibrain_choice(prompt = paste0("Are we finished with ", cbf,"? yes/no "))
  }
  ### Save to Google Sheet
  save_soma_to_gsheet(neurons = neurons, gs = gs, selected_file = selected_file)
  save_checked_soma_to_gsheet(neurons = neurons, gs = gs, selected_file = selected_file)

  message("Task updated! ")
}

# Hey Nik, this is the key function
# For correcting sets of neurons, groups by cell body
# fibre at once. It already works I think.
# Try running on some example neurons, read maybe using
# neuprint_read_neurons.
# hidden
correct_group <- function(neurons, brain = NULL){
  ### Cycle
  continue = TRUE
  while (isTRUE(continue)) {
    ### Extract primary neurites, as labelled
    neuron.points = do.call(plyr::rbind.fill, lapply(neurons, function(x) x$d))
    points = nat::xyzmatrix(neuron.points)
    pnts = primary_neurite_cable(neurons, OmitFailures = TRUE)
    ### Plot
    plot3d(pnts, lwd = 2, col = hemibrainr::hemibrain_bright_colors["purple"], soma = FALSE)
    plot3d(neurons, lwd = 2, col = hemibrainr::hemibrain_bright_colors["paleorange"], soma = FALSE)
    plot3d_somas(neurons)
    ### Should we continue?
    continue <- hemibrain_choice(prompt = c("Select new soma points? yes/no "))
    if(!continue){
      break
    }
    ### Select region to reroot
    make.selection = TRUE
    while(make.selection){
      angle <- readline(prompt = "Rotate data in order to best select somas (Enter to continue) ")
      message("Draw box over true soma region")
      selection <- rgl::select3d()
      selected <- selection(points)
      selected.points <- neuron.points[selected,]
      bis <- unique(selected.points$bodyid)
      bis <- bis[!is.na(bis)]
      message("Points in ", length(bis)," neurons selected")
      make.selection <- !hemibrain_choice(prompt = c("Happy with selection? yes/no "))
    }
    ### Reroot neurons in selection
    if(!is.issue(bis)){
      for(bi in bis){
        x <- neurons[[bi]]
        eps <- nat::endpoints(x)
        ep.sel <- selection(nat::xyzmatrix(x)[eps,])
        ep.sel <- eps[ep.sel][1]
        if(is.issue(ep.sel)){
          eps <- nat::branchpoints(x)
          ep.sel <- selection(nat::xyzmatrix(x)[eps,])
          ep.sel <- eps[ep.sel][1]
        }
        if(!is.issue(ep.sel)){
          soma.id <- x$d$PointNo[match(ep.sel, 1:nrow(x$d))]
          y <- nat::as.neuron(nat::as.ngraph(x$d), origin = soma.id)
          y <- hemibrain_carryover_labels(x=x,y=y)
          y <- hemibrain_carryover_tags(x=x,y=y)
          y$tags$soma.edit <- TRUE
          y$tags$soma <- soma.id
          y = hemibrain_neuron_class(y)
          neurons[[bi]] <- y
        }
      }
    }
  }
  neurons
}

# Hey nik, this is a highly unfinished function to try and
# auto-detect the soma. It might not be a good idea to try this really.
# I have already tried to auto-detect the folder on the drive, with the neurons
# that can be read using hemibrain_neurons()
# hidden
cbf_somagroup <- function(neurons){
  if(sum(neurons[,]$soma)!=length(neurons)){
    # For every neuron with an unlabelled soma
    pnts = primary_neurite_cable(neurons, OmitFailures = TRUE)

    neurons.soma = neurons[neurons[,"soma"]]
    somapos.soma <- as.data.frame(catmaid::soma(neurons.soma), stringsAsFactors = FALSE)

    neurons.nosoma = neurons[!neurons[,"soma"]]
    somapos.soma <- as.data.frame(catmaid::soma(neurons.nosoma), stringsAsFactors = FALSE)

    ## Find closest neuron with a labelled soma, within range
    ## Move soma tag to closes leaf to that closest soma
  }
}

# hidden
type_batches <- function(ids, gs, batch_size = 10){
  pcts = unique(purify(gs[match(ids, gs$bodyid),"type"]))
  pcts = pcts[order(pcts)]
  batches = split(pcts, ceiling(seq_along(pcts)/batch_size))
  batches
}

