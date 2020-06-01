### Core Function

#' Manually, or semi-automatically check and edit soma possitions for hemibrain neurons
#' @param bodyids list of bodyids
#' @param cbf string of a cell body fibre
#' @param brain hemibrain surface mesh
#' @param selected_file identifier for google sheet
#' @examples
#' \dontrun{
#' if (!require(RColorBrewer)) utils::install.packages("RColorBrewer")
#' if (!require(dbscan)) utils::install.packages("dbscan")
#' #' dep <- hemibrain_adjust_saved_somas(bodyids, brain = hemibrain.surf, selected_file = 1YjkVjokXL4p4Q6BR-rGGGKWecXU370D1YMc1mgUYr8E, db = NULL)
#' }
#' @export
#' @rdname hemibrain_adjust_saved_split
##### lets tear this function apart
hemibrain_adjust_saved_somas <- function(bodyids = NULL,
                                         cbf = NULL,
                                         brain = hemibrainr::hemibrain.surf,
                                         selected_file = "1YjkVjokXL4p4Q6BR-rGGGKWecXU370D1YMc1mgUYr8E",
                                         db = NULL) {
  ### Get GoogleSheet Database
  gs = gsheet_manipulation(
    FUN = googlesheets4::read_sheet,
    ss = selected_file,
    sheet = "roots",
    return = TRUE
  )
  # how would you like to go through neurons? one at a time, or based on cell body fibre
  # or Auto, if I ever add it...
  mode <-
    must_be(prompt = "Which mode do you want to use? single neurons (s) / Cell body fibre from neurons (n) / Identified cbf (c): ", answers = c("s", "n", "c"))
  # reset3d(brain=brain)
  if (mode == "s") {
    if (is.null(bodyids)) {
      message("Please provide neuron bodyids as input, exiting pipeline")
      stop()
    }
    neurons = pipeline_read_neurons(
      batch = bodyids,
      db = NULL,
      clean = FALSE,
      motivate = FALSE
    )
    neurons = correct_singles(neurons, brain = brain)
    save_soma_to_gsheet(neurons = neurons,
                        gs = gs,
                        selected_file = selected_file)
    message("Task updated! ")
  }
  if (mode == "n") {
    if (is.null(bodyids)) {
      message("Please provide neuron bodyids as input, exiting pipeline")
      stop()
    }
    message("Collecting Cell fibre bodies from provided neurons...")
    # get cell fibre body list from neurons
    cbfs = cbf_from_bodyid(bodyids)
    # how many are their, and list them
    message(c(
      "The provided bodyids are from ",
      length(cbfs),
      " different cell fibre bodies"
    ))
    if ("unknown" %in% cbfs) {
      message(
        "This includes some neurons labeled with an 'unknown' cbf, these neurons will be checked individually"
      )
    }
    # loop over each, printing which you are on
    for (c in cbfs) {
      message(c("correcting cbf ", c))
      correct_cbf(cbf = c,
                  gs = gs,
                  selected_file = selected_file)
    }
  }
  if (mode == "c") {
    if (is.null(cbf)) {
      # input a cell fibre body
      cbfs = cbf_list()
      message(
        " You need to provide a valid cell body fibre here. If you're not sure, check hemibrainr:::cbf_list() for a list"
      )
      cbf = must_be(prompt = "Please input a Cell body fibre to correct: ",
                    answers = cbf_list())
    }
    # if yes, fix...
    message(c("Fixing Cell fibre body: ", cbf))
    correct_cbf(cbf = cbf,
                gs = gs,
                selected_file = selected_file)
  }
}


#### Correct single neurons

# hidden
#' @importFrom nat neuronlist nopen3d
#' @importFrom rgl clear3d points3d spheres3d
correct_singles <- function(neurons, brain = NULL) {
  list = 1
  if (is.neuron(neurons)) {
    neurons = neuronlist(neurons)
    list = 0
  }
  nopen3d()
  correcting = TRUE
  while (correcting) {
    for (n in neurons) {
      n.points =  nat::xyzmatrix(n)
      end_points = nat::endpoints(n)
      end_points = n.points[end_points, ]
      clear3d()
      plot3d(brain, col = "grey70", alpha = 0.1)
      #
      points3d(
        n.points,
        lwd = 2,
        col = hemibrainr::hemibrain_bright_colors["yellow"],
        soma = FALSE
      )
      plot3d_somas(n)
      fix = hemibrain_choice(prompt = "Does the soma need fixing? Current possition in Green, if present yes/no ")
      if (isTRUE(fix)) {
        make.selection = TRUE
      }
      else {
        make.selection = FALSE
      }
      while (make.selection) {
        angle = readline(prompt = "position view in order to best select soma (Enter to continue) ")
        message("Click and drag a box over an end point, in red, to select new soma location")
        clear3d()
        plot3d(brain, col = "grey70", alpha = 0.1)
        points3d(
          end_points,
          lwd = 2,
          col = hemibrainr::hemibrain_bright_colors["red"],
          soma = FALSE
        )
        points3d(
          n.points,
          lwd = 2,
          col = hemibrainr::hemibrain_bright_colors["yellow"],
          soma = FALSE
        )

        selection <- rgl::select3d()

        selected = selection(end_points)
        selected.point = end_points[selected, ]
        if (length(selected.point) != 3) {
          selected.point = selected.point[1, ]
        }
        clear3d()
        plot3d(brain, col = "grey70", alpha = 0.1)
        points3d(
          n.points,
          lwd = 2,
          col = hemibrainr::hemibrain_bright_colors["yellow"],
          soma = FALSE
        )
        spheres3d(selected.point, radius = 300, col = 'blue')
        plot3d_somas(n)
        make.selection = !hemibrain_choice(prompt = c(
          "Happy with selection? (old soma in Green, new soma in Blue) yes/no "
        ))

        # reroot neuron
        eps <- nat::endpoints(n)
        #
        ep.sel <- selection(nat::xyzmatrix(n)[eps, ])
        ep.sel <- eps[ep.sel][1]
        soma.id <- n$d$PointNo[match(ep.sel, 1:nrow(n$d))]
        # create ourneuron as a graph, with the new origin point
        y <- nat::as.neuron(nat::as.ngraph(n$d), origin = 26)
        # carryover labels and tags
        y <- hemibrain_carryover_labels(x = n, y = y)
        y <- hemibrain_carryover_tags(x = n, y = y)
        y$tags$soma = soma.id
        y$tags$soma.edit = TRUE
        y = hemibrain_neuron_class(y)
        neurons[[toString(n$bodyid)]] = y
      }
    }
    clear3d()
    plot3d(neurons, soma = TRUE, WithConnectors = FALSE)
    plot3d_somas(neurons)
    correcting = !hemibrain_choice(prompt = c(
      "Final check, are you happy with the new soma possitions? yes/no "
    ))
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
                                  clean = TRUE) {
  readfail = FALSE
  if (!is.null(db)) {
    message("Reading locally saved neurons ...")
    someneuronlist = tryCatch(
      db[as.character(batch)],
      error = function(e)
        NULL
    )
    if (is.null(someneuronlist)) {
      message("Errors reading from given neuronlist, reading batch from neuPrint instead ...")
      readfail = TRUE
    }
  }
  if (is.null(db) | readfail) {
    if (motivate) {
      plot_inspirobot()
    }
    message("Reading and manipulating neurons from neuPrint ...")
    someneuronlist = hemibrain_read_neurons(
      x = as.character(batch),
      savedir = FALSE,
      remove.bad.synapses = FALSE,
      microns = FALSE,
      clean = FALSE,
      OmitFailures = TRUE
    )
  }
  if (clean) {
    someneuronlist = hemibrain_clean_skeleton(someneuronlist, rval = "neuron")
  }
  someneuronlist = add_field_seq(someneuronlist, someneuronlist[, "bodyid"], field =
                                   "bodyid")
  someneuronlist
}

# hidden
save_soma_to_gsheet = function(neurons = NULL,
                               cbf = NULL,
                               gs = gs,
                               selected_file = selected_file) {
  ### get index of bodyids in gs
  if(!is.null(cbf)){
    c = cbf
    bodyids = subset(gs, cbf == c)$bodyid
  }
  ind = which(gs$bodyid %in% bodyids)
  update = gs[ind,]
  update$soma.checked = "TRUE"
  # for each neuron in the neuron list:
  for(n in neurons){
    # if the root point id doesn't match between the update and the neuron list
    if(n$soma != update[which(update$bodyid == n$bodyid),]$position) {
      # update the values in update with the ones from the neuron list
      update[which(update$bodyid == n$bodyid),]$position = n$soma
      update[which(update$bodyid == n$bodyid),]$X = n$d[n$soma,]$X
      update[which(update$bodyid == n$bodyid),]$Y = n$d[n$soma,]$Y
      update[which(update$bodyid == n$bodyid),]$Z = n$d[n$soma,]$Z
      update[which(update$bodyid == n$bodyid),]$soma.edit = "TRUE"
    }
  }
  range = paste0("A",ind[1],":W",ind[length(ind)])
  gsheet_manipulation(
    FUN = googlesheets4::range_write,
    ss = selected_file,
    range = range,
    data = update,
    sheet = "roots",
    col_names = FALSE
  )
}

#### Cell fibre body related functions

#' @importFrom rgl clear3d spheres3d legend3d
correct_cbf = function(cbf = cbf, gs, selected_file, neurons_from_gs = TRUE) {
  if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
    stop(
      "Please install RColorBrewer using:\n",
      call. = FALSE,
      "install.packages('RColorBrewer')"
    )
  }
  if(!isTRUE(neurons_from_gs)){
    neurons = neurons_in_cbf(cbf)
  } else {
    c = cbf
    neurons = hemibrain_read_neurons(subset(gs, cbf == c)$bodyid)
  }

  db = dbscan_neurons(neurons)

  ###
  # If only returns 0 in cluster...
  if (length(unique(db$cluster)) == 1 & unique(db$cluster) == 0) {
    message("All", length(db$cluster),  " somas have been labeled as noise")
    fix = hemibrain_choice(prompt = "Do you wish to correct each manually? yes|no")
    if (isTRUE(fix)) {
      correct_singles(neurons)
    } else {
      message("You lazy person you... exiting pipeline")
      stop()
    }
  } else {
    # how many clusters and noise points are their?
    message("dbscan returns ", length(unique(db$cluster)), " cluster(s)")
    if (0 %in% db$cluster) {
      message("This includes one cluster of ",
              sum(db$cluster == 0),
              " incorrect 'noise' somas")
    }
    ###
    # if only 1 cluster exists after removing 0's, is it correct?
    if (length(unique(db$cluster)) <= 2) {
      somas = soma_locations(neurons)
      # remove the 'noise' somas
      noise = somas[which(db$cluster == 0), ]
      somas = somas[which(db$cluster == 1), ]
      # plot and check if cluster is correct
      clear3d()
      plot3d(neurons, col = "grey70")
      spheres3d(somas[, c('X', 'Y', 'Z')], radius = 300, col = 'blue')
      spheres3d(noise[, c('X', 'Y', 'Z')], radius = 300, col = 'red')
      cluster_correct = hemibrain_choice(prompt = c(
        "Has dbscan identified the correct soma cluster (in blue) yes|no "
      ))
      ###
      # if yes, mark somas as correct and update sheet

      # if present correct 'noise' somas
      if (isTRUE(cluster_correct)) {
        if (length(noise$PointNo) > 0) {
          correct_noise = hemibrain_choice(prompt = c("Do you want to correct the incorrect somas? yes|no"))
          if (isTRUE(correct_noise)) {
            neurons[which(db$cluster == 0)] = correct_singles(neurons[which(db$cluster == 0)])
          } else {
            message("fine, be lazy then")
          }
        }
      }
    } else {
      ### if multiple potential clusters were returned
      # get colour vector of distinct colours
      qual_col_pals = RColorBrewer::brewer.pal.info[brewer.pal.info$category == 'qual', ]
      col = unlist(mapply(
        RColorBrewer::brewer.pal,
        qual_col_pals$maxcolors,
        rownames(qual_col_pals)
      ))
      # plot neurons with each cluster coloured, and select which is correct
      clear3d()
      plot3d(neurons, col = "grey70")
      count = 0
      for (c in unique(db$cluster)) {
        count = count + 1
        cluster = neurons[which(db$cluster == c)]
        cluster = soma_locations(cluster)
        spheres3d(cluster[, c('X', 'Y', 'Z')], radius = 500, col = col[count])
      }
      legend3d(
        "topright",
        legend = paste('Cluster', c(as.character(1:count))),
        pch = 10,
        col = col[1:count],
        cex = 1,
        inset = c(0.02)
      )
      clusters = as.numeric(
        must_be(
          prompt = "Which cluster is just the correctly identified somas? If none are correct, press n. Otheriwse, provide the number of the correct cluster. ",
          answers = c(as.character(1:count), "n")
        )
      )
      if (clusters == "n") {
        # if no single cluster contains just somas
        message("Time to correct each soma individually, bad luck...")
        neurons = correct_singles(neurons)
      } else {
        # you have selected this cluster
        message(
          c(
            " You have chosen cluster ",
            as.character(clusters),
            " as the correct cluster of somas..."
          )
        )
        message(c(
          "sadly, now you have to fix the other ",
          as.character(sum(db$cluster != clusters)),
          " somas individually."
        ))
        #
        neurons[which(db$cluster == clusters)] = correct_singles(neurons[which(db$cluster == clusters)])
      }
    }
  }
  save_soma_to_gsheet(neurons = neurons,
                      cbf = cbf,
                      gs = gs,
                      selected_file = selected_file)
  message("Task updated! ")
}

# get a list of all annotated cell body fibres
#' @importFrom neuprintr neuprint_fetch_custom
cbf_list = function() {
  neu_cbf_list = neuprint_fetch_custom(cypher = "MATCH (n :`hemibrain_Neuron`) RETURN DISTINCT n.cellBodyFiber AS val")
  cbf_list = c()
  for (cbf in neu_cbf_list$data) {
    cbf_list = c(cbf_list, cbf[[1]][[1]])
  }
  cbf_list
}

# get neurons for a given cbf
neurons_in_cbf = function(cbf) {
  ns = neuprintr::neuprint_search(search = cbf,
                                  field = "cellBodyFiber",
                                  all_segments = TRUE)
  neurons = pipeline_read_neurons(
    batch = ns$bodyid,
    db = NULL,
    clean = FALSE,
    motivate = FALSE
  )
  neurons
}

# get cbf for a given neuron/neuronlist
cbf_from_bodyid = function(bodyids) {
  ### Get cell body fibre information
  meta = neuprint_get_meta(bodyids)
  cbfs = unique(meta$cellBodyFiber)
  if (NA %in% cbfs) {
    cbfs = c(cbfs, "unknown")
  }
  cbfs
}

### Soma related miscelany
# get soma coords for neuron list
soma_locations = function(neurons) {
  # get soma possitions, and identify missing
  somas = data.frame(matrix(NA, nrow = length(neurons),
                            ncol = length(colnames(neurons[[1]]$d))))
  colnames(somas) = colnames(neurons[[1]]$d)
  for (n in 1:length(neurons)) {
    if (sum(is.na(neurons[[n]]$d[neurons[[n]]$soma, ])) == 0) {
      if (!"bodyid" %in% names(neurons[[n]]$d[neurons[[n]]$soma, ])) {
        add = neurons[[n]]$d[neurons[[n]]$soma, ]
        add$bodyid = neurons[[n]]$bodyid
        somas[n, ] = add
      } else {
        somas[n, ] = neurons[[n]]$d[neurons[[n]]$soma, ]
      }
    }
  }
  # remove some of the columns...
  somas = somas[, !colnames(somas) %in% c("Label", "W", "Parent")]
  somas
}

# Correct missing somas in neuron list
fix_missing_soma = function(neurons) {
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
dbscan_neurons = function(neurons,
                          eps = 1500,
                          minPts = 5) {
  if (!requireNamespace("dbscan", quietly = TRUE)) {
    stop("Please install dbscan using:\n",
         call. = FALSE,
         "install.packages('dbscan')")
  }
  # check missing somas
  neurons = fix_missing_soma(neurons)

  # get soma points
  somas = soma_locations(neurons)

  # dbscan on soma points
  set.seed(123)
  db = dbscan::dbscan(somas[, c('X', 'Y', 'Z')], eps = eps, minPts = minPts)
}
