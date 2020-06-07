#' Title heibrain_adjust_saved_somas
#'
#' @param bodyids Optional, unless using neurons method. Otherwise, list of neuron bodyids
#' @param c Optional. For use when correcting Cell body fibres. If not input, function will ask you to enter one later
#' @param brain Optional. by default will use the hemibrain surface. Other neuropil surfaces can be provided here however
#' @param selected_file Identifier for gsheet you wish to read from and write to. By default, Curated_splitpoints soma sheet
#' @param db Optional. if provided, local directory to read neurons from.
#'
#' @return
#'
#' @examples
#' @export
#' @rdname hemibrain_adjust_saved_somas
hemibrain_adjust_saved_somas = function(bodyids = NULL,
                                        c = NULL,
                                        brain = NULL,
                                        selected_file = "1YjkVjokXL4p4Q6BR-rGGGKWecXU370D1YMc1mgUYr8E",
                                        db = NULL) {
  if (is.null(brain)) {
    brain = hemibrainr::hemibrain.surf
  }
  ### Get GoogleSheet Database
  gs = gsheet_manipulation(
    FUN = googlesheets4::read_sheet,
    ss = selected_file,
    sheet = "somas",
    return = TRUE
  )
  # choose method:
  message("which mode would you like to use?")
  mode = must_be(prompt = "Neurons (n) | Cell fibre body (c)", answers = c("n", "c"))
  # neuron method implementation
  if (mode == "n") {
    neuron_method(
      bodyids = bodyids,
      brain = brain,
      gs = gs,
      ss = selected_file
    )
    # or cbf method implemntation
  } else if (mode == "c") {
    cbf_method(brain = brain,
               gs = gs,
               ss = selected_file)

  }
}

# method_functions ---------------------------------------------------
neuron_method = function(bodyids = NULL,
                         brain = NULL,
                         gs = NULL,
                         ss = NULL) {
  if (is.null(bodyids)) {
    message("Please provide neuron bodyids as input, exiting pipeline")
    stop()
  }

  # create update sheet and indicies
  ind = which(gs$bodyid %in% bodyids)
  update = gs[ind,]
  # would you like to try DBSCAN on your neurons somas?
  data = list(neurons = list(), update = update, gs_somas = data.matrix(update[c('X','Y','Z')]))
  single_or_batch = must_be(
    "Would you like to correct somas individualy (s), or try DBSCAN on them (d)?",
    answers = c("s", "d")
  )
  if (single_or_batch == "d") {
    # DBSCAN neurons...
    data = correct_DBSCAN(data = data,
                          brain = brain)
    # if no, just go through each individualy
  } else {
    # read in all neurons here to data object
    data = correct_singles(data = data,
                           brain = brain)
  }
  # write update
  data = batch_somaupdate(data = data)
  write_somaupdate(update = data$update,
                   ind = ind,
                   ss = ss)
}

cbf_method = function(cbf = NULL,
                      brain = NULL,
                      gs = NULL,
                      ss = NULL,
                      neurons_from_gs = TRUE) {
  # input a cbf
  if (is.null(cbf)) {
    # input a cell fibre body
    cbfs = cbf_list(from = "google_sheet",
                    gs = gs)
    message(" You need to provide a valid cell body fibre here.")
    c = must_be(prompt = "Please input a Cell body fibre to correct: ",
                answers = cbfs)
  }
  message(c("Fixing Cell body fibre", c))
  bodyids = subset(gs, cbf == c)$bodyid
  # create update sheet and indicies
  ind = which(gs$bodyid %in% bodyids)
  update = gs[ind,]
  # dbscan neurons
  data = list(neurons = list(), update = update, gs_somas = data,matrix(update[c('X','Y','Z')]))
  data = correct_DBSCAN(data = data,
                        brain = brain)
  # write update
  data = batch_somaupdate(data = data)
  write_somaupdate(update = data$update,
                   ind = ind,
                   ss = ss)
}

# method_implementation ---------------------------------------------------

#' @importFrom nat neuronlist nopen3d
#' @importFrom rgl clear3d points3d spheres3d
correct_singles <- function(data = NULL,
                            brain = NULL,
                            subset = NULL) {
  list = 1
  if (is.neuron(data$neurons)) {
    data$neurons = neuronlist(data$neurons)
    list = 0
  }
  if (!is.null(subset)) {
    if (length(data$neurons) == 0){
      N_all = pipeline_read_neurons(data$update$bodyid[subset])
    } else {
      N_all = data$neurons[subset]
    }
  } else {
    if (length(data$neurons == 0)){
      N_all = pipeline_read_neurons(data$update$bodyid)
    } else {
      N_all = data$neurons
    }
  }
  nopen3d()
  correcting = TRUE
  while (correcting) {
    for (n in N_all) {
      n.points =  nat::xyzmatrix(n)
      end_points = nat::endpoints(n)
      end_points = n.points[end_points,]
      clear3d()
      plot3d(brain, col = "grey70", alpha = 0.1)
      #
      points3d(
        n.points,
        lwd = 2,
        col = hemibrainr::hemibrain_bright_colors["cyan"],
        soma = FALSE
      )
      plot3d_somas(n)
      fix = hemibrain_choice(prompt = "Does the soma need fixing? Current possition in Green, if present. yes|no ")
      if (isTRUE(fix)) {
        make.selection = TRUE
      } else {
        make.selection = FALSE
      }
      while (make.selection) {
        c = hemibrain_choice(prompt = "Do you think the cbf is correct? yes|no ")
        if (!isTRUE(c)) {
          message("making note of possibly incorrect cbf")
          data$update[which(data$update$bodyid == n$bodyid), ]$wrong.cbf = "TRUE"
        }
        f = hemibrain_choice(prompt = "can the soma be easily identified? yes|no ")
        if (!isTRUE(f)) {
          message("passing neuron, making note that soma can't be fixed this way")
          data$update[which(data$update$bodyid == n$bodyid), ]$unfixed = "TRUE"
          make.selection = FALSE
          next
        }
        angle = readline(prompt = "position view in order to best select soma (Enter to continue) ")
        message("Click and drag a box over an end point, in red, to select new soma location")
        clear3d()
        plot3d(brain, col = "grey70", alpha = 0.1)
        points3d(end_points,
                 lwd = 2,
                 col = hemibrain_bright_colors["red"],
                 soma = FALSE)
        points3d(n.points,
                 lwd = 2,
                 col = hemibrain_bright_colors["cyan"],
                 soma = FALSE)

        selection <- rgl::select3d()

        selected = selection(end_points)
        selected.point = end_points[selected,]
        if (length(selected.point) != 3) {
          selected.point = selected.point[1,]
        }
        clear3d()
        plot3d(brain, col = "grey70", alpha = 0.1)
        points3d(n.points,
                 lwd = 2,
                 col = hemibrain_bright_colors["cyan"],
                 soma = FALSE)
        spheres3d(selected.point, radius = 300, col = 'blue')
        plot3d_somas(n)
        make.selection = !hemibrain_choice(prompt = c(
          "Happy with selection? (old soma in Green, new soma in Blue) yes/no "
        ))
        # reroot neuron
        y = reroot_from_selection(n, selection)
        data$neurons[[toString(n$bodyid)]] = y
      }
    }
    clear3d()
    plot3d(data$neurons, soma = TRUE, WithConnectors = FALSE)
    plot3d_somas(data$neurons)
    correcting = !hemibrain_choice(prompt = c(
      "Final check, are you happy with the new soma possitions? yes/no "
    ))
  }
  if (list == 0) {
    data$neurons = hemibrain_neuron_class(data$neurons)
  }
  data
}

#' @importFrom rgl clear3d spheres3d legend3d
correct_DBSCAN = function(data = NULL,
                          brain = NULL) {
  if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
    stop(
      "Please install RColorBrewer using:\n",
      call. = FALSE,
      "install.packages('RColorBrewer')"
    )
  }
  # Run dbscan on neurons (soma possitions, without readin in gsheet)
  db = dbscan_neurons(somas = data$gs_somas,
                      brain = brain)

  # If only returns 1 cluster...
  if (length(unique(db$cluster)) == 1) {
    if (unique(db$cluster) == 0) {
      message("All",
              length(db$cluster),
              " somas have been labeled as noise")
      # read in neurons
      data$neurons = pipeline_read_neurons(update$bodyid, clean = FALSE)
      # plot and check if cluster is correct
      clear3d()
      plot3d(brain, col = "grey70", alpha = 0.1)
      plot3d(data$neurons, col = "grey70")
      spheres3d(data$gs_somas, radius = 300, col = 'red')
      fix = hemibrain_choice(prompt = "Do you wish to correct each manually? yes|no ")
      if (isTRUE(fix)) {
        correct_singles(data = data,
                        brain = brain)
      } else {
        message("You lazy person you... exiting pipeline")
        stop()
      }
    } else if (unique(db$cluster) == 1) {
      message("all neurons seem to have a labeled soma, and form a single cluster. Good times")
      # plot and check if cluster is correct
      clear3d()
      plot3d(brain, col = "grey70", alpha = 0.1)
      # plot3d(data$neurons, col = "grey70")
      spheres3d(somas[, c('X', 'Y', 'Z')], radius = 300, col = 'blue')
      cluster_correct = hemibrain_choice(prompt = c(
        "Has dbscan identified the correct soma cluster (in blue) yes|no "
      ))
      message("Would you like to plot the neurons along with the soma, to double check?")
      double_check = hemibrain_choice(prompt = "If you are confident the somas are correct, say no... (yes|no): ")
      if(isTRUE(double_check)){
        data$neurons = pipeline_read_neurons(data$update$bodyids)
        # plot and check if cluster is correct
        clear3d()
        plot3d(brain, col = "grey70", alpha = 0.1)
        plot3d(data$neurons, col = "grey70")
        spheres3d(somas[, c('X', 'Y', 'Z')], radius = 300, col = 'blue')
        cluster_correct = hemibrain_choice(prompt = c(
          "Has dbscan identified the correct soma cluster (in blue) yes|no "
        ))
      }
      if (!isTRUE(cluster_correct)) {
        fix = hemibrain_choice(
          prompt = c(
            "Do you want to manually correct each of the neurons in the cluster? yes|no "
          )
        )
        if (isTRUE(fix)) {
          data = correct_singles(data = data, brain = brain)
        } else {
          message("You lazy person you... exiting pipeline")
          stop()
        }
      }
    }
    # if multiple clusters are returned
  } else {
    # how many clusters and noise points are their?
    message("dbscan returns ", length(unique(db$cluster)), " cluster(s)")
    if (0 %in% db$cluster) {
      message("This includes one cluster of ",
              sum(db$cluster == 0),
              " incorrect 'noise' somas")
    }
    # if only a single cluster and noise
    if ((length(unique(db$cluster)) == 2) & (0 %in% db$cluster)) {
      somas = data$gs_somas
      # remove the 'noise' somas
      noise = somas[which(db$cluster == 0),]
      somas = somas[which(db$cluster == 1),]
      # plot and check if cluster is correct
      clear3d()
      plot3d(brain, col = "grey70", alpha = 0.1)
      # plot3d(data$neurons, col = "grey70")
      spheres3d(somas[, c('X', 'Y', 'Z')], radius = 300, col = 'blue')
      spheres3d(noise[, c('X', 'Y', 'Z')], radius = 300, col = 'red')
      cluster_correct = hemibrain_choice(prompt = c(
        "Has dbscan identified the correct soma cluster (in blue) yes|no "
      ))
      message("Would you like to plot the neurons along with the soma, to double check?")
      double_check = hemibrain_choice(prompt = "If you are confident the somas are correct, say no... (yes|no): ")
      if(isTRUE(double_check)){
        data$neurons = pipeline_read_neurons(data$update$bodyids)
        # plot and check if cluster is correct
        clear3d()
        plot3d(brain, col = "grey70", alpha = 0.1)
        plot3d(data$neurons, col = "grey70")
        spheres3d(somas[, c('X', 'Y', 'Z')], radius = 300, col = 'blue')
        cluster_correct = hemibrain_choice(prompt = c(
          "Has dbscan identified the correct soma cluster (in blue) yes|no "
        ))
      }
      ###
      # if present correct 'noise' somas
      if (isTRUE(cluster_correct)) {
        if (length(noise$X) > 0) {
          correct_noise = hemibrain_choice(prompt = c("Do you want to correct the incorrect somas? yes|no "))
          if (isTRUE(correct_noise)) {
            data = correct_singles(
              data = data,
              brain = brain,
              subset = which(db$cluster == 0)
            )
          } else {
            message("fine, be lazy then")
          }
        }
      }
      # if multiple potential clusters are returned
    } else {
      # get colour vector of distinct colours
      qual_col_pals = RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
      col = unlist(mapply(
        RColorBrewer::brewer.pal,
        qual_col_pals$maxcolors,
        rownames(qual_col_pals)
      ))
      # plot neurons with each cluster coloured, and select which is correct
      clear3d()
      plot3d(brain, col = "grey70", alpha = 0.1)
      # plot3d(data$neurons, col = "grey70")
      count = 0
      for (c in unique(db$cluster)) {
        count = count + 1
        cluster = list(neurons = data$neurons[which(db$cluster == c)])
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
      # Choose a cluster of correct somas
      message("Which cluster is just the correctly identified somas?")
      clusters = as.numeric(
        must_be(
          prompt =  "If none are correct, press n. Otheriwse, provide the number of the correct cluster. ",
          answers = c(as.character(1:count), "n")
        )
      )
      any_more = hemibrain_choice(prompt = "Are any other clusters also likely all somas? yes|no ")
      while (isTRUE(any_more)) {
        additional_clusters = as.numeric((
          must_be(prompt = "choose another cluster. ",
                  answers = as.character(1:count))
        ))
        clusters = c(clusters, additional_clusters)
        any_more = hemibrain_choice(prompt = "Are any other clusters also all somas? yes|no ")
      }
      message("Would you like to plot the neurons along with the soma, to double check?")
      double_check = hemibrain_choice(prompt = "If you are confident the somas are correct, say no... (yes|no): ")
      if(isTRUE(double_check)){
        data$neurons = pipeline_read_neurons(data$update$bodyids)
        # plot and check if cluster is correct
        clear3d()
        plot3d(brain, col = "grey70", alpha = 0.1)
        plot3d(data$neurons, col = "grey70")
        spheres3d(somas[, c('X', 'Y', 'Z')], radius = 300, col = 'blue')
        cluster_correct = hemibrain_choice(prompt = c(
          "Has dbscan identified the correct soma cluster(s) (in blue) yes|no "
        ))
      }
      if (clusters == "n") {
        # if no single cluster contains just somas
        message("Time to correct each soma individually, bad luck...")
        data$neurons = correct_singles(data = data, brain = brain)
      } else {
        message(c(
          "sadly, now you have to fix the other ",
          as.character(sum(db$cluster != clusters)),
          " somas individually."
        ))
        #
        data = correct_singles(
          data = data,
          brain = brain,
          subset = which(db$cluster != clusters)
        )
      }
    }
  }
  data
}


# read/write --------------------------------------------------------------
chunks = function(ind) {
  bool = diff(ind) == 1
  # get the indicies of FALSE values in bool
  jump = which(bool == FALSE)
  # for each FALSE...
  chunks = list(ind[1:jump[1]])
  for (j in 1:(length(jump) - 1)) {
    sub = list(ind[jump[j] + 1]:ind[jump[j + 1]])
    chunks = c(chunks, sub)
  }
  sub = list(ind[(jump[length(jump)] + 1):length(ind)])
  chunks = c(chunks, sub)
}

batch_somaupdate = function(data) {
  for (n in data$neurons) {
    data$update$soma.checked = "TRUE"
    # if the root point id doesn't match between the update and the neuron list
    if (n$soma != data$update[which(data$update$bodyid == n$bodyid), ]$position) {
      # update the values in update with the ones from the neuron list
      data$update[which(data$update$bodyid == n$bodyid), ]$position = n$soma
      data$update[which(data$update$bodyid == n$bodyid), ]$X = n$d[n$soma, ]$X
      data$update[which(data$update$bodyid == n$bodyid), ]$Y = n$d[n$soma, ]$Y
      data$update[which(data$update$bodyid == n$bodyid), ]$Z = n$d[n$soma, ]$Z
      data$update[which(data$update$bodyid == n$bodyid), ]$soma.edit = "TRUE"
    }
  }
  data
}

write_somaupdate = function(update = NULL,
                            ind = NULL,
                            ss = NULL) {
  # if ind is consectutive
  if (isTRUE(all(diff(ind) == 1))) {
    range = paste0("A", ind[1] + 1, ":K", ind[length(ind)] + 1)
    gsheet_manipulation(
      FUN = googlesheets4::range_write,
      ss = ss,
      range = range,
      data = update,
      sheet = "somas",
      col_names = FALSE
    )
  } else {
    chunks = chunks(ind)
    for (i in 1:length(chunks)) {
      ind = chunks[[i]]
      range = paste0("A", ind[1] + 1, ":K", ind[length(ind)] + 1)
      gsheet_manipulation(
        FUN = googlesheets4::range_write,
        ss = ss,
        range = range,
        data = update,
        sheet = "somas",
        col_names = FALSE
      )
    }
  }
}

# dbscan_utils ------------------------------------------------------------
# run dbscan on a set of neurons
dbscan_neurons = function(data = NULL,
                          somas = NULL,
                          eps = 1500,
                          minPts = 5,
                          brain = brain) {
  if (!requireNamespace("dbscan", quietly = TRUE)) {
    stop("Please install dbscan using:\n",
         call. = FALSE,
         "install.packages('dbscan')")
  }
  # check missing somas
  if(!is.null(data)){
    data = fix_missing_soma(data = data,
                            brain = brain)
  }
  # get soma points
  if (is.null(somas)){
    somas = soma_locations(data = data)
    somas = somas[, c('X', 'Y', 'Z')]
  }
  # dbscan on soma points
  set.seed(123)
  db = dbscan::dbscan(x = data.matrix(somas),
                      eps = eps,
                      minPts = minPts)
}

# get soma coords for neuron list
soma_locations = function(data = NULL) {
  # get soma possitions, and identify missing
  somas = data.frame(matrix(
    NA,
    nrow = length(data$neurons),
    ncol = length(colnames(data$neurons[[1]]$d))
  ))
  colnames(somas) = colnames(data$neurons[[1]]$d)
  for (n in 1:length(data$neurons)) {
    if (sum(is.na(data$neurons[[n]]$d[data$neurons[[n]]$soma,])) == 0) {
      if (!"bodyid" %in% names(data$neurons[[n]]$d[data$neurons[[n]]$soma,])) {
        add = data$neurons[[n]]$d[data$neurons[[n]]$soma,]
        add$bodyid = data$neurons[[n]]$bodyid
        somas[n,] = add
      } else {
        somas[n,] = data$neurons[[n]]$d[data$neurons[[n]]$soma,]
      }
    }
  }
  # remove some of the columns...
  somas = somas[,!colnames(somas) %in% c("Label", "W", "Parent")]
  somas
}

# Correct missing somas in neuron list
fix_missing_soma = function(data = NULL,
                            brain = NULL) {
  somas = soma_locations(data = data)

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
      data = correct_singles(data = data,
                             brain = brain,
                             subset = which(is.na(somas[1])))
    } else {
      message(" Fine, be lazy - we will use the StartPoint of the neurons with no soma")
      ind = which(is.na(somas[1]))
      for (i in ind)  {
        # use the neurons strt point
        data$neurons[[i]]$soma = data$neurons[[i]]$StartPoint
      }

    }
  }
  data
}
# neuron_utils ------------------------------------------------------------

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

#' @importFrom nat as.neuron endpoints xyzmatrix as.ngraph
reroot_from_selection = function(neuron = NULL,
                                 selection = NULL) {
  # reroot neuron
  eps <- nat::endpoints(neuron)
  #
  ep.sel <- selection(nat::xyzmatrix(neuron)[eps,])
  ep.sel <- eps[ep.sel][1]
  soma.id <- neuron$d$PointNo[match(ep.sel, 1:nrow(neuron$d))]
  # create ourneuron as a graph, with the new origin point
  y <- nat::as.neuron(nat::as.ngraph(neuron$d), origin = soma.id)
  # carryover labels and tags
  y <- hemibrain_carryover_labels(x = neuron, y = y)
  y <- hemibrain_carryover_tags(x = neuron, y = y)
  y$tags = as.list(y$tags)
  y$tags["soma"] = soma.id
  y$tags["soma.edit"] = TRUE
  y$soma = y$tags$soma
  y = hemibrain_neuron_class(y)
}
# cbf_utils ---------------------------------------------------------------

#' @importFrom neuprintr neuprint_fetch_custom
cbf_list = function(from = "google_sheet",
                    gs = NULL) {
  if (from == "google_sheet") {
    if (is.null(gs)) {
      message("No Google Sheet provided, falling back to neuprint meta info")
      neu_cbf_list = neuprint_fetch_custom(cypher = "MATCH (n :`hemibrain_Neuron`) RETURN DISTINCT n.cellBodyFiber AS val")
      cbf_list = c()
      for (cbf in neu_cbf_list$data) {
        cbf_list = c(cbf_list, cbf[[1]][[1]])
      }
      cbf_list
    } else {
      cbf_list = unique(gs$cbf)
    }
  } else if (from == "neuprint") {
    neu_cbf_list = neuprint_fetch_custom(cypher = "MATCH (n :`hemibrain_Neuron`) RETURN DISTINCT n.cellBodyFiber AS val")
    cbf_list = c()
    for (cbf in neu_cbf_list$data) {
      cbf_list = c(cbf_list, cbf[[1]][[1]])
    }
    cbf_list
  }
}

### fix so works from gs
# get neurons for a given cbf
neurons_in_cbf = function(c = NULL,
                          gs = NULL) {
  if (is.null(gs)) {
    ns = neuprintr::neuprint_search(search = c,
                                    field = "cellBodyFiber",
                                    all_segments = TRUE)
    neurons = pipeline_read_neurons(
      batch = ns$bodyid,
      db = NULL,
      clean = FALSE,
      motivate = FALSE
    )
  } else {
    bodyids = subset(gs, gs$cbf == c)$bodyid
    neurons = pipeline_read_neurons(
      batch = bodyids,
      db = NULL,
      clean = FALSE,
      motivate = FALSE
    )
  }
  neurons
}
