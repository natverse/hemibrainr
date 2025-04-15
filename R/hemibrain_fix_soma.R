#' Adjust and correct hemibrain neurons somas.
#'
#' @description Currently provides three "modes": correct single neurons, correct cell body fibers (cbf), and correct from \href{https://docs.google.com/spreadsheets/d/1YjkVjokXL4p4Q6BR-rGGGKWecXU370D1YMc1mgUYr8E/edit#gid=1524900531}{Google Sheet}.
#'   Correct singles will allow you to manually correct each neuron in a provided list of bodyids. If requested, DBSCAN can be used
#'       to try to predict potential soma positions against already correctly identified somas, based on the provided Google sheet.
#'   Correct cell body fibers asks for a cbf to be inputted based on those labelled in the Google sheet. DBSCAN is then used to
#'       cluster the already annotated soma positions on the \href{https://docs.google.com/spreadsheets/d/1YjkVjokXL4p4Q6BR-rGGGKWecXU370D1YMc1mgUYr8E/edit#gid=1524900531}{Google Sheet}t within the cbf. This will allow you to identify correct
#'       soma positions, and then manually fix incorrect somas. for each neuron with an incorrect soma, a potential new soma will be
#'       suggested based on the DBSCAN clustering result.
#'   Correct Google sheet is useful when screening large numbers of neurons. Once a large group of neurons has been split into
#'       morphological clusters, using NBLAST for example, this mode will request a cluster id and apply DBSCAN to somas within
#'       the cluster. You can then quickly screen to identify correct soma clusters, and label incorrect somas to be corrected separately.
#'
#' @param bodyids Optional, unless using neurons method. Otherwise, list of neuron bodyids
#' @param c Optional. For use when correcting Cell body fibers. If not input, function will ask you to enter one later
#' @param brain Optional. by default will use the hemibrain surface. Other neuropil surfaces can be provided here however
#' @param selected_file Identifier for gsheet you wish to read from and write to. By default, Curated_splitpoints soma sheet
#' @param db Optional. if provided, local directory to read neurons from.
#' @param plot_sample logical. TRUE by default, choose if you wish to plot a random subset of neurons within a DBSCAN cluster
#' @param eps the distance in nanometres used by DBSCAN to form clusters. 1500 by default
#' @param minPts The minimum number of points needed to form a cluster using DBSCAN. 5 by default
#' @param neurons_from_gsheet logical, TRUE by default. If true, will collect neurons based on bodyids in the Google sheet,
#' otherwise, will search for them based on CBF data in neuprint.
#' @param for_Imaan extra little bit for Imaan... FALSE by default
#'
#' @return Updates \href{https://docs.google.com/spreadsheets/d/1YjkVjokXL4p4Q6BR-rGGGKWecXU370D1YMc1mgUYr8E/edit#gid=1524900531}{Google Sheet} with soma information
#'
#'
#' @export
#' @rdname hemibrain_adjust_saved_somas
hemibrain_adjust_saved_somas = function(bodyids = NULL,
                                        c = NULL,
                                        brain = NULL,
                                        selected_file = "1YjkVjokXL4p4Q6BR-rGGGKWecXU370D1YMc1mgUYr8E",
                                        db = NULL,
                                        plot_sample = TRUE,
                                        eps = NULL,
                                        minPts = NULL,
                                        neurons_from_gsheet = TRUE,
                                        for_Imaan = FALSE) {
  if (is.null(eps)) {
    eps = 1500
  }
  if (is.null(minPts)) {
    minPts = 5
  }


  if (is.null(brain)) {
    brain = hemibrainr::hemibrain.surf
  }

  # choose method:
  message("which mode would you like to use?")
  mode = must_be(prompt = "Neurons (n) | Cell fibre body (c) | gsheet (g): ", answers = c("n", "c", "g"))


  gs = gsheet_manipulation(
    FUN = googlesheets4::read_sheet,
    ss = selected_file,
    sheet = "somas",
    return = TRUE
  )
  # sort out stupid data types
  gs[which(gs$cbf == "unknown"),]$clusters = as.integer(gs[which(gs$cbf == "unknown"),]$clusters)
  # sometimes bodyids are a character with a psace in front, so fix
  gs$bodyid = trimws(gs$bodyid)
  gs$fixed = trimws(gs$fixed)
  gs$position = as.integer(gs$position)
  gs$X = as.integer(gs$X)
  gs$Y = as.integer(gs$Y)
  gs$Z = as.integer(gs$Z)
  # check for NA in coords
  gs = check_coord_nans(gs = gs, selected_file = selected_file)

  data = create_SomaData(
    gs = gs,
    c = c,
    mode = mode,
    eps = eps,
    minPts = minPts,
    brain = brain,
    plot_sample = plot_sample,
    neurons_from_gsheet = neurons_from_gsheet,
    ss = selected_file,
    bodyids = bodyids,
    for_Imaan = for_Imaan
  )

  # neuron method implementation
  if (mode == "n") {
    if (for_Imaan == TRUE) {
      ini = must_be(prompt = "Is this Nik (ND), or Imaan (IT): ",
                    answers = c("ND", "IT"))

      message(
        "Good day to you Imaan! hope you're having a good day! Todays task,
        should you choose to accept it, will be to start working on a final
        double check of the neurons which have been looked at but have been
        labelled as unfixed"
      )
      message(c("There are currently ", length(which(
        data$gs$init == ini
      )), " of these guys..."))
      batch_size =
        must_be(prompt =  "How many of these would you like to have a look at? ",
                answers = c(1:length(
                  which(data$gs$fixed == FALSE & data$gs$soma.checked == TRUE)
                )))
      data$bodyids = gs[which(data$gs$init == ini),]$bodyid[1:as.integer(batch_size)]
    }




    neuron_method(data = data,
                  gs = gs,
                  ss = selected_file)
    # or cbf method implementation
  } else if (mode == "c") {
    cbf_method(data = data,
               gs = gs,
               ss = selected_file)

  } else if (mode == "g") {
    gsheet_method(data = data,
                  gs = gs,
                  ss = selected_file)
  }

}

# Data_Objects -------------------------------------------------------------

# Create empty data object used by majority of functions as input...
create_SomaData = function(gs = NULL,
                           c = NULL,
                           mode = NULL,
                           eps = NULL,
                           minPts = NULL,
                           bodyids = NULL,
                           brain = NULL,
                           plot_sample = NULL,
                           neurons_from_gsheet = NULL,
                           ss = NULL,
                           for_Imaan = NULL) {
  # create data
  data = list(
    gs = gs,
    neurons = list(),
    subset = list(),
    update = list(),
    db = NULL,
    gs_somas = list(),
    c = c,
    mode = mode,
    ind = list(),
    eps = eps,
    minPts = minPts,
    brain = brain,
    plot_sample = plot_sample,
    neurons_from_gsheet = neurons_from_gsheet,
    ss = ss,
    bodyids = bodyids,
    for_Imaan = for_Imaan
  )

  if (!is.null(data$gs)) {
    data$gs$position = as.integer(data$gs$position)
    data$gs$X = as.integer(data$gs$X)
    data$gs$Y = as.integer(data$gs$Y)
    data$gs$Z = as.integer(data$gs$Z)
  }
  data
}

generate_update = function(data = NULL,
                           gs = gs) {
  if (data$mode == "g") {
    data$c = as.integer(data$c)
    data$bodyids = subset(gs, gs$clusters == data$c)$bodyid
  } else if (data$mode == "c") {
    data$bodyids = subset(gs, gs$cbf == data$c)$bodyid
  }

  # create update sheet and indicies
  data$ind = which(gs$bodyid %in% as.character(data$bodyids))
  data$update = gs[data$ind,]
  if (typeof(data$update$X) == "character") {
    data$update$X = as.integer(data$update$X)
  }
  if (typeof(data$update$Y) == "character") {
    data$update$Y = as.integer(data$update$Y)
  }
  if (typeof(data$update$Z) == "character") {
    data$update$Z = as.integer(data$update$Z)
  }
#  data$update$soma.checked = FALSE
#  data$update$unfixed = FALSE
#  data$update$wrong.cbf = FALSE
#  data$update$soma.edit = FALSE
  # fix data types
  data$update$position = as.integer(data$update$position)
  data$update$X = as.integer(data$update$X)
  data$update$Y = as.integer(data$update$Y)
  data$update$Z = as.integer(data$update$Z)
  # create soma position bit
  data$gs_somas = data.matrix(data$update[, c('X', 'Y', 'Z')])
  data
}
# method_functions ---------------------------------------------------
neuron_method = function(data = NULL,
                         gs = NULL,
                         ss = NULL) {
  if (is.null(data$bodyids)) {
    message("Please provide neuron bodyids as input, exiting pipeline")
    stop()
  }

  data = generate_update(gs = gs,
                         data = data)

  single_or_batch = must_be(
    "Would you like to correct somas individualy (s), or try DBSCAN on them (d)?",
    answers = c("s", "d")
  )
  if (single_or_batch == "d") {
    # DBSCAN neurons...
    data = correct_DBSCAN(data = data)
    # if no, just go through each individualy
  } else {
    # read in all neurons here to data object
    data = correct_singles(data = data)
  }
  # write update
  # data = batch_somaupdate(data = data)
  data$update$soma.checked = "TRUE"
  if (isTRUE(data$for_Imaan == TRUE)){
    data$update$init = "Fixed"
  }
  write_somaupdate(data = data)
}

cbf_method = function(data = NULL,
                      gs = NULL,
                      ss = NULL) {
  # input a cbf
  if (is.null(data$c)) {
    # input a cell fibre body
    cbfs = cbf_list(from = "google_sheet",
                    gs = gs)
    message(" You need to provide a valid cell body fibre here.")
    data$c = must_be(prompt = "Please input a Cell body fibre to correct: ",
                     answers = cbfs)
  }
  message(c("Fixing Cell body fibre ", c))


  data = generate_update(gs = gs,
                         data = data)

  data = correct_DBSCAN(data = data)
  # write update
  # data = batch_somaupdate(data = data)
  data$update$soma.checked = "TRUE"
  write_somaupdate(data = data)
}

gsheet_method = function(gs = NULL,
                         ss = NULL,
                         data = NULL) {
  # select a cluster (1:40)
  message(" Which cluster of neurons without a CBF would you like to correct?")
  data$c = as.integer(must_be(prompt = "Please input a number from 1:40 to correct: ",
                              answers = c(1:44)))
  #
  data = generate_update(gs = gs,
                         data = data)
  # dbscan neurons
  data = correct_gsheet(data = data)
  # write update
  # data = batch_somaupdate(data = data)
  data$update$soma.checked = "TRUE"
  write_somaupdate(data = data)
}

# method_implementation ---------------------------------------------------

#' @importFrom nat neuronlist nopen3d
#' @importFrom rgl clear3d points3d spheres3d
correct_singles <- function(data = NULL,
                            subset_ind = NULL) {
  if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
    stop(
      "Please install RColorBrewer using:\n",
      call. = FALSE,
      "install.packages('RColorBrewer')"
    )
  }


  list = 1
  if (is.neuron(data$neurons)) {
    data$neurons = neuronlist(data$neurons)
    list = 0
  }
  if (!is.null(subset_ind)) {
    if (length(data$neurons) == 0) {
      N_all = pipeline_read_neurons(data$bodyids[subset_ind])
    } else {
      N_all = data$neurons[subset_ind]
    }
  } else {
    if (length(data$neurons) == 0) {
      N_all = pipeline_read_neurons(data$bodyids)
    } else {
      N_all = data$neurons
    }
  }
  nopen3d()
  correcting = TRUE
  while (correcting) {
    for (n in N_all) {
      n.points =  nat::xyzmatrix(n)
      end_ind = nat::endpoints(n)
      end_points = n.points[end_ind, ]
      #
      clear3d()
      plot3d(
        n,
        lwd = 1.5,
        col = hemibrainr::hemibrain_bright_colors["cyan"],
        soma = FALSE,
        WithConnectors = FALSE,
        WithNodes = FALSE
      )
      plot3d(data$brain, col = "grey70", alpha = 0.1)
      spheres3d(data$gs[data$gs$bodyid == n$bodyid,c('X','Y','Z')], radius = 300, col = 'green')
      fix = hemibrain_choice(prompt = "Does the soma need fixing? Current possition in Green, if present. yes|no ")
      if (isTRUE(fix)) {
        make.selection = TRUE
      } else {
        make.selection = FALSE
      }

      while (make.selection) {
        if (data$mode == "c") {
          c = hemibrain_choice(prompt = "Do you think the cbf is correct? yes|no ")
          if (!isTRUE(c)) {
            message("making note of possibly incorrect cbf")
            data$update[which(data$update$bodyid == n$bodyid),]$wrong.cbf = "TRUE"
          }
        }

        f = hemibrain_choice(prompt = "can the soma, or path the soma, be easily identified? yes|no ")
        if (!isTRUE(f)) {
          message("So, what is wrong with this bad boy?")
          ans = must_be(
            prompt = "Is there no soma (n), is the neuron Bilateral (b),
            is this just a fragment(f), maybe it is truncated (t), or is it just weird(w)? ",
            answers = c("n", "b", "f", "w")
          )
          message("passing neuron, and adding note...")
          if (ans == "n") {
            data$update[which(data$update$bodyid == n$bodyid),]$fixed = "No Soma"
          } else if (ans == "b") {
            data$update[which(data$update$bodyid == n$bodyid),]$fixed = "Bilateral"
          } else if (ans == "f") {
            data$update[which(data$update$bodyid == n$bodyid),]$fixed = "Fragment"
          } else if (ans == "w") {
            data$update[which(data$update$bodyid == n$bodyid),]$fixed = "Weird"
          } else if (ans == "t") {
            data$update[which(data$update$bodyid == n$bodyid),]$fixed = "Truncated"
          }
          make.selection = FALSE
          next
        } else {
          soma = hemibrain_choice(prompt = "Is only a tract to the soma visible? yes|no")
          if (isTRUE(soma)) {
            data$update[which(data$update$bodyid == n$bodyid),]$fixed = "Tract"
          }
        }

        # if data$db exists, have a guess at a soma possition for n
        # and add a message to say so
        if (is.null(data$db)) {
          message("creating a global DBSCAN clustering...")
          data$db = dbscan::dbscan(x = data$gs[which(data$gs$soma.checked == TRUE &
                                                       data$gs$fixed == TRUE), c("X", "Y", "Z")],
                                   eps = data$eps,
                                   minPts = 3)
        }
        if (!is.null(data$db)) {
          message("trying to automatically suggest a soma...")
          sugestion = suggest_soma(data = data,
                                   points = end_points,
                                   ind = end_ind)
          if (length(sugestion$ind) > 1) {
            # multiple potential somas returned
            qual_col_pals = RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual', ]
            col = unlist(mapply(
              RColorBrewer::brewer.pal,
              qual_col_pals$maxcolors,
              rownames(qual_col_pals)
            ))
            spheres3d(sugestion[, c("X", "Y", "Z")], radius = 300, col = col[1:length(sugestion)])
            legend3d(
              "topright",
              legend = paste('Sugested soma', c(as.character(
                1:length(sugestion)
              ))),
              pch = 10,
              col = col[1:length(sugestion)],
              cex = 1,
              inset = c(0.02)
            )
            message("Which sugested soma is the correct soma?")
            soma =
              must_be(prompt =  "If none are correct, press n. Otheriwse, provide the number of the correct soma. ",
                      answers = c(as.character(1:length(
                        sugestion
                      )), "n"))
            if (soma != "n") {
              sugestion = c()
            } else {
              sugestion = sugestion[which(sugestion$sugg_clust == soma), ]
            }
          }
          if (length(sugestion$ind) == 0) {
            message("No viable soma sugested,  lets do this the old fashioned way")
          } else if (length(sugestion$ind) == 1) {
            message("plotting sugested somas in blue")
            spheres3d(sugestion[, c("X", "Y", "Z")], radius = 300, col = "blue")
            sug_correct = hemibrain_choice(prompt = "Is the sugested soma possition correct? yes|no ")
            if (isTRUE(sug_correct)) {
              # reroot neuron
              # y = reroot_from_suggestion(n, sugestion)

              # N_all[[toString(n$bodyid)]] = y
              # update the values in update with the ones from the neuron list
              data$update[which(data$update$bodyid == n$bodyid),]$position = sugestion$ind
              data$update[which(data$update$bodyid == n$bodyid),]$X = sugestion$X
              data$update[which(data$update$bodyid == n$bodyid),]$Y = sugestion$Y
              data$update[which(data$update$bodyid == n$bodyid),]$Z = sugestion$Z
              data$update[which(data$update$bodyid == n$bodyid),]$soma.edit = "TRUE"

              make.selection = FALSE
            } else {
              message("bad luck,  lets do this the old fashioned way")
            }
          }
          if (isTRUE(make.selection)) {
            angle = readline(prompt = "position view in order to best select soma (Enter to continue) ")
            message("Click and drag a box over an end point, in red, to select new soma location")
            clear3d()
            plot3d(data$brain, col = "grey70", alpha = 0.1)
            points3d(
              end_points,
              lwd = 2,
              col = hemibrain_bright_colors["red"],
              soma = FALSE,
              WithConnectors = FALSE
            )
            plot3d(
              n,
              lwd = 1.5,
              col = hemibrain_bright_colors["cyan"],
              soma = FALSE,
              WithConnectors = FALSE,
              WithNodes = FALSE
            )
            # not needed ?
            testthat::try_again(100, selection <- rgl::select3d())
            selected = selection(end_points)
            while (sum(selected) == 0) {
              message("Please select one of the points in red ")
              testthat::try_again(100, selection <- rgl::select3d())
              selected = selection(end_points)
            }
            selected.point = end_points[selected, ]
            if (length(selected.point) != 3) {
              selected.point = selected.point[1, ]
            }
            clear3d()
            plot3d(data$brain, col = "grey70", alpha = 0.1)
            plot3d(
              n,
              lwd = 2,
              col = hemibrain_bright_colors["cyan"],
              soma = FALSE,
              WithConnectors = FALSE,
              WithNodes = FALSE
            )
            spheres3d(selected.point,
                      radius = 300,
                      col = 'blue')
            spheres3d(data$gs[data$gs$bodyid == n$bodyid,c('X','Y','Z')], radius = 300, col = 'green')
            make.selection = !hemibrain_choice(prompt = c(
              "Happy with selection? (old soma in Green, new soma in Blue) yes/no "
            ))
            # reroot neuron
            # y = reroot_from_selection(n, selection)
            # N_all[[toString(n$bodyid)]] = y

            # update the values in update with the ones from the neuron list
            data$update[which(data$update$bodyid == n$bodyid),]$position = which(selected == TRUE)
            data$update[which(data$update$bodyid == n$bodyid),]$X = selected.point["X"]
            data$update[which(data$update$bodyid == n$bodyid),]$Y = selected.point["Y"]
            data$update[which(data$update$bodyid == n$bodyid),]$Z = selected.point["Z"]
            data$update[which(data$update$bodyid == n$bodyid),]$soma.edit = "TRUE"

          }
        }
      }
    }
    clear3d()
    plot3d(N_all, WithConnectors = FALSE, WithNodes = FALSE)

    spheres3d(data$update[,c('X','Y','Z')], radius = 300, col = 'green')

    correcting = !hemibrain_choice(prompt = c(
      "Final check, are you happy with the new soma possitions? yes/no "
    ))
    if (list == 0) {
      data$neurons = hemibrain_neuron_class(data$neurons)
    }
}
  data
}

#' @importFrom rgl clear3d spheres3d legend3d
correct_DBSCAN = function(data = NULL) {
  if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
    stop(
      "Please install RColorBrewer using:\n",
      call. = FALSE,
      "install.packages('RColorBrewer')"
    )
  }
  # Run dbscan on neurons (soma possitions, without reading in gsheet)
  data$db = dbscan_neurons(data = data)

  # If only returns 1 cluster...
  if (length(unique(data$db$cluster)) == 1) {
    if (unique(data$db$cluster) == 0) {
      message("All",
              length(data$db$cluster),
              " somas have been labelled as noise")
      # read in neurons
      data$neurons = pipeline_read_neurons(data$update$bodyid)
      # plot and check if cluster is correct
      clear3d()
      plot3d(data$brain, col = "grey70", alpha = 0.1)
      plot3d(data$neurons, col = "grey70")
      spheres3d(data$gs_somas, radius = 300, col = 'red')
      fix = hemibrain_choice(prompt = "Do you wish to correct each manually? yes|no ")
      if (isTRUE(fix)) {
        correct_singles(data = data)
      } else {
        message("You lazy person you... exiting pipeline")
        stop()
      }
    } else if (unique(data$db$cluster) == 1) {
      message("all neurons seem to have a labelled soma, and form a single cluster. Good times")
      # plot and check if cluster is correct
      clear3d()


      # plot subset of neurons
      if (isTRUE(data$plot_sample)) {
        if (length(data$gs_somas[which(data$db$cluster == 1)]) < 10) {
          data$subset = pipeline_read_neurons(sample(data$update[which(data$db$cluster == 1),]$bodyid, 2))
          plot3d(data$subset,
                 col = "grey70",
                 WithConnectors = FALSE)
        } else {
          data$subset = pipeline_read_neurons(sample(data$update[which(data$db$cluster == 1),]$bodyid, 5))
          plot3d(data$subset,
                 col = "grey70",
                 WithConnectors = FALSE)
        }
      }
      plot3d(data$brain, col = "grey70", alpha = 0.1)
      # plot3d(data$neurons, col = "grey70")
      spheres3d(data$gs_somas, radius = 300, col = 'blue')
      cluster_correct = hemibrain_choice(prompt = c(
        "Has dbscan identified the correct soma cluster (in blue) yes|no "
      ))
      message("Would you like to plot the neurons along with the soma, to double check?")
      double_check = hemibrain_choice(prompt = "If you are confident the somas are correct, say no... (yes|no): ")
      if (isTRUE(double_check)) {
        data$neurons = pipeline_read_neurons(data$update$bodyid)
        # plot and check if cluster is correct
        clear3d()
        plot3d(data$brain, col = "grey70", alpha = 0.1)
        plot3d(data$neurons, col = "grey70")
        spheres3d(data$gs_somas, radius = 300, col = 'blue')
        cluster_correct = hemibrain_choice(
          prompt = c(
            "Double check: Has dbscan identified the correct soma cluster (in blue) yes|no "
          )
        )
      }
      if (!isTRUE(cluster_correct)) {
        fix = hemibrain_choice(
          prompt = c(
            "Do you want to manually correct each of the neurons in the cluster? yes|no "
          )
        )
        if (isTRUE(fix)) {
          data = correct_singles(data = data)
        } else {
          message("You lazy person you... exiting pipeline")
          stop()
        }
      }
    }
    # if multiple clusters are returned
  } else {
    # how many clusters and noise points are their?
    message("dbscan returns ", length(unique(data$db$cluster)), " cluster(s)")
    if (0 %in% data$db$cluster) {
      message(
        "This includes one cluster of ",
        sum(data$db$cluster == 0),
        " incorrect 'noise' somas"
      )
    }
    # if only a single cluster and noise
    if ((length(unique(data$db$cluster)) == 2) &
        (0 %in% data$db$cluster)) {
      somas = data$gs_somas

      #
      noise = somas[which(data$db$cluster == 0), ]
      if (ncol(noise) == 1) {
        noise = t(noise)
      }
      somas = as.data.frame(somas[which(data$db$cluster == 1), ])
      # plot and check if cluster is correct
      clear3d()

      # plot subset of neurons
      if (isTRUE(data$plot_sample)) {
        if (length(data$gs_somas[which(data$db$cluster == 1)]) < 10) {
          data$subset = pipeline_read_neurons(sample(data$update[which(data$db$cluster == 1),]$bodyid, 2))
          plot3d(data$subset,
                 col = "grey70",
                 WithConnectors = FALSE)
        } else {
          data$subset = pipeline_read_neurons(sample(data$update[which(data$db$cluster == 1),]$bodyid, 5))
          plot3d(data$subset,
                 col = "grey70",
                 WithConnectors = FALSE)
        }
      }
      # plot3d(data$neurons, col = "grey70")
      plot3d(data$brain, col = "grey70", alpha = 0.1)
      spheres3d(somas, radius = 300, col = 'blue')
      spheres3d(noise, radius = 300, col = 'red')
      cluster_correct = hemibrain_choice(prompt = c(
        "Has dbscan identified the correct soma cluster (in blue) yes|no "
      ))
      message("Would you like to plot the neurons along with the soma, to double check?")
      double_check = hemibrain_choice(prompt = "If you are confident the somas are correct, say no... (yes|no): ")
      if (isTRUE(double_check)) {
        data$neurons = pipeline_read_neurons(data$update$bodyid)
        # plot and check if cluster is correct
        clear3d()
        plot3d(data$brain, col = "grey70", alpha = 0.1)
        plot3d(data$neurons, col = "grey70")
        spheres3d(data$gs_somas, radius = 300, col = 'blue')
        cluster_correct = hemibrain_choice(
          prompt = c(
            "Double check: Has dbscan identified the correct soma cluster (in blue) yes|no "
          )
        )
      }
      ###
      # if present correct 'noise' somas
      if (isTRUE(cluster_correct)) {
        if (length(nrow(noise)) > 0) {
          correct_noise = hemibrain_choice(prompt = c("Do you want to correct the incorrect somas? yes|no "))
          if (isTRUE(correct_noise)) {
            data = correct_singles(data = data,
                                   subset_ind = which(data$db$cluster == 0))
          } else {
            message("fine, be lazy then")
          }
        }
      }
      # if multiple potential clusters are returned
    } else {
      # get colour vector of distinct colours
      qual_col_pals = RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual', ]
      col = unlist(mapply(
        RColorBrewer::brewer.pal,
        qual_col_pals$maxcolors,
        rownames(qual_col_pals)
      ))
      # plot neurons with each cluster coloured, and select which is correct
      clear3d()
      plot3d(data$brain, col = "grey70", alpha = 0.1)

      # plot random subset of neurons within each ID'd cluster (not noise)
      if (isTRUE(data$plot_sample)) {
        for (c in unique(data$db$cluster)) {
          if (c != 0) {
            # plot subset of neurons
            if (length(data$gs_somas[which(data$db$cluster == c)]) < 10) {
              clust = pipeline_read_neurons(sample(data$update[which(data$db$cluster == c),]$bodyid, 2))
              plot3d(clust,
                     col = "grey70",
                     WithConnectors = FALSE)
            } else {
              clust = pipeline_read_neurons(sample(data$update[which(data$db$cluster == c),]$bodyid, 5))
              plot3d(clust,
                     col = "grey70",
                     WithConnectors = FALSE)
            }
          }
        }
      }
      count = 0
      for (c in unique(data$db$cluster)) {
        if (c != 0) {
          count = count + 1
          cluster = data$gs_somas[which(data$db$cluster == c), ]
          spheres3d(cluster, radius = 500, col = col[count])
        }
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
      clusters =
        must_be(prompt =  "If none are correct, press n. Otheriwse, provide the number of the correct cluster. ",
                answers = c(as.character(1:count), "n"))
      if (clusters != "n") {
        clusters = as.numeric(clusters)
      }
      any_more = hemibrain_choice(prompt = "Are any other clusters also likely all somas? yes|no ")
      while (isTRUE(any_more)) {
        additional_clusters = as.numeric((
          must_be(prompt = "choose another cluster. ",
                  answers = as.character(1:count))
        ))
        clusters = c(clusters, additional_clusters)
        any_more = hemibrain_choice(prompt = "Are any other clusters also all somas? yes|no ")
      }
      message("Would you like to plot ALL neurons along with the soma, to double check?")
      double_check = hemibrain_choice(prompt = "If you are confident the somas are correct, say no... (yes|no): ")
      if (isTRUE(double_check)) {
        data$neurons = pipeline_read_neurons(data$update$bodyid)
        # plot and check if cluster is correct
        clear3d()
        plot3d(data$brain, col = "grey70", alpha = 0.1)
        plot3d(data$neurons, col = "grey70")
        # colour this by the somas!
        count = 0
        for (c in unique(data$db$cluster)) {
          count = count + 1
          cluster = data$gs_somas[which(data$db$cluster == c), ]
          spheres3d(cluster, radius = 500, col = col[count])
        }
        cluster_correct = hemibrain_choice(prompt = c(
          "Has dbscan identified the correct soma cluster(s) yes|no "
        ))
      }
      if (clusters == "n") {
        # if no single cluster contains just somas
        message("Time to correct each soma individually, bad luck...")
        data = correct_singles(data = data)
      } else {
        message(c(
          "sadly, now you have to fix the other ",
          as.character(
            sum(data$db$cluster %in% clusters),
            " somas individually."
          )
        ))
        #
        data = correct_singles(data = data,
                               subset_ind = which(data$db$cluster %in% (clusters)))
      }
    }
  }
  data
}

#' @importFrom rgl clear3d spheres3d legend3d
correct_gsheet = function(data = NULL) {
  if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
    stop(
      "Please install RColorBrewer using:\n",
      call. = FALSE,
      "install.packages('RColorBrewer')"
    )
  }
  # Run dbscan on neurons (soma possitions, without reading in gsheet)
  data$db = dbscan_neurons(data = data)

  # If only returns 1 cluster...
  if (length(unique(data$db$cluster)) == 1) {
    if (unique(data$db$cluster) == 0) {
      message(
        "All",
        length(data$db$cluster),
        " somas have been labelled as noise, making a note, and move on. It's over for us"
      )
      ###
      data$update$fixed = "FALSE"
      data$update$soma.checked = "TRUE"
    } else if (unique(data$db$cluster) == 1) {
      message("all neurons seem to have a labelled soma, and form a single cluster. Good times")
      # plot and check if cluster is correct
      clear3d()
      # plot subset of neurons
      if (isTRUE(data$plot_sample)) {
        if (length(data$gs_somas[which(data$db$cluster == 1)]) < 10) {
          data$subset = pipeline_read_neurons(sample(data$update[which(data$db$cluster == 1),]$bodyid, 2))
          plot3d(data$subset,
                 col = "grey70",
                 WithConnectors = FALSE)
        } else {
          data$subset = pipeline_read_neurons(sample(data$update[which(data$db$cluster == 1),]$bodyid, 5))
          plot3d(data$subset,
                 col = "grey70",
                 WithConnectors = FALSE)
        }
      }
      plot3d(data$brain, col = "grey70", alpha = 0.1)
      # plot3d(data$neurons, col = "grey70")
      spheres3d(data$gs_somas, radius = 300, col = 'blue')
      cluster_correct = hemibrain_choice(prompt = c(
        "Has dbscan identified the correct soma cluster (in blue) yes|no "
      ))
      message("Would you like to plot the neurons along with the soma, to double check?")
      double_check = hemibrain_choice(prompt = "If you are confident the somas are correct, say no... (yes|no): ")
      if (isTRUE(double_check)) {
        data$neurons = pipeline_read_neurons(data$update$bodyid)
        # plot and check if cluster is correct
        clear3d()
        plot3d(data$brain, col = "grey70", alpha = 0.1)
        plot3d(data$neurons, col = "grey70")
        spheres3d(data$gs_somas, radius = 300, col = 'blue')
        cluster_correct = hemibrain_choice(
          prompt = c(
            "Double check: Has dbscan identified the correct soma cluster (in blue) yes|no "
          )
        )
      }
      if (!isTRUE(cluster_correct)) {
        ### note that the cluster is "wrong" - label all as noise
        data$update$fixed = "FALSE"
        data$update$soma.checked = "TRUE"
      } else {
        data$update$soma.checked = "TRUE"
      }
    }
    # if multiple clusters are returned
  } else {
    # how many clusters and noise points are their?
    message("dbscan returns ", length(unique(data$db$cluster)), " cluster(s)")
    if (0 %in% data$db$cluster) {
      message(
        "This includes one cluster of ",
        sum(data$db$cluster == 0),
        " incorrect 'noise' somas"
      )
    }
    # if only a single cluster and noise
    if ((length(unique(data$db$cluster)) == 2) &
        (0 %in% data$db$cluster)) {
      somas = data$gs_somas
      #
      noise = as.data.frame(somas[which(data$db$cluster == 0), ])
      if (length(ncol(noise)) == 1) {
        noise = t(noise)
      }
      somas = as.data.frame(somas[which(data$db$cluster == 1), ])
      # plot and check if cluster is correct
      clear3d()

      # plot subset of neurons
      if (isTRUE(data$plot_sample)) {
        if (length(data$gs_somas[which(data$db$cluster == 1)]) < 10) {
          data$subset = pipeline_read_neurons(sample(data$update[which(data$db$cluster == 1),]$bodyid, 2))
          plot3d(data$subset,
                 col = "grey70",
                 WithConnectors = FALSE)
        } else {
          data$subset = pipeline_read_neurons(sample(data$update[which(data$db$cluster == 1),]$bodyid, 5))
          plot3d(data$subset,
                 col = "grey70",
                 WithConnectors = FALSE)
        }
      }
      # plot3d(data$neurons, col = "grey70")
      plot3d(data$brain, col = "grey70", alpha = 0.1)
      spheres3d(somas, radius = 300, col = 'blue')
      spheres3d(noise, radius = 300, col = 'red')
      cluster_correct = hemibrain_choice(prompt = c(
        "Has dbscan identified the correct soma cluster (in blue) yes|no "
      ))
      message("Would you like to plot the neurons along with the soma, to double check?")
      double_check = hemibrain_choice(prompt = "If you are confident the somas are correct, say no... (yes|no): ")
      if (isTRUE(double_check)) {
        data$neurons = pipeline_read_neurons(data$update$bodyid)
        # plot and check if cluster is correct
        clear3d()
        plot3d(data$brain, col = "grey70", alpha = 0.1)
        plot3d(data$neurons, col = "grey70")
        spheres3d(somas, radius = 300, col = 'blue')
        spheres3d(noise, radius = 300, col = 'red')
        cluster_correct = hemibrain_choice(
          prompt = c(
            "Double check: Has dbscan identified the correct soma cluster (in blue) yes|no "
          )
        )
      }
      ###
      # if present correct 'noise' somas
      if (isTRUE(cluster_correct)) {
        if (length(nrow(noise)) > 0) {
          ### mark relevant somas as noise
          data$update$fixed[data$db$cluster == 0] = "FALSE"
        }
        data$update$soma.checked = "TRUE"
      }
      # if multiple potential clusters are returned
    } else {
      # get colour vector of distinct colours
      qual_col_pals = RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual', ]
      col = unlist(mapply(
        RColorBrewer::brewer.pal,
        qual_col_pals$maxcolors,
        rownames(qual_col_pals)
      ))
      # plot neurons with each cluster coloured, and select which is correct
      clear3d()
      plot3d(data$brain, col = "grey70", alpha = 0.1)

      # plot random subset of neurons within each ID'd cluster (not noise)
      if (isTRUE(data$plot_sample)) {
        for (c in unique(data$db$cluster)) {
          if (c != 0) {
            # plot subset of neurons
            if (length(data$gs_somas[which(data$db$cluster == c)]) < 10) {
              clust = pipeline_read_neurons(sample(data$update[which(data$db$cluster == c),]$bodyid, 2))
              plot3d(clust,
                     col = "grey70",
                     WithConnectors = FALSE)
            } else {
              clust = pipeline_read_neurons(sample(data$update[which(data$db$cluster == c),]$bodyid, 5))
              plot3d(clust,
                     col = "grey70",
                     WithConnectors = FALSE)
            }
          }
        }
      }
      count = 0
      for (c in unique(data$db$cluster)) {
        if (c != 0) {
          count = count + 1
          cluster = data$gs_somas[which(data$db$cluster == c), ]
          spheres3d(cluster, radius = 500, col = col[count])
        }
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
      clusters =
        must_be(prompt =  "If none are correct, press n. Otheriwse, provide the number of the correct cluster. ",
                answers = c(as.character(1:count), "n"))
      if (clusters != "n") {
        clusters = as.numeric(clusters)
      }
      any_more = hemibrain_choice(prompt = "Are any other clusters also likely all somas? yes|no ")
      while (isTRUE(any_more)) {
        additional_clusters = as.numeric((
          must_be(prompt = "choose another cluster. ",
                  answers = as.character(1:count))
        ))
        clusters = c(clusters, additional_clusters)
        any_more = hemibrain_choice(prompt = "Are any other clusters also all somas? yes|no ")
      }
      message("Would you like to plot ALL neurons along with the soma, to double check?")
      double_check = hemibrain_choice(prompt = "If you are confident the somas are correct, say no... (yes|no): ")
      if (isTRUE(double_check)) {
        data$neurons = pipeline_read_neurons(data$update$bodyid)
        # plot and check if cluster is correct
        clear3d()
        plot3d(data$brain, col = "grey70", alpha = 0.1)
        plot3d(data$neurons, col = "grey70")
        # colour this by the somas!
        count = 0
        for (c in unique(data$db$cluster)) {
          count = count + 1
          cluster = data$gs_somas[which(data$db$cluster == c), ]
          spheres3d(cluster, radius = 500, col = col[count])
        }
        cluster_correct = hemibrain_choice(prompt = c(
          "Has dbscan identified the correct soma cluster(s) yes|no "
        ))
      }
      if ("n" %in% clusters) {
        # if no single cluster contains just somas
        data$update$soma.checked = "TRUE"
        data$update$fixed = "FALSE"
        # note all as noise
      } else {
        # note the remainder as noise
        message(c("Noting  ",
                  as.character(sum(
                    !(data$db$cluster %in% clusters)
                  ),
                  " as unfixed.")))
        data$update$fixed[!(data$db$cluster %in% clusters)] = "FALSE"
        data$update$soma.checked = "TRUE"
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
  for (j in seq_along(jump[-1])) {
    sub = list(ind[jump[j] + 1]:ind[jump[j + 1]])
    chunks = c(chunks, sub)
  }
  sub = list(ind[(jump[length(jump)] + 1):length(ind)])
  chunks = c(chunks, sub)
}

check_coord_nans = function(gs = NULL, selected_file = NULL) {
  ind = which(is.na(gs[, c("X", "Y", "Z")]))
  if (length(ind) != 0) {
    for (i in ind) {
      curr = gs[i, ]
      # read in neuron
      n = pipeline_read_neurons(batch = curr$bodyid)[[1]]
      # update curr
      curr[, c("X", "Y", "Z")] = n$d[n$soma, c("X", "Y", "Z")]
      # sort data types
      curr$bodyid = as.integer(curr$bodyid)
      curr$position = as.integer(curr$position)
      curr$soma.edit = as.logical(curr$soma.edit)
      curr$soma.checked = as.logical(curr$soma.checked)
      curr$wrong.cbf = as.logical(curr$wrong.cbf)
      curr$fixed = as.logical(curr$fixed)
      # update gs
      gs[i, ] = curr
      # write to gs
      range = paste0("A", i + 1, ":M", i + 1)
      gsheet_manipulation(
        FUN = googlesheets4::range_write,
        ss = selected_file,
        range = range,
        data = curr,
        sheet = "somas",
        col_names = FALSE
      )
    }
  }
  gs
}
batch_somaupdate = function(data = NULL) {
  for (n in data$neurons) {
    data$update$soma.checked = "TRUE"
    # if the root point id doesn't match between the update and the neuron list
    if (n$soma != data$update[which(data$update$bodyid == n$bodyid),]$position) {
      # update the values in update with the ones from the neuron list
      data$update[which(data$update$bodyid == n$bodyid),]$position = n$soma
      data$update[which(data$update$bodyid == n$bodyid),]$X = n$d[n$soma,]$X
      data$update[which(data$update$bodyid == n$bodyid),]$Y = n$d[n$soma,]$Y
      data$update[which(data$update$bodyid == n$bodyid),]$Z = n$d[n$soma,]$Z
      data$update[which(data$update$bodyid == n$bodyid),]$soma.edit = "TRUE"
    }
  }
  data
}

write_somaupdate = function(data = data) {
  sheet = "somas"
  last = ":M"
  # if ind is consecutive
  if (isTRUE(all(diff(data$ind) == 1))) {
    range = paste0("A", data$ind[1] + 1, last, data$ind[length(data$ind)] + 1)
    gsheet_manipulation(
      FUN = googlesheets4::range_write,
      ss = data$ss,
      range = range,
      data = data$update,
      sheet = sheet,
      col_names = FALSE
    )
  } else {
    chunks = chunks(data$ind)
    for (i in 1:length(chunks)) {
      ind = chunks[[i]]
      range = paste0("A", ind[1] + 1, last, ind[length(ind)] + 1)
      gsheet_manipulation(
        FUN = googlesheets4::range_write,
        ss = data$ss,
        range = range,
        data = data$update[which(rownames(data$update) %in% ind), ],
        sheet = sheet,
        col_names = FALSE
      )
    }
  }
}

# dbscan_utils ------------------------------------------------------------
# run dbscan on a set of neurons
dbscan_neurons = function(data = NULL) {
  if (!requireNamespace("dbscan", quietly = TRUE)) {
    stop("Please install dbscan using:\n",
         call. = FALSE,
         "install.packages('dbscan')")
  }
  # check missing somas - remove for now?
  #  if (is.null(data)) {
  #    data = fix_missing_soma(data = data,
  #                            brain = brain)
  #  }
  # get soma points
  #  if (is.null(somas)) {
  #    somas = soma_locations(data = data)
  #    somas = data.matrix(somas[, c('X', 'Y', 'Z')])
  #  }

  # dbscan on soma points
  set.seed(123)
  db = dbscan::dbscan(x = data$gs_somas,
                      eps = data$eps,
                      minPts = data$minPts)
}

suggest_soma = function(data = NULL,
                        points = NULL,
                        ind = NULL) {
  if (!requireNamespace("dbscan", quietly = TRUE)) {
    stop("Please install dbscan using:\n",
         call. = FALSE,
         "install.packages('dbscan')")
  }
  # suggest which are part of a cluster
  ## dbscan::predict.dbscan_fast
  suggest = stats::predict(object = data$db,
                                         data = data$gs_somas,
                                         newdata = points)
  if ((length(unique(suggest)) == 1) && (unique(suggest) == 0)) {
    message("no pottential soma found")
    sugg_points = c()
    return()
  }
  sugg_points = data.frame(
    ind = ind,
    X = points[, 1],
    Y = points[, 2],
    Z = points[, 3],
    sugg_clust = suggest,
    stringsAsFactors = FALSE
  )
  sugg_points = sugg_points[which(sugg_points$sugg_clust != 0), ]
  # cut down to the just the points furthest from boundary in each cluster
  sugg_points$dist = nat::pointsinside(sugg_points[, c("X", "Y", "Z")], data$brain, rval = "distance")
  # get the minimum distance value for each cluster
  to_keep = c()
  for (i in unique(sugg_points$sugg_clust)) {
    to_keep = c(to_keep, which(sugg_points$dist == min(sugg_points$dist)))
  }
  sugg_points = sugg_points[to_keep, ]
}

# get soma coords for neuron list
soma_locations = function(data = NULL) {
  # get soma possitions, and identify missing
  somas = data.frame(matrix(
    NA,
    nrow = length(data$neurons),
    ncol = length(colnames(data$neurons[[1]]$d))
  ), stringsAsFactors = FALSE)
  colnames(somas) = colnames(data$neurons[[1]]$d)
  for (n in 1:length(data$neurons)) {
    if (sum(is.na(data$neurons[[n]]$d[data$neurons[[n]]$soma, ])) == 0) {
      if (!"bodyid" %in% names(data$neurons[[n]]$d[data$neurons[[n]]$soma, ])) {
        add = data$neurons[[n]]$d[data$neurons[[n]]$soma, ]
        add$bodyid = data$neurons[[n]]$bodyid
        somas[n, ] = add
      } else {
        somas[n, ] = data$neurons[[n]]$d[data$neurons[[n]]$soma, ]
      }
    }
  }
  # remove some of the columns...
  somas = somas[,!colnames(somas) %in% c("Label", "W", "Parent")]
  somas
}

# Correct missing somas in neuron list
fix_missing_soma = function(data = NULL) {
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
                             subset_ind = which(is.na(somas[1])))
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
                                  motivate = FALSE,
                                  clean = FALSE) {
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
  ep.sel <- selection(nat::xyzmatrix(neuron)[eps, ])
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

reroot_from_suggestion = function (neuron = NULL,
                                   suggestion = NULL) {
  y <-
    nat::as.neuron(nat::as.ngraph(neuron$d), origin = suggestion$ind)
  # carryover labels and tags
  y <- hemibrain_carryover_labels(x = neuron, y = y)
  y <- hemibrain_carryover_tags(x = neuron, y = y)
  y$tags = as.list(y$tags)
  y$tags["soma"] = suggestion$ind
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
