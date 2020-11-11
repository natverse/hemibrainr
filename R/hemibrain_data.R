#' Bodyids for neuron classes
#'
#' Bodyids for neuron classes. Hard saved equivalent for \code{\link{class2ids}}.
#' However, you should use \code{\link{class2ids}} for the most up-to-date information.
#'
#' @details neuPrint bodyids for neurons that the FlyConnectome groups at the University of Cambridge Department of Zoology,
#' have determined to fall into different classes.
#'
#' @seealso \code{\link{class2ids}}
#' @examples
#' \donttest{
#' head(rn.ids)
#' }
#'
#' @name classed.ids
#' @rdname classed.ids
"rn.ids"
#' @rdname classed.ids
"orn.ids"
#' @rdname classed.ids
"hrn.ids"
#' @rdname classed.ids
"pn.ids"
#' @rdname classed.ids
"upn.ids"
#' @rdname classed.ids
"mpn.ids"
#' @rdname classed.ids
"vppn.ids"
#' @rdname classed.ids
"dan.ids"
#' @rdname classed.ids
"mbon.ids"
#' @rdname classed.ids
"alln.ids"
#' @name classed.ids
#' @rdname classed.ids
"ton.ids"
#' @name classed.ids
#' @rdname classed.ids
"lhn.ids"
#' @name classed.ids
#' @rdname classed.ids
"dn.ids"
#' @name classed.ids
#' @rdname classed.ids
"kc.ids"
#' @name classed.ids
#' @rdname classed.ids
"apl.ids"
#' @name classed.ids
#' @rdname classed.ids
"cent.ids"
#' @name classed.ids
#' @rdname classed.ids
"lc.ids"
#' @name classed.ids
#' @rdname classed.ids
"kc.ids"

#' The hemibrain mesh
#'
#' A 3D mesh of the overall neuropil structures in the hemibrain volume as a
#' \code{\link[nat]{hxsurf}} object. Created by Philipp Schlegel from the union
#' of individual ROIs in FlyEM's hemibrain project (see details). Mesh available
#' both in raw voxels and in microns.
#'
#' @details The hemibrain mesh is watertight and contains no internal structure
#'   or non-manifold vertices. It is the right thing to use if you want to see
#'   if things are inside the hemibrain neuropil. Philipp Schlegel used the
#'   following steps to construct the mesh:
#'
#'   \enumerate{
#'
#'   \item Taking all the individual ROIs
#'
#'   \item Voxelizing them
#'
#'   \item Combining the voxels into a single image stack
#'
#'   \item Filling the gaps in between (binary operations + some manual fixing
#'   in Fiji) \item Using a marching cubes to generate a mesh from the stack
#'
#'   }
#'
#'   Note that hemibrain coordinate system does not have the anterior-posterior
#'   axis aligned with the Z axis.
#'
#' @source \url{https://neuprint.janelia.org/}
#' @seealso \code{\link{hemibrain_roi_meshes}}, \code{\link{hemibrain_al.surf}}
#' @examples
#' \donttest{
#' library(nat)
#' \dontrun{
#' plot3d(hemibrain.surf, alpha=.3, col='grey')
#'
#' # can also convert to rgl mesh3d object for more control
#' hemibrain <- as.mesh3d(hemibrain.surf)
#' rgl::wire3d(hemibrain)
#' rgl::shade3d(hemibrain, alpha=.3)
#' rgl::shade3d(hemibrain, alpha=.3, col='grey')
#' }
#' }
#' @name hemibrain.surf
#' @rdname hemibrain.surf
"hemibrain.surf"
#' @rdname hemibrain.surf
"hemibrain_microns.surf"

#' Meshes for Antennal Lobe (AL) glomeruli
#'
#' 3D meshes of 51 olfactory + 7 thermosensory antennal lobe glomeruli for the
#' the hemibrain volume. They are represented as a \code{\link[nat]{hxsurf}}
#' object, in which 58 'regions' are defined, one for each glomerulus. Meshes
#' are available both in raw voxels and in microns. You can see which glomeruli are truncated, and by
#' roughly how much by examining \code{hemibrain_glomeruli_summary}.
#'
#' @details note that hemibrain coordinate system does has the
#'   anterior-posterior axis aligned with the Y axis (rather than the Z axis,
#'   which is more commonly observed).
#'
#' @source \url{https://neuprint.janelia.org/}
#' @seealso \code{\link{hemibrain_roi_meshes}}, \code{\link{hemibrain.surf}}
#' @examples
#' \donttest{
#' nvertices(hemibrain_al.surf)
#' materials(hemibrain_al.surf)
#'
#' \dontrun{
#' plot3d(hemibrain_al.surf, alpha=.3)
#'
#' # can also plot single glomeruli
#' plot3d(subset(hemibrain_al.surf,"DM1"))
#'
#' }
#' }
#' @name hemibrain_al.surf
#' @rdname hemibrain_al.surf
"hemibrain_al.surf"
#' @rdname hemibrain_al.surf
"hemibrain_al_microns.surf"
#' @rdname hemibrain_al.surf
"hemibrain_glomeruli_summary"

#' Axon-dendrite split points calculated for all hemibrain neurons
#'
#' Axon-dendrite split points calculated for the whole hemibrain dataset,
#' hemibrain:v1.0.1, using a variety of parameters. The file used to generate this can be found in this
#' package, as \code{data-raw/hemibrain_split_all_neurons.R}. The arguments given
#' to \code{\link{flow_centrality}} are in the name of the data frame.
#' Note that, because the proportion of recovered synapses is unevenly distributed across brain
#' regions, after compartment detection, sometimes the dendrite was chosen as the neurite
#' closest to a neuron's primary branch point (data frames with 'distance' in their name), see reference.
#'
#' @references Bates, A. S., P. Schlegel, R. J. V. Roberts, N. Drummond, I. F.
#'   M. Tamimi, R. Turnbull, X. Zhao, et al. 2020. “Complete Connectomic
#'   Reconstruction of Olfactory Projection Neurons in the Fly Brain.” bioRxiv.
#'   \href{https://doi.org/10.1101/2020.01.19.911453}{doi:10.1101/2020.01.19.911453}.
#'
#' @source \url{https://www.biorxiv.org/content/10.1101/2020.01.19.911453v1}
#' @seealso \code{\link{hemibrain_flow_centrality}},
#'   \code{\link{flow_centrality}},
#'   \code{\link{hemibrain_use_splitpoints}}
#' @examples
#' \donttest{
#' head(hemibrain_splitpoints_polypre_centrifugal_distance)
#' \dontrun{
#' View(hemibrain_splitpoints_polypre_centrifugal_distance)
#' }
#' }
#' @name hemibrain_precomputed_splitpoints
#' @rdname hemibrain_precomputed_splitpoints
"hemibrain_all_splitpoints"
#' @rdname hemibrain_precomputed_splitpoints
"hemibrain_splitpoints_polypre_centrifugal_distance"
#' @rdname hemibrain_precomputed_splitpoints
"hemibrain_splitpoints_pre_centrifugal_distance"
#' @rdname hemibrain_precomputed_splitpoints
"hemibrain_splitpoints_polypre_centrifugal_synapses"

#' Neuron metrics for all neurons in the hemibrain dataset
#'
#' This \code{data.frame} give metrics for neurons, such as cable lengths and pre/post synapse numbers and densities,
#' broken down by putative compartment. E.g. axons versus dendrite. This has been calculated for neurons across the hemibrain dataset,
#' hemibrain:v1.0.1. The file used to generate this can be found in this package,
#' as \code{data-raw/hemibrain_split_all_neurons.R}. Neurons were split as described by \code{\link{hemibrain_precomputed_splitpoints}}. Note
#' that the first two columns, \code{pre} and \code{post} are generated by neuPrint. They differ from similar numbers in other columns of
#' this \code{data.frame} because 'pre' is later given as the numbers of outgoing connections (not presynapses) and some erroneous pre/post synapses will
#' have been removed by \code{\link{hemibrain_skeleton_check}}. Overlap locality is calculated as the 'overlap score' (Frechter et al. 2019, eLife) between
#' a neuron's axon and dendrite.
#'
#' @references Frechter, Shahar, Alexander Shakeel Bates, Sina Tootoonian, Michael-John Dolan, James D. Manton, Arian Rokkum Jamasb, Johannes Kohl, Davi Bock, and Gregory Sxe Jefferis. 2019.
#' “Functional and Anatomical Specificity in a Higher Olfactory Centre.” eLife 8 (May). https://doi.org/10.7554/eLife.44590.
#'
#' @source \url{https://elifesciences.org/articles/44590}
#' @seealso \code{\link{hemibrain_flow_centrality}}, \code{\link{hemibrain_precomputed_splitpoints}}, \code{\link{axonic_cable}}, \code{\link{dendritic_cable}}
#' @examples
#' \donttest{
#' head(hemibrain_metrics)
#' hist(hemibrain_metrics$overlap_locality,breaks = 1000)
#' }
#' @name hemibrain_metrics
#' @rdname hemibrain_metrics
"hemibrain_metrics"
#' @rdname hemibrain_metrics
"hemibrain_metrics_polypre_centrifugal_distance"
#' @rdname hemibrain_metrics
"hemibrain_metrics_polypre_centrifugal_synapses"

#' Hemilineages names from T. Lee, K. Ito and V. Hartenstein
#'
#' Hemilineages names and putative transmitter.
#' Transmitter guesses are based on observing at least one cell type in the literature from the hemilineage in question, known to express the given transmitter
#'
#' @name hemilineages
#' @rdname hemilineages
"hemibrain_hemilineages"

#' Hemibrain neurons in different Olfactory layers
#'
#' Results of an infection model used by Philipp Schlegel.
#'
#' @name hemibrain_olfactory_layers
#' @rdname hemibrain_olfactory_layers
"hemibrain_olfactory_layers"

#' Soma / best 'root point' locations for hemibrain neurons
#'
#' Somas have been manually tagged by the FlyEM project.
#' However, in some cases the roots are wrong, or somas are outside of the volume. The Fly Connectome team at the University of Cambridge has
#' manually tagged better root points for these neurons. Soma locations, updated by this work, and saved as \code{hemibrain_somas} in this package.
#' A user can read somas from here, or get a bleeding edge version from
#' the \href{https://docs.google.com/spreadsheets/d/1YjkVjokXL4p4Q6BR-rGGGKWecXU370D1YMc1mgUYr8E/edit#gid=1524900531}{Google Sheet}. This information is used
#' by \code{\link{hemibrain_read_neurons}} to re-root hemibrain neurons read from neuprint.
#'
#' @details An interactive root point fixing pipeline \code{\link{hemibrain_adjust_saved_somas}} was used to identify correct and
#' incorrect soma positions within a group. Specifically, the DBSCAN algorithm is used. All neurons within the hemibrain data set
#' are either grouped by their CBF, or, when no CBF annotation is available, an all by all NBLAST clustering was done for neurons without a
#' CBF and this was used to group neurons into similar morphological groups. DBSCAN within one of these groups then clusters the positions of the somas,
#' returning one or more clusters of actual soma positions, and a cluster of ‘noise’ positions which do not fit any cluster.
#' The user can then choose to manually fix either the noise somas, or the entire cluster, or subset of clusters if multiple potential soma clusters are returned.
#'
#' In the instance a manual error has been made when labelling a soma position, the somas can be fixed using:
#'
#' \code{hemibrainr::adjust_saved_somas(bodyid = #####)}
#'
#' And selecting the neuron method when prompted and run through the pipeline instructions provided.
#'  When asked if you wish to try using DBSCAN on the neurons say no. The pipeline will then try
#'  to suggest a soma position based on a global clustering of those neurons which are already identified as correct.
#'
#' @return a \code{data.frame} with the 'best case' root points for each large hemibrain neuron. This is typically the soma, where one has been tagged by the FlyEM project.
#' However, for many neurons their soma has been truncated and a leaf node on the cell body fiber has been chosen. The columns to this \code{data.frame} denote:
#' \itemize{
#'
#'   \item{"bodyid"} { - the unique body ID identifier for a neuron.}
#'
#'   \item{"position"}{ - the node ID of the rootpoint. This is
#'   for the skeleton associated with \code{bodyid}, as read from neuPrint.}
#'
#'   \item{"X"}{ - x coordinate for the root point.}
#'
#'   \item{"Y"}{ - y coordinate for the root point.}
#'
#'   \item{"Z"}{ - z coordinate for the root point.}
#'
#'   \item{"soma.edit"}{ - logical. If the soma is 300 voxels away from the soma/root for the naive neuprint skeleton,
#'   we consider the soma to have been 'edited' for the better.
#'   I.e. these are the somas that were fixed by
#'   Cambridge Drosophila Connectomics using \code{\link{hemibrain_adjust_saved_somas}}.}
#'
#'   \item{"fixed"}{ - notes if there was a problem identifying the soma, eg, is the neuron only a fragment, or just the tract to the soma is visible.}
#'
#'   \item{"cellBodyFiber"}{ - The cell body fiber for a neuron, as read from neuPrint and annotated by a team under \href{https://www.janelia.org/people/kei-ito}{Kei Ito}.
#'   Where the cell body fiber is unknown, possible because it is an older neuron born in the embryo, course clustering was done to group them into loose groups.
#'   Here, \code{"unknown_number"} denotes the cluster.}
#'
#'}
#'
#' @name hemibrain_somas
#' @rdname hemibrain_somas
"hemibrain_somas"

#' Manually verified neuron matches between the hemibrain and FAFB dataset
#'
#' Manually verified neuron matches between the hemibrain and FAFB dataset, as well as some light level data from the package natverse-compatible \code{lhns}.
#' The Fly Connectome team at the University of Cambridge has
#' manually made cell-to-cell matches based on an initial NBLAST calculation, and given matches a qualitative score.
#' There can be no
#'   match (none), a poor match (poor) an okay match (okay) or an exact match (good). As a
#'   rule of thumb, a poor match could be a neuron from a very similar same cell
#'   type or a highly untraced neuron that may be the correct cell type. An okay
#'   match should be a neuron that looks to be from the same morphological cell
#'   type but there may be some discrepancies in its arbour. An exact match is a
#'   neuron that corresponds well between FAFB and the hemibrain data.
#' the \href{https://docs.google.com/spreadsheets/d/1OSlDtnR3B1LiB5cwI5x5Ql6LkZd8JOS5bBr-HTi0pOw/edit#gid=0}{Google Sheet} or using \code{\link{hemibrain_matches}}
#' if they have an authenticated email.
#'
#' @name hemibrain_matched
#' @rdname hemibrain_matched
"hemibrain_matched"
#' @rdname hemibrain_matched
"lm_matched"

#' Meta information of certain neuron classes
#'
#' For olfactory projection neurons (PNS), antennal lobe local neurons (ALLN),
#' non-MB third order olfactory neurons (TONs), lobula connection visual neurons (LC),
#' Kenyon cells (KC) and MB output neurons (MBONs).
#' Collated by the Flyconnectome group.
#'
#' @name info
#' @rdname info
"pn.info"
#' @rdname info
"mbon.info"
#' @rdname info
"ton.info"
#' @rdname info
"lc.info"
#' @rdname info
"alln.info"
#' @rdname info
"kc.info"
#' @rdname info
"cent.info"
