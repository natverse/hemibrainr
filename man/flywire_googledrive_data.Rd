% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/flywire_googledrive.R, R/flywire_synapses.R
\name{flywire_googledrive_data}
\alias{flywire_googledrive_data}
\alias{flywire_meta}
\alias{flywire_failed}
\alias{flywire_contributions}
\alias{flywire_ids}
\alias{flywire_elist}
\alias{flywire_connections}
\alias{flywire_synapses}
\title{Read precomputed flywire data from the hemibrainr Google Drive}
\usage{
flywire_meta(local = FALSE, folder = "flywire_neurons/", sql = FALSE, ...)

flywire_failed(local = FALSE, folder = "flywire_neurons/", sql = FALSE, ...)

flywire_contributions(
  local = FALSE,
  folder = "flywire_neurons/",
  sql = TRUE,
  ...
)

flywire_ids(local = FALSE, folder = "flywire_neurons/", sql = FALSE, ...)

flywire_elist(local = FALSE, folder = "flywire_neurons/", sql = TRUE, ...)

flywire_connections(
  local = FALSE,
  folder = "flywire_neurons/",
  sql = TRUE,
  ...
)

flywire_synapses(
  local = FALSE,
  folder = "flywire_neurons/",
  simplified = FALSE,
  sql = TRUE,
  ...
)
}
\arguments{
\item{local}{\code{FALSE} or path. By default (\code{FALSE}) data is read from \code{options()$remote_connectome_data}),
but the user can specify an alternative path.}

\item{folder}{A subfolder on the hemibrain team drive or your local data folder
containing the data object to read.}

\item{sql}{logical. Whether not to read the desired data from an \code{SQL} database. This can save you from having to load
a lot of information into memory, and it works well with \code{dplyr} pipes. If \code{FALSE} the relevant \code{.csv} is read from
the connected google drive (see \code{\link{hemibrainr_set_drive}}) and (see \code{\link{hemibrainr_rclone}}) and loaded into memory,
which can take some time.}

\item{...}{if \code{sql=TRUE}, methods passed to \code{dplyr::tbl}.}

\item{simplified}{flywire synapses, where nearby Buhmann predicted synapses (often multiple per real synapses) are algorithmically combined into approximate single synapses
using \code{hemibrainr:::flywire_synapse_simplify}.}
}
\value{
a \code{data.frame}. Depending on which synapse function was called, it can contain the columns:

\itemize{

  \item{"flywire_xyz"} { - coordinates of a point in the neuron in flywire voxel space. XYZ, separated by a semicolon.}

  \item{"root_id"}{ - the unique ID associated with this flywire neuron. This ID changes every time a neuron is, even slightly, modified. So it is an unstable identifier.
  This is why \code{flywire_xyz} is sometimes used.}

  \item{"fw.x"}{ - the x coordinate of a point in the flywire neuron, in flywire voxel space..}

  \item{"fw.y"}{ - the y coordinate of a point in the flywire neuron, in flywire voxel space..}

  \item{"fw.z"}{ - the z coordinate of a point in the flywire neuron, in flywire voxel space..}

  \item{"user_name"}{ - the name of the user who made the number of edits given in this row.}

  \item{"edits"}{ - the number of edits (merges, splits, etc.) made by a user for the given \code{root_id}.}

  \item{"proportion"}{ - the proportion of total edits for this neuron, that the given user made.}

  \item{"dataset"}{ - the dataset this neuron is from, i.e. flywire.}

}
}
\description{
Read precomputed data available on the hemibrain Google Team
  Drive. (see \code{\link{hemibrainr_set_drive}}) and (see
  \code{\link{hemibrainr_rclone}}). This includes body IDs for all flywire
  neurons ((\code{flywire_ids})), and user contributions towards their
  creation (\code{flywire_contributions}), as well as Flywire related NBLAST
  scores retrieved using \code{\link{hemibrain_nblast}}.
}
\examples{
\donttest{
\dontrun{

# All flywire IDs for neurons that have a split precomputed
fw.ids = flywire_ids()

# For these flywire IDs, their meya data:
fw.meta = flywire_meta()

# For flywire IDs, which users contributed what:
fw.edits = flywire_contributions()

}}
}
\seealso{
\code{\link{hemibrain_splitpoints}},
  \code{\link{hemibrain_flow_centrality}},
  \code{\link{hemibrainr_googledrive_data}},
  \code{\link{hemibrain_metrics}}
}
