% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hemibrain_googledrive.R
\name{hemibrain_nblast}
\alias{hemibrain_nblast}
\title{Read precomputed NBLASTs from the hemibrainr Google Drive}
\usage{
hemibrain_nblast(
  nblast = c("hemibrain", "flywire", "flywire-sameside", "hemibrain-flywire",
    "hemibrain-fafb14", "hemibrain-flycircuit", "flywire-mirror", "flywire-spine",
    "flywire-pnt", "hemibrain-primary.neurites", "hemibrain-primary.dendrites",
    "hemibrain-axons", "hemibrain-dendrites", "hemibrain-spines", "hemibrain-tracts",
    "hemibrain-arbour", "hemibrain-simplified"),
  local = FALSE,
  fafbsegpy = TRUE
)
}
\arguments{
\item{nblast}{the NBLAST matrix you would like to retrieve, e.g.
\code{"arbours"} gives you a normalised all by all NBLAST matrix of all
branching arbour.
\code{\link{dotprops}} object.}

\item{local}{\code{FALSE} or path. By default (\code{FALSE}) data is read from \code{options()$remote_connectome_data}),
but the user can specify an alternative path.}

\item{fafbsegpy}{logical. If \code{TRUE} the NBLASTs are fetched from Philipp Schlegel's nightly  \href{https://github.com/flyconnectome/flywire_nblast_pipeline/tree/main/data}{flywire NBLAST pipeline}.}
}
\value{
a \code{matrix} with named rows and columns describing an NBLAST result.
}
\description{
Read precomputed data available on the hemibrain Google Team
  Drive. This includes body IDs for all hemibrain neurons ((\code{hemibrain_neuron_bodyids})),
  all synapses (\code{hemibrain_synapses}),
  neuron-neuron connections (\code{hemibrain_connections}) and an
  edgelist  (\code{hemibrain_elist}) for all hemibrain neurons, broken down by axon and dendrite
  assignments. NBLAST matrices for all neurons against all neurons in the
  data set are also available, including ones broken down by neuron
  compartment.
}
\details{
NBLASTs were made in \code{JRCFIB2018F} space. Hemibrain neurons first had their 'twigs' pruned:

\code{fib.twigs5 = nlapply(all.neurons.flow, prune_twigs, twig_length=5000, .parallel = TRUE, OmitFailures = TRUE)}.

Where \code{all.neurons.flow} can be called with \code{hemibrain_neurons}. Neurons were moved from \code{JRC2018Fraw}
to \code{JRC2018F} by:

\code{all.neurons.flow.microns = hemibrainr:::scale_neurons.neuronlist(fib.twigs5, .parallel = TRUE, OmitFailures = TRUE)}

For 'compartment' NBLASTs neurons' axons, dendrites, primary neurite
tracts, etc., were extracted. E.g.

\code{all.neurons.pnt = primary_neurite_cable(x = all.neurons.flow.microns, .parallel = TRUE, OmitFailures = TRUE)}
\code{all.neurons.arbour = arbour_cable(x = all.neurons.flow.microns, .parallel = TRUE, OmitFailures = TRUE)}
\code{all.neurons.tract = tract_cable(x = all.neurons.flow.microns, .parallel = TRUE, OmitFailures = TRUE)}

For the 'simplified' hemibrain neuron NBLAST, the function \code{nat::simplify_neuron} was used as so:

\code{all.neurons.simp = nat::nlapply(X = all.neurons.flow.microns, FUN = nat::simplify_neuron, n = 1, invert = FALSE, .parallel = TRUE, OmitFailures = TRUE)}

Neurons, simplified neurons and compartments were turned into vector cloud using \code{nat::dotprops} as so:

\code{hemibrain.microns.dps =dotprops(all.neurons.flow.microns, k=5, resample=1, .parallel=T, OmitFailures = T)}

NBLASTs were run as so:

\code{hemibrain.aba.mean=nat.nblast::nblast_allbyall(hemibrain.twigs5.dps,
.progress='text',.parallel=TRUE, normalisation='mean')}

NBLASTs were also run against flywire and flycircuit neurons.
Flywire neurons were pulled from the \href{https://ngl.flywire.ai/?local_id=1191adb1d6f85f5f8be0fedc675460f6}{flywire project} as meshes and skeletonised
using \code{fafbseg::skeletor}.

The code for this pipeline can be found as a \href{https://github.com/flyconnectome/fafbpipeline}{flyconnectome GitHub repository}
}
\examples{
\donttest{
\dontrun{

# A normalised neuron-neuron NBLASST for all hemibrain neurons
hemibrain.nblast = hemibrain_nblast(nblast = "hemibrain")

# A normalised neuron-neuron NBLASST for all hemibrain neurons's axons
hemibrain.axon.nblast = hemibrain_nblast(nblast = "hemibrain-axons")

# And for simplified hemibrain neurons
hemibrain.simp.nblast = hemibrain_nblast(nblast = "hemibrain-simplified")

}}
}
\seealso{
\code{\link{hemibrain_splitpoints}},
  \code{\link{hemibrain_flow_centrality}},
  \code{\link{hemibrain_metrics}},
  \code{\link{hemibrainr_googledrive_data}}
}
