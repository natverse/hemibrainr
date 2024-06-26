% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/neuron_side.R
\name{soma_side}
\alias{soma_side}
\alias{soma_side.neuronlist}
\alias{soma_side.neuron}
\alias{synapse_side_index}
\alias{synapse_side_index.neuronlist}
\alias{synapse_side_index.neuron}
\title{Calculate which side of the brain a neuron is on}
\usage{
soma_side(
  x,
  flywire_nuclei = FALSE,
  brain = nat.flybrains::JRC2018F,
  bound = NULL,
  transform = FALSE
)

\method{soma_side}{neuronlist}(
  x,
  flywire_nuclei = FALSE,
  brain = nat.flybrains::JRC2018F,
  bound = NULL,
  transform = FALSE
)

\method{soma_side}{neuron}(
  x,
  nuclei = NULL,
  brain = nat.flybrains::JRC2018F,
  bound = NULL,
  transform = FALSE
)

synapse_side_index(
  x,
  brain = nat.flybrains::JRC2018F,
  bound = NULL,
  transform = FALSE
)

\method{synapse_side_index}{neuronlist}(
  x,
  brain = nat.flybrains::JRC2018F,
  bound = NULL,
  transform = FALSE
)

\method{synapse_side_index}{neuron}(
  x,
  brain = nat.flybrains::JRC2018F,
  bound = NULL,
  transform = FALSE
)
}
\arguments{
\item{x}{a \code{nat::neuronlist} or \code{nat::neuron} object}

\item{brain}{wthe brainspace what \code{x} is in.}

\item{bound}{the dividing line that separates left and right in e \code{JRC2018F}.}

\item{transform}{whether to transform \code{x} (when \code{TRUE}) when the given brain is not \code{JRC2018F}}
}
\value{
a \code{nat::neuronlist} object with either 'side' or synapse_side_index updated in its metadata.
}
\description{
These functions make use of the symmetrised brainspace \code{JRC2018F} to determine the 'side' of the brain on which a neuron's root node, typically the soma, falls, i.e. \code{soma_side}.
You can also calculate a 'synapse side index', defined as:  (number of synapses on the right - number on the left) / total synapses.
}
\examples{
\donttest{
library(hemibrainr)
library(nat.jrcbrains)

# Read a neuron from the hemibrain
neurons = neuprint_read_neurons(1796818119)

# Split into axon and dendrite
neurons.flow = flow_centrality(neurons)

# Calculate soma side
neurons.flow = soma_side(neurons.flow, brain = "JRCFIB2018F")

# Calculate synapse side index
neurons.flow = synapse_side_index(neurons, brain = "JRCFIB2018F")

# See side data
neurons.flow[,]

}
}
\seealso{
\code{\link{hemibrain_reroot}}
}
