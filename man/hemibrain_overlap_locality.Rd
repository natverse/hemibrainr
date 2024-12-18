% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hemibrain_locality.R
\name{hemibrain_overlap_locality}
\alias{hemibrain_overlap_locality}
\title{Calculate the overlap score between a neurons axon dendrite}
\usage{
hemibrain_overlap_locality(x, resample = 10, delta = 62.5, ...)
}
\arguments{
\item{x}{a \code{nat::neuronlist} or \code{nat::neuron} object. It is assumed
that this neuron has been read in by
\code{neuprintr::neuprint_read_neurons} or possibly
\code{catmaid::read.neurons.catmaid}.}

\item{resample}{stepsize to which to resample neurons. If set to \code{NULL},
neurons are not resampled.}

\item{delta}{the distance at which a synapse might occur}

\item{...}{Additional arguments passed to methods or eventually to
\code{nat::\link{nlapply}}}
}
\value{
a named vector of overlap_locality scores
}
\description{
Calculates an overlap score using an overlap score between a neuron's axon and dendrite,
if axon and dendrite have been marked out in the neuron object, e.g. by using \code{flow_centrality}.
}
\seealso{
\code{\link{hemibrain_compartment_metrics}}
}
