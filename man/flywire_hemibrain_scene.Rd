% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ngl-scenes.R
\name{flywire_hemibrain_scene}
\alias{flywire_hemibrain_scene}
\title{Create a flywire scene that shows your chosen flywire and hemibrain neuron meshes}
\usage{
flywire_hemibrain_scene(
  bodyid,
  flywire.id,
  hb.rainbow = FALSE,
  fw.rainbow = FALSE,
  ...
)
}
\arguments{
\item{bodyid}{character,  a vector of one or more valid hemibrain body IDs}

\item{flywire.id}{character, a vector of one or more valid flywire root IDs}

\item{hb.rainbow}{logical, if \code{TRUE} hemibrain neurons are shown in rainbow colours, else in shades of blue}

\item{fw.rainbow}{logical, if \code{TRUE} flywire neurons are shown in rainbow colours, else in shades of red}

\item{...}{Additional arguments passed to \code{fafbseg::ngl_encode_url}.}
}
\value{
A URL
}
\description{
Return a sample Neuroglancer scene URL for FlyWire + hemibrain  dataset
}
\examples{
\donttest{
\dontrun{
flywire.id = c("720575940611249010","720575940626212584")
bodyid = c("1540129092","1513363614")
scene = flywire_hemibrain_scene(bodyid, flywire.id)
scene
}}
}
\seealso{
\code{\link{hemibrain_read_neurons}}
}
