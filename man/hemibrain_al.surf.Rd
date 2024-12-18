% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hemibrain_data.R
\docType{data}
\name{hemibrain_al.surf}
\alias{hemibrain_al.surf}
\alias{hemibrain_al_microns.surf}
\alias{hemibrain_glomeruli_summary}
\alias{odour_scenes}
\alias{flywire_al.surf}
\title{Meshes for Antennal Lobe (AL) glomeruli}
\format{
An object of class \code{hxsurf} of length 4.

An object of class \code{hxsurf} of length 4.

An object of class \code{data.frame} with 58 rows and 11 columns.

An object of class \code{data.frame} with 83 rows and 7 columns.

An object of class \code{hxsurf} of length 4.
}
\source{
\url{https://neuprint.janelia.org/}
}
\usage{
hemibrain_al.surf

hemibrain_al_microns.surf

hemibrain_glomeruli_summary

odour_scenes

flywire_al.surf
}
\value{
\code{hemibrain_glomeruli_summary} is a \code{data.frame} with the
 columns:

 \itemize{

 \item{"expected_cit"}{the citation for the number of expected RNs},

 \item{"expected_RN_female_1h"}{.}

 \item{"expected_RN_female_SD"}{.}

 \item{"missing"}{.}

 \item{"odour_scenes"}{the general 'odour scene(s)' which this glomerulus may
 help signal (Mansourian 2015, Bates 2020)}

 \item{"key_ligand"}{the ligand that excites the cognate ALLRN or receptor
 the most, based on pooled data from multiple studies (Münch 2016)).}

 \item{"valence"}{the presumed valence of this odour channel (Badel 2016).} }
}
\description{
3D meshes of 51 olfactory + 7 thermosensory antennal lobe glomeruli for the
the hemibrain volume. They are represented as a \code{\link[nat]{hxsurf}}
object, in which 'regions' are defined, one for each glomerulus. Meshes are
available both in raw voxels and in microns. You can see which glomeruli are
truncated, and by roughly how much by examining. The glomeruli are made from PN
dendrites.
\code{hemibrain_glomeruli_summary}.
}
\details{
note that hemibrain coordinate system  has the anterior-posterior
 axis aligned with the Y axis (rather than the Z axis, which is more commonly
 observed).

 These meshes were generated by Sri Jagannathan and Tomke Stuerner based on
 the location of PN dendritic synapses.
}
\examples{
\donttest{
nvertices(hemibrain_al.surf)
materials(hemibrain_al.surf)

\dontrun{
plot3d(hemibrain_al.surf, alpha=.3)

# can also plot single glomeruli
plot3d(subset(hemibrain_al.surf,"DM1"))

}
}
}
\seealso{
\code{\link{hemibrain_roi_meshes}}, \code{\link{hemibrain.surf}}
}
\keyword{datasets}
