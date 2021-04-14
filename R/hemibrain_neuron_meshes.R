# convert cloudvolume python mesh to an R mesh3d object
# method vf just uses the vertex and face arrays
# ply writes out to Stanford ply format and reads back in again
cvmesh2mesh <- function(x, method=c('vf', 'ply'), ...) {
  method=match.arg(method)
  if(method=='vf') {
    verts=t(x$vertices)
    stopifnot(nrow(verts)==3)
    faces=t(x$faces+1)
    stopifnot(max(faces)<=ncol(verts))
    stopifnot(min(faces)>0)
    stopifnot(nrow(faces)==3)
    m=rgl::tmesh3d(vertices = verts, indices = faces, homogeneous = F)
  } else {
    bytes=x$to_ply()
    tf=tempfile(fileext = paste0('.', method))
    on.exit(unlink(tf))
    writeBin(bytes, con = tf)
    m=Rvcg::vcgPlyRead(tf, ...)
  }
  m
}

#' Read hemibrain neurite or mitochondrial meshes into an R neuronlist object
#'
#' @description \code{hemibrain_neuron_meshes} reads standard neuron meshes
#'
#' @param bodyids One more hemibrain body ids (including as expressions passed
#'   to \code{\link{neuprint_ids}})
#' @param lod the level of detail (highest resolution is 0, default of 2 gives a
#'   good overall morphology while 3 is also useful and smaller still)
#' @param cloudvolume.url Optional a specific URL from which to download mesh
#'   data.
#' @param conn A neuprint connection to the hemibrain dataset
#' @param dataset Optionally specify which version of the hemibrain dataset to
#'   use (when \code{NULL} uses the most recent)
#' @param ... Additional arguments passed to \code{\link{nlapply}}
#'
#' @return A \code{\link{neuronlist}} object
#' @export
#'
#' @examples
#' \donttest{
#' da2.lod3=hemibrain_neuron_meshes("DA2_lPN", lod=3)
#' da2.lod3
#' }
#' \dontrun{
#' da2=hemibrain_neuron_meshes("DA2_lPN")
#' # lod 2 ~ 5x bigger than lod 3 but does look nicer
#' da2
#' plot3d(da2)
#'
#' da2mito=hemibrain_mito_meshes("DA2_lPN")
#'
#' # if you just wanted one neuron (as a mesh3d object) rather than neuronlist
#' dl4=hemibrain_neuron_meshes("name:DL4_adPN_R")[[1]]
#' }
#' @importFrom fafbseg flywire_cloudvolume
hemibrain_neuron_meshes <- function(bodyids, lod=2L, cloudvolume.url=NULL, conn=NULL, dataset=NULL, ...) {
  dataset = neuprintr:::check_dataset(dataset, conn = conn)
  bodyids = neuprintr::neuprint_ids(bodyids, conn = conn, dataset = dataset)
  if(is.null(cloudvolume.url)) {
    cloudvolume.url = sprintf(
      "precomputed://gs://neuroglancer-janelia-flyem-%s/segmentation",
      sub(":", "/", dataset, fixed = T)
    )
  }
  # in case they are integer64
  bodyids=as.character(bodyids)
  nn=bodyids;names(nn) <- bodyids
  res=nlapply(nn, hemibrain_neuron_mesh, lod=lod, cloudvolume.url=cloudvolume.url, ...)
  res[,] <- neuprint_get_meta(names(res), conn=conn, dataset = dataset)
  return(res)
}

hemibrain_neuron_mesh <- function(bodyids, lod, cloudvolume.url, method='vf') {
  cv=fafbseg::flywire_cloudvolume(cloudvolume.url = cloudvolume.url)
  # convert to python longint if necessary
  check_package_available("reticulate")
  li=reticulate::py_eval(bodyids, convert = F)
  lod=as.integer(lod)
  cm=cv$mesh$get(li, lod = lod)
  if(!any(bodyids %in% names(cm))) {
    stop("Failed to read segid: ", bodyids)
  }
  # we get back a dict where the key is a python 64 bit int
  cmesh=reticulate::py_get_item(cm, li)
  m=cvmesh2mesh(cmesh, method=method)
  m
}


#' @rdname hemibrain_neuron_meshes
#' @description \code{hemibrain_mito_meshes} reads mitochondrial meshes
hemibrain_mito_meshes <- function(bodyids, lod=2L, ...) {
  cvu="precomputed://gs://neuroglancer-janelia-flyem-hemibrain/v1.2/mito-objects-grouped"
  hemibrain_neuron_meshes(bodyids, lod=lod, cloudvolume.url = cvu, ...)
}

