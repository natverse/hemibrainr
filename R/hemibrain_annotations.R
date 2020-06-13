name2cbf <- function(x) stringr::str_match(x, '([AP][DV][ML][0-9]+)')[,2]

#' Return consistent cell body fiber annotation for neurons
#'
#' @details New style cell body fiber names are of the form \code{ADL01} and
#'   match the regular expression \code{'([AP][DV][ML][0-9]+)'}. As of hemibrain
#'   v1.0.1 these seem to be variably located in the cellBodyFiber and instance
#'   fields and in some cases are missing from both.
#'
#' @param x A vector of bodyids, a query or a data.frame compatible with
#'   \code{\link{neuprint_ids}}
#' @param meta Whether to return just the name of the cell body fibers (when
#'   \code{meta=FALSE}, the default) or a data frame with additional metadata,
#'   as returned by \code{\link{neuprint_get_meta}}.
#' @param ... Additional arguments passed to \code{\link{neuprint_get_meta}}
#'
#' @return a vector of \bold{new style} cell body fiber names or a data.frame
#'   (as returned by \code{\link{neuprint_get_meta}} with an extra column
#'   \code{cbf} containing the preferred cell body fiber names).
#' @export
#' @seealso \code{\link{neuprint_get_meta}}, \code{\link{neuprint_ids}}
#' @examples
#' \donttest{
#' # these are all one lineage but some missing tract info
#' table(hemibrain_cbf("adPN"), useNA = 'ifany')
#' # Do these come from two different hemilineages?
#' table(hemibrain_cbf("lPN"))
#' }
hemibrain_cbf <- function(x, ..., meta=FALSE) {
  if(!is.data.frame(x)) {
    x <- neuprint_get_meta(x, ...)
  }
  # cbf from name should be preferred if it exists
  cbf.name=name2cbf(x$name)
  cbf.cbf=name2cbf(x$cellBodyFiber)
  x$cbf=ifelse(is.na(cbf.name), cbf.cbf, cbf.name)

  if(isTRUE(meta)) x else x$cbf
}

#' Return body ids for large classes of neurons
#'
#' Will use our google sheet annotation spreadsheets or neuPrint queries as
#' appropriate
#'
#' @param x A class of neuron (e.g. ORN, PN, DAN, MBON LHON, TOON, TOLN, HRN
#'   etc.)
#' @param ... Additional arguments passed to internal functions
#' @param unlist return all the bodyids
#' @param unique Whether to ensure that no duplicate ids are returned when
#'   \code{unlist=TRUE}.
#' @param refresh Whether to refresh cached class information from remote
#'   sources. When \code{FALSE} (the default) this is automatically refreshed
#'   every hour.
#'
#' @return When \code{unlist=FALSE} a list with elements named by the elements
#'   of \code{x} (classes of neurons). When \code{unlist=TRUE}, a vector without
#'   names of all bodyids.
#' @export
#'
#' @examples
#' # NB these examples will not be run by default as they currently depend on
#' # access to some private google sheets
#' \dontrun{
#' class2ids("uPN")
#' length(class2ids("ORN", possible = FALSE))
#' length(class2ids("ORN", possible = TRUE))
#'
#' class2ids("DAN")
#' class2ids("MBON")
#' }
class2ids <- function(x, ..., unlist=TRUE, unique=TRUE, refresh=FALSE) {
  res=list()
  aln=c("RN", "ORN", "HRN", "TRN", "PN", "uPN", "mPN", "ALLN", "VPPN")
  aln.toget=intersect(x, aln)
  if(length(aln.toget))
    res=sapply(aln.toget, alns, simplify = F, refresh=refresh, ...)
  x=setdiff(x, aln)

  mbn=c("DAN","MBON")
  mbn.toget=intersect(x, mbn)
  if(length(mbn.toget))
    res=c(res, sapply(mbn.toget, mbns, simplify = F, ...))
  x=setdiff(x, mbn)

  dnn=c("DN")
  dnn.toget=intersect(x, dnn)
  if(length(dnn.toget))
    res=c(res, sapply(dnn.toget, dnns, simplify = F, ...))
  x=setdiff(x, dnn)

  tonn=c("TON")
  tonn.toget=intersect(x, tonn)
  if(length(tonn.toget))
    res=c(res, sapply(tonn.toget, tonns, simplify = F, ...))
  x=setdiff(x, tonn)

  if(length(x)) {
    stop("Classes:", paste(x, collapse=" "), " not yet implemented")
  }

  if(unlist) {
    res=unlist(res, use.names = F)
    if(unique) res=unique(res)
  }
  res
}

aldf <- function(...)
  cached_read_sheet("124eTYqQ8evTGm_z75V8jNVmfBI763_s4h1EAPVMiSvI", ...)

alns <- function(x="RN", possible=TRUE, refresh=FALSE) {
  aldf=aldf(forget=refresh)
  if(x=='RN') {
    aldf$bodyid[which(aldf$class=="RN")]
  } else if (x=="HRN") {
    if(possible) {
      aldf$bodyid[grepl("HRN_", aldf$our_type)]
    } else {
      aldf$bodyid[grepl("^HRN_", aldf$our_type)]
    }
  } else if (x=="TRN") {
    if(possible) {
      aldf$bodyid[grepl("TRN_", aldf$our_type)]
    } else {
      aldf$bodyid[grepl("^TRN_", aldf$our_type)]
    }
  } else if (x=="ORN") {
    if(possible) {
      aldf$bodyid[grepl("ORN_", aldf$our_type)]
    } else {
      aldf$bodyid[grepl("^ORN_", aldf$our_type)]
    }
  } else if(x=="ALLN") {
    aldf$bodyid[which(aldf$class=="ALLN")]
  } else if(x=="PN") {
    aldf$bodyid[which(aldf$class=="PN")]
  } else if(x=="uPN") {
    upns=aldf$class == 'PN' &
      stringr::str_detect(aldf$our_type, "multi", negate = T) &
      stringr::str_detect(aldf$our_type, "\\+", negate = T)
    aldf$bodyid[which(upns)]
  } else if(x=="mPN") {
    mpns=aldf$class == 'PN' &
      stringr::str_detect(aldf$our_type, "multi", negate = F) |
      stringr::str_detect(aldf$our_type, "\\+", negate = F)
    aldf$bodyid[which(mpns)]
  } else if(x=="VPPN") {
    vppns=aldf$class == 'PN' &
      stringr::str_detect(aldf$our_type, "VP", negate = F)
    aldf$bodyid[which(vppns)]
  } else stop("I do not recognise: ", x, " as a class of AL neuron!")
}

#' Find the antennal lobe glomerulus for hemibrain body ids
#'
#' @param bodyids Vector of bodyids or a search string passed to
#'   \code{\link{neuprint_ids}}
#' @param exclude.multi When \code{TRUE} multiglomerular neurons will be
#'   returned as \code{NA}.
#'
#' @return A character vector of glomeruli named with the \code{bodyids}
#' @export
#'
#' @examples
#' \dontrun{
#' table(glomerulus(class2ids('PN')), useNA='ifany')
#' table(glomerulus(class2ids('PN'), exclude.multi = TRUE), useNA='ifany')
#' }
#' @importFrom neuprintr neuprint_ids
glomerulus <- function(bodyids, exclude.multi=FALSE) {

  aldf <- aldf()
  bodyids <- neuprint_ids(bodyids)

  res=character(length = length(bodyids))
  names(res)=bodyids

  rns=intersect(bodyids, class2ids("RN"))
  if(length(rns)) {
    gloms=sub(".*RN_(.*)","\\1", aldf[['our_type']][match(rns, aldf$bodyid)])
    res[rns]=gloms
  }

  pns=intersect(bodyids, class2ids("PN"))
  if(length(pns)) {
    gloms=stringr::str_match(aldf[['our_type']][match(pns, aldf$bodyid)],
                        "(.*)\\s+[a-z]+PN")[,2]
    if(exclude.multi) {
      gloms[grepl("multi", gloms)]=NA_character_
    }
    res[pns]=gloms
  }
  res
}

mbdf <- function(sheet=c("MBONs", "DANs"), ...) {
  sheet=match.arg(sheet)
  cached_read_sheet("1NvQjqt4sSZR_rqiuDZoIsuaPB1mTV9x-R-tAtfyCd2c", sheet=sheet, ...)
}

mbns <- function(x=c("MBON", "DAN"), possible=TRUE, refresh=FALSE) {
  # FIXME decide what possible means for MBONs/DANs
  x=match.arg(x)
  mbdf=mbdf(x, forget=refresh)
  if(x == "MBON")
    neuprint_ids(mbdf[mbdf$classification=='safe',])
  else
    neuprint_ids(mbdf)
}

dndf <- function(...) {
  sheet=c("dns")
  cached_read_sheet("14FtIkZV2XZHUdQFEwhDHundO0SbfWd_HKo0aJwT5kd0", sheet=sheet, ...)
}

dnns <- function(possible=TRUE, refresh=FALSE) {
  dndf=dndf(forget=refresh)
  neuprint_ids(dndf)
}

tondf <- function(sheet=c("lhns"), ...) {
  sheet=match.arg(sheet)
  cached_read_sheet("1OSlDtnR3B1LiB5cwI5x5Ql6LkZd8JOS5bBr-HTi0pOw", sheet=sheet, ...)
}

tonns <- function(x=c("lhns"), possible=TRUE, refresh=FALSE) {
 x=match.arg(x)
 tondf=tondf(x, forget=refresh)
 neuprint_ids(tondf)
}

# hidden
cached_read_sheet <- function(id, forget=FALSE, ...) {
  if(!requireNamespace('googlesheets4', quietly = TRUE))
    stop("Please install the suggested package googlesheets4!\n",
         '  install.packages("googlesheets4")')
  if(forget)
    memoise::forget(cached_read_sheet_memo)
  cached_read_sheet_memo(id, ...)
}

# hidden
cached_read_sheet_memo <- memoise::memoise(function(id, sheet=NULL) {
  googlesheets4::read_sheet(id, sheet=sheet)
}, ~memoise::timeout(3600))

