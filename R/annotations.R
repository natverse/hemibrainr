#' Return body ids for large classes of neurons
#'
#' Will use our google sheet annotation spreadsheets or neuPrint queries as
#' appropriate
#'
#' @param x A class of neuron (e.g. ORN, PN, DAN, MBON LHON, TOON, TOLN, HRN
#'   etc.)
#' @param ...
#' @param unlist return all the bodyids
#' @param unique
#'
#' @return When \code{unlist=FALSE} a list with elements named by the elements
#'   of \code{x} (classes of neurons). When \code{unlist=TRUE}, a vector without
#'   names of all bodyids.
#' @export
#'
#' @examples
#' class2ids("uPN")
#' length(class2ids("ORN", possible = FALSE))
#' length(class2ids("ORN", possible = TRUE))
class2ids <- function(x, ..., unlist=TRUE, unique=TRUE) {
  res=list()
  aln=c("RN", "ORN", "HRN", "TRN", "PN", "uPN", "mPN", "ALLN", "VPPN")
  aln.toget=intersect(x, aln)
  if(length(aln.toget))
    res=sapply(aln.toget, alns, simplify = F, ...)
  x=setdiff(x, aln)

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

alns <- function(x="RN", possible=TRUE) {
  aldf=aldf()
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
  }
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
#' \donttest{
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
    gloms=sub(".*RN_(.*)","\\1", aldf[match(rns, aldf$bodyid),'our_type'])
    res[rns]==gloms
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

cached_read_sheet <- function(id, forget=FALSE) {
  if(forget)
    memoise::forget(cached_read_sheet_memo)
  cached_read_sheet_memo(id)
}

cached_read_sheet_memo <- memoise::memoise(function(id) {
  googlesheets4::read_sheet(id)
}, ~memoise::timeout(3600))

