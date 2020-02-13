#' Search neuprint for body ids
#'
#' @param x A set of bodyids or a query
#' @param field Which field to query
#' @inheritParams neuprintr::neuprint_fetch_custom
#' @param exact whether to insist that the query string *exactly* matches the
#'   whole of a field.
#' @param regex whether to treat x as a regular expression (default FALSE)
#'
#' @return character vector of bodyids
#' @export
#' @seealso \code{\link[neuprintr]{neuprint_search}}
#'
#' @examples
#' \donttest{
#' neuprint_ids("MBON")
#' neuprint_ids("MBON01", exact=TRUE)
#' neuprint_ids("MBON[0-9]+", regex=TRUE)
#' }
neuprint_ids <- function(x, field='type', exact=FALSE, regex=FALSE, conn=conn, ...) {
  x <- if(is.character(x) && length(x)==1 && !valid_id(x)){
    if(isFALSE(regex))
      x <- if(exact) paste0("^",x,"$") else paste0(".*",x,".*")
    neuprintr::neuprint_search(x, meta = F, field = field, ...)
  } else {
    if (is.data.frame(x)) {
      nx=tolower(names(x))
      if('bodyid' %in% nx) x[[match('bodyid', nx)]]
        else stop("No bodyid column in `x`!")
    } else {
      # directly specified
      x
    }
  }
  x <- as.character(x)
  # if(!length(x))
  #   stop("No valid ids!")
  if(!all(valid_id(x)))
    stop("Some ids are invalid!")
  x
}

valid_id <- function(x) {
  if(!length(x)) return(FALSE)
  if(is.factor(x)) x <- as.character(x)
  bx=bit64::as.integer64(x)
  # unfortunately really big integers overflow to the biggest possible int64
  biggestid=bit64::as.integer64('9223372036854775807')
  is.finite(bx) & bx>0 & bx < biggestid
}
