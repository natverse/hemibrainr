################################################################################
################################ Colours #######################################
################################################################################

#' Some nice colours used by Schlegel and Bates et al. 2020 (in prep.) and Li et
#' al. 2020 (in prep)
#'
#' @description \code{hemibrain_colours} (also \code{hemibrain_blender_colours})
#'   is a 6 colour palette from Li et al. 2020 (in prep). Although this has been
#'   described as colour blind friendly it is really not very friendly for green
#'   colour blind or green anomalous individuals and was actually designed to
#'   achieve similar perceptual brightness in Blender renderings.
#'
#'   \code{hemibrain_fig_colours} is an 8 colour palette that has been used for
#'   some of the figures. This actually works better for colour blind
#'   individiuals.
#'
#'   \code{hemibrain_bright_colours} are nice brighter colours inspired by the
#'   LaCroix colour palette.
#'
#' @param n the number of colours to generate
#' @references \url{https://en.wikipedia.org/wiki/Color_blindness}
#'
#' @source \url{https://github.com/johannesbjork/LaCroixColoR}
#'
#' @return a named vector of colours
#' @examples
#' \donttest{
#'
#' # Li et al dark colours (Blender)
#' par(mar = rep(0, 4))
#' pie(rep(1, length(hemibrain_colours)), col = hemibrain_colours)
#'
#' # Li et al figure colours
#' par(mar = rep(0, 4))
#' pie(rep(1, length(hemibrain_fig_colours)), col = hemibrain_fig_colours)
#'
#' # Bright colours
#' par(mar = rep(0, 4))
#' pie(rep(1, length(hemibrain_bright_colours)), col = hemibrain_bright_colours)
#'
#' }
#' @name hemibrain_colours
#' @export
#' @importFrom grDevices colorRampPalette
hemibrain_blender_colours <-
hemibrain_blender_colors <-
hemibrain_colours <-
hemibrain_colors <-
  c(
    red = "#A53600",
    magenta = "#B32DB5",
    midblue = "#0072B2",
    darkgold = "#908827",
    green = "#348E53",
    blue = "#053CFF"
  )

#' @rdname hemibrain_colours
hemibrain_fig_colours <-
hemibrain_fig_colors <-
  c(
    "#0072b2",
    "#56b4e9",
    "#a3cce9",
    "#cc79a7",
    "#009e73",
    "#d55e00",
    "#e69f00",
    "#f0e442"
  )

#' @rdname hemibrain_colours
#' @export
hemibrain_colour_ramp = grDevices::colorRampPalette(hemibrain_colours)

#' @rdname hemibrain_colours
#' @export
hemibrain_bright_colours <-
  c(
    purple = "#C70E7B",
    pink = "#FC6882",
    blue = "#007BC3",
    cyan = "#54BCD1",
    darkorange = "#EF7C12",
    paleorange = "#F4B95A",
    darkgreen = "#009F3F",
    green = "#8FDA04",
    brown = "#AF6125",
    palebrown = "#F4E3C7",
    mauve = "#B25D91",
    lightpink = "#EFC7E6",
    orange = "#EF7C12",
    midorange = "#F4B95A",
    darkred = "#C23A4B",
    darkyellow = "#FBBB48",
    yellow = "#EFEF46",
    palegreen = "#31D64D",
    navy = "#132157",
    cerise = "#EE4244",
    red = "#D72000",
    marine = "#1BB6AF"
  )

#' @rdname hemibrain_colours
#' @export
hemibrain_bright_colors <- hemibrain_bright_colours

#' @rdname hemibrain_colours
#' @export
hemibrain_bright_colour_ramp = grDevices::colorRampPalette(hemibrain_bright_colours)
