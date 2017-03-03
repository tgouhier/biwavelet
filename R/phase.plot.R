#' Plot phases with arrows
#'
#' @param x X-coordinates
#' @param y Y-coordinates
#' @param phases Phases
#' @param arrow.len Size of the arrows. Default is based on plotting region.
#' @param arrow.lwd Width/thickness of arrows.
#' @param arrow.col Arrow line color.
#'
#' @author Tarik C. Gouhier (tarik.gouhier@@gmail.com)
#'
#' Huidong Tian provided a much better implementation of the phase.plot function
#' that allows for more accurate phase arrows.
#'
#' Original code based on WTC MATLAB package written by Aslak Grinsted.
#'
#' @note
#' Arrows pointing to the right mean that \code{x} and \code{y} are in phase.
#'
#' Arrows pointing to the left mean that \code{x} and \code{y} are in anti-phase.
#'
#' Arrows pointing up mean that \code{y} leads \code{x} by \eqn{\pi/2}.
#'
#' Arrows pointing down mean that \code{x} leads \code{y} by \eqn{\pi/2}.
#'
#' @examples
#' # Not run: phase.plot(x, y, phases)
#'
#' @export
phase.plot <- function(x, y, phases,
                       arrow.len = min(par()$pin[2] / 30, par()$pin[1] / 40),
                       arrow.col = "black",
                       arrow.lwd = arrow.len * 0.3) {

  a.row <- seq(1, NROW(phases), round(NROW(phases) / 30))
  a.col <- seq(1, NCOL(phases), round(NCOL(phases) / 40))

  phases[-a.row, ] <- NA
  phases[, -a.col] <- NA

  for (i in seq_len(NROW(phases))) {
    for (j in seq_len(NCOL(phases))) {
      if (!is.na(phases[i, j])) {
        arrow(x[j], y[i], l = arrow.len, w = arrow.lwd,
              alpha = phases[i, j], col = arrow.col)
      }
    }
  }
}

#' Helper function for \code{\link{phase.plot}} (not exported)
#' @param x X-coordinate of the arrow.
#' @param y Y-coordinate of the arrow.
#' @param l Length of the arrow.
#' @param w Width of the arrow.
#' @param alpha Angle of the arrow in radians (0 .. 2*pi).
#' @param col Color of the arrow.
#'
#' @examples
#' plot.new()
#' arrow(0,0, alpha = 0)
arrow <- function(x, y, l = 0.1, w = 0.3 * l, alpha, col = "black") {
  l2 <- l / 3
  w2 <- w / 6
  l3 <- l / 2
  x1 <- l * cos(alpha)
  y1 <- l * sin(alpha)
  x2 <- w * cos(alpha + pi / 2)
  y2 <- w * sin(alpha + pi / 2)
  x7 <- w * cos(alpha + 3 * pi / 2)
  y7 <- w * sin(alpha + 3 * pi / 2)
  x3 <- l2 * cos(alpha) + w2 * cos(alpha + pi / 2)
  y3 <- l2 * sin(alpha) + w2 * sin(alpha + pi / 2)
  x6 <- l2 * cos(alpha) + w2 * cos(alpha + 3 * pi / 2)
  y6 <- l2 * sin(alpha) + w2 * sin(alpha + 3 * pi / 2)
  x4 <- l3 * cos(alpha + pi) + w2 * cos(alpha + pi / 2)
  y4 <- l3 * sin(alpha + pi) + w2 * sin(alpha + pi / 2)
  x5 <- l3 * cos(alpha + pi) + w2 * cos(alpha + 3 * pi / 2)
  y5 <- l3 * sin(alpha + pi) + w2 * sin(alpha + 3 * pi / 2)
  X <- (par()$usr[2] - par()$usr[1]) / par()$pin[1] * c(x1,x2,x3,x4,x5,x6,x7)
  Y <- (par()$usr[4] - par()$usr[3]) / par()$pin[2] * c(y1,y2,y3,y4,y5,y6,y7)
  polygon(x + X, y + Y, col = col, ljoin = 1, border = NA)
}

#' This is an alternative helper function that plots arrows.
#' It uses \code{\link{text}()} to print a character using a default font.
#' This way, it is possible to render different types of arrows.
#'
#' @author Viliam Simko
#'
#' @param x X-coordinate of the arrow.
#' @param y Y-coordinate of the arrow.
#' @param angle Angle in radians.
#' @param size Similar to \code{arrow.len} parameter. Notice that we don't need
#'   the \code{arrow.lwd} anymore
#' @param col Color of the arrow.
#' @param chr Character representing the arrow. You should provide the character
#'   as escaped UTF-8.
#' @importFrom graphics text
#' @examples
#' # Not run: arrow2(x[j], y[i], angle = phases[i, j],
#' # Not run:        col = arrow.col, size = arrow.len)
arrow2 <- function(x, y, angle, size = .1, col = "black",
                   chr =  intToUtf8(0x279B)) {
  # speed optimized: 180/pi =~= 57.29578
  # note: size is 10x smaller to be compatible with the old implementation
  text(x,y, labels = chr, col = col, cex = 10 * size, srt = 57.29578 * angle)
}
