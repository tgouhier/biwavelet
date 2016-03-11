#' Plot phases with arrows
#' 
#' @param x x-coordinates
#' @param y y-coordinates
#' @param phases phases
#' @param arrow.len size of the arrows. Default is based on plotting region (min(par()$pin[2]/30,par()$pin[1]/40).
#' @param arrow.lwd width/thickness of arrows. Default is arrow.len * 0.3.
#' @param arrow.col arrow line color. Default is \code{black}.
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

#' Helper function for phase.plot
#' @param x TODO
#' @param y TODO
#' @param l TODO
#' @param w TODO
#' @param alpha TODO
#' @param col TODO
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
