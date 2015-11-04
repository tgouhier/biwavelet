#' Rectangular grid in two dimensions
#' 
#' Replicates the grid vectors \code{xv} and \code{yv} to generate a full grid
#' @author Tarik C. Gouhier (tarik.gouhier@@gmail.com)
#' 
#' @param xv vector of numeric values
#' @param yv vector of numeric values
#' 
#' @return Returns a list containing the full grid components \code{xv} and \code{yv}:
#' \item{x}{replicated values of xv}
#' \item{y}{replicated values of yv}
#' 
#' @examples
#' xv <- runif(10)
#' yv <- runif(5)
#' g <- meshgrid(xv, yv)
meshgrid <- function (xv, yv) {
  list(
    x = outer(yv * 0, xv, FUN="+"),
    y = outer(yv, xv * 0, FUN="+")
  )
}
