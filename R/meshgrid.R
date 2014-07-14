meshgrid <-
function (xv, yv) {
  return(list(
    x=outer(yv * 0, xv, FUN="+"),
    y=outer(yv, xv * 0, FUN="+")
    ))
}
