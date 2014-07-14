phase.plot <-
function (x, y, phases, arrow.size=0.05, arrow.size.head=0.05, arrow.lwd=2, arrow.col="black") {
  g=meshgrid(x, y)
  x=g$x
  y=g$y
  asp<-par("pin")[2]/par("pin")[1] 
  dx=cos(phases)*arrow.size
  dy=sin(phases)*arrow.size*asp

  suppressWarnings(
    arrows(x, y, x + dx, y + dy, length=arrow.size.head, lwd=arrow.lwd, 
           col=arrow.col, code=2))
}
