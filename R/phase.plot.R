phase.plot <- function (x, y, phases, arrow.len=min(par()$pin[2]/30,par()$pin[1]/40), 
                        arrow.col="black", arrow.lwd=arrow.len*0.3) {
  
  a.row =seq(1,NROW(phases),round(NROW(phases)/30))
  a.col =seq(1,NCOL(phases),round(NCOL(phases)/40))
  ratio = par()$pin[2]/par()$pin[1]
  phases[-a.row,] = NA
  phases[,-a.col] = NA
  
  for (i in 1:NROW(phases)) {
    for (j in 1:NCOL(phases)) {
      if (!is.na(phases[i,j])) {	
        arrow(x[j], y[i], l=arrow.len, w = arrow.lwd, alpha = phases[i,j], 
              col=arrow.col)
      }
    }
  }
}

arrow <- function (x, y, l = 0.1, w = 0.3*l, alpha, col='black') {
  l2 <- l/3; w2 <- w/6; l3 <- l/2
  x1 <- l*cos(alpha); y1 <- l*sin(alpha)
  x2 <- w*cos(alpha+pi/2); 	y2 <- w*sin(alpha+pi/2)
  x7 <- w*cos(alpha+3*pi/2); 	y7 <- w*sin(alpha+3*pi/2)	
  x3 <- l2*cos(alpha)+w2*cos(alpha+pi/2); y3 <- l2*sin(alpha)+w2*sin(alpha+pi/2)
  x6 <- l2*cos(alpha)+w2*cos(alpha+3*pi/2); y6 <- l2*sin(alpha)+w2*sin(alpha+3*pi/2)
  x4 <- l3*cos(alpha+pi)+w2*cos(alpha+pi/2); y4 <- l3*sin(alpha+pi)+w2*sin(alpha+pi/2)
  x5 <- l3*cos(alpha+pi)+w2*cos(alpha+3*pi/2); y5 <- l3*sin(alpha+pi)+w2*sin(alpha+3*pi/2)
  X <- (par()$usr[2]-par()$usr[1])/par()$pin[1]*c(x1,x2,x3,x4,x5,x6,x7)
  Y <- (par()$usr[4]-par()$usr[3])/par()$pin[2]*c(y1,y2,y3,y4,y5,y6,y7)
  polygon(x+X, y+Y, col = col, ljoin = 1, border = NA)
}
