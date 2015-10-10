plot.biwavelet <- function (x, ncol=64, fill.cols=NULL, xlab="Time", ylab="Period", 
                            tol=1, plot.cb=FALSE, plot.phase=FALSE,
                            type=c("power.corr.norm", "power.corr", "power.norm", "power", "wavelet", "phase"), 
                            plot.coi=TRUE, lwd.coi=1, col.coi="white", lty.coi=1, alpha.coi=0.5,
                            plot.sig=TRUE, lwd.sig=4, col.sig="black", lty.sig=1,
                            bw=FALSE,
                            legend.loc=NULL, 
                            legend.horiz=FALSE,
                            arrow.len=min(par()$pin[2]/30,par()$pin[1]/40), 
                            arrow.lwd=arrow.len*0.3, 
                            arrow.cutoff=0.9, 
                            arrow.col="black", 
                            xlim = NULL, ylim = NULL, zlim = NULL,
                            xaxt = "s", yaxt = "s", form='%Y', ...) {
  
  if (is.null(fill.cols)) {
    if (bw) {
      fill.cols <- c("black", "white")
    }
    else {
      fill.cols <- c("#00007F", "blue", "#007FFF", 
                     "cyan","#7FFF7F", "yellow", 
                     "#FF7F00", "red", "#7F0000")
    }
  }
  
  col.pal <- colorRampPalette(fill.cols)
  fill.colors = col.pal(ncol)
  
  yrange <- ylim
  y.ticks = 2^(floor(log2(min(x$period, yrange))):floor(log2(max(x$period, yrange))))
  types=c("power.corr.norm", "power.corr", "power.norm", "power", "wavelet", "phase")
  type=match.arg(tolower(type), types)
  
  if (type=="power.corr" | type=="power.corr.norm") {
    if (x$type=="wtc" | x$type=="xwt") {
      x$power=x$power.corr
      x$wave=x$wave.corr
    }
    else
      x$power=x$power.corr
  }
  
  if (type=="power.norm" | type=="power.corr.norm") {
    if (x$type == "xwt") {
      zvals=log2(x$power)/(x$d1.sigma*x$d2.sigma)
      if (is.null(zlim))
        zlim=range(c(-1, 1) * max(zvals))
      zvals [zvals < zlim[1]]=zlim[1]
      # locs=pretty(range(zvals), n=5)
      locs=pretty(range(zlim), n=5)
      leg.lab=2^locs
    }
    else if (x$type == "wtc" | x$type == "pwtc") {
      zvals=x$rsq
      zvals[!is.finite(zvals)]=NA
      if (is.null(zlim))
        zlim=range(zvals, na.rm=TRUE)
      zvals [zvals < zlim[1]]=zlim[1]
      # locs=pretty(range(zvals, na.rm=TRUE), n=5)
      locs=pretty(range(zlim), n=5)
      leg.lab=locs
    }      
    else {
      zvals=log2(abs(x$power / x$sigma2))
      if (is.null(zlim))
        zlim=range(c(-1, 1) * max(zvals))
      zvals [zvals < zlim[1]]=zlim[1]  
      # locs=pretty(range(zvals), n=5)
      locs=pretty(range(zlim), n=5)
      leg.lab=2^locs
    }
  }
  else if (type=="power"| type=="power.corr") {
    zvals=log2(x$power)
    if (is.null(zlim))
      zlim=range(c(-1, 1)*max(zvals))
    zvals [zvals < zlim[1]]=zlim[1]
    # locs=pretty(range(zvals), n=5)
    locs=pretty(range(zlim), n=5)
    leg.lab=2^locs
  }
  else if (type=="wavelet") {
    zvals=(Re(x$wave))
    if (is.null(zlim))
      zlim=range(zvals)
    # locs=pretty(range(zvals), n=5)
    locs=pretty(range(zlim), n=5)
    leg.lab=locs
  }
  else if (type=="phase") {
    zvals=x$phase
    if (is.null(zlim))
      zlim=c(-pi, pi)
    # locs=pretty(range(zvals), n=5)
    locs=pretty(range(zlim), n=5)
    leg.lab=locs        
  }
  else {
    stop("type must be power, power.norm, power.corr, power.corr.norm, wavelet or phase")
  }
  if (is.null(xlim))
    xlim=range(x$t)
  yvals=log2(x$period)
  if (is.null(ylim))
    ylim=range(yvals)
  else
    ylim=log2(ylim)
  image(x$t,
        yvals, 
        t(zvals), 
        zlim=zlim,
        xlim=xlim,
        ylim=rev(ylim),
        xlab=xlab, 
        ylab=ylab,
        yaxt="n",
        xaxt="n",
        col=fill.colors, ...)
  box()
  if (class(x$xaxis)[1] == "Date" | class(x$xaxis)[1] == "POSIXct") {
    if (xaxt!="n") {
      xlocs=pretty(x$t)+1
      axis(side=1, at=xlocs, labels=format(x$xaxis[xlocs], form))
    }
  }
  else {
    if (xaxt!="n") {
      xlocs=axTicks(1)
      axis(side=1, at=xlocs)        
    }
  }
  if (yaxt!="n") {
    axis.locs=axTicks(2)
    yticklab=format(2^axis.locs, dig=1)
    axis(2, at=axis.locs, labels=yticklab)
  }
  
  ## Add color bar
  if (plot.cb) {
    image.plot(x$t, 
               yvals, 
               t(zvals), 
               zlim=zlim,
               ylim=rev(range(yvals)),
               xlab=xlab, 
               ylab=ylab,
               col=fill.colors,
               smallplot=legend.loc,
               horizontal=legend.horiz,
               legend.only=TRUE, 
               axis.args=
                 list(at=locs, labels=format(leg.lab, dig=2)), 
               xpd=NA)
    #box()
  }
  # COI
  if (plot.coi) {
    # lines(x$t, log2(x$coi), lty=lty.coi, lwd=lwd.coi, col=col.coi)
    polygon (x = c(x$t, rev(x$t)), lty = lty.coi, lwd = lwd.coi,
             y = c(log2(x$coi), rep(max(log2(x$coi), na.rm = TRUE), length(x$coi))), 
             col = adjustcolor(col.coi, alpha.f = alpha.coi), border = col.coi)
  }
  # sig.level contour (default is 95%)
  if (plot.sig & length(x$signif) > 1) {
    if (x$type %in% c("wt", "xwt")) {
      contour(x$t, yvals, t(x$signif), level=tol, col=col.sig, lwd=lwd.sig, 
              add=TRUE, drawlabels=FALSE)
    }
    else {
      tmp=x$rsq/x$signif
      contour(x$t, yvals, t(tmp), level=tol, col=col.sig, lwd=lwd.sig, 
              add=TRUE, drawlabels=FALSE)
    }
  }
  # Plot phases
  if (plot.phase) {
    a=x$phase
    # Remove phases where power is weak
    locs=which (zvals < quantile(zvals, arrow.cutoff))
    a[locs]=NA
    
    phase.plot(x$t, log2(x$period), a, 
               arrow.len=arrow.len, arrow.lwd=arrow.lwd, arrow.col=arrow.col)
  }
  box()
}
