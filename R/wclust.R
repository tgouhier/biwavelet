wclust <- function (w.arr) {
  s=dim(w.arr)
  nW=s[1]
  dist.matrix=matrix(NA, nrow=nW, ncol=nW)
  k=1
  nWaves=seq(from=1, to=nW, by=1)
  prog.bar=txtProgressBar(min = 0, length(nWaves)^2, style = 3)
  for (n in nWaves) {
    for (j in nWaves) {
      dist.matrix[n,j]=wdist(w.arr[n, ,], w.arr[j, ,])
      k=k+1
      setTxtProgressBar(prog.bar, k)
    }
  }
  close(prog.bar)
  return (list(diss.mat=dist.matrix, dist.mat=as.dist(dist.matrix)))
}
