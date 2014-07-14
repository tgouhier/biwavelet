wdist <- function (wt1, wt2, cutoff=0.99) {
  wcov=Re(wt1)%*%(t(Re(wt2)))
  wsvd=svd(wcov)
  ## Cutoff point: find first value greater than cutoff (select min of 3 freqs)
  nfreqs=max(3, which(cumsum(sqrt(wsvd$d))/sum(sqrt(wsvd$d)) >= cutoff)[1])
  u=wsvd$u[1:nfreqs,]
  v=wsvd$v[1:nfreqs,]
  
  Lnk=t(u[, 1:nfreqs]) %*% Re(wt1[1:nfreqs,])
  Ljk=t(v[, 1:nfreqs]) %*% Re(wt2[1:nfreqs,])
  
  ## Distances 1
  D1=rowSums(atan(abs(
    (Lnk[,1:(NCOL(Lnk)-1)] - Ljk[,1:(NCOL(Ljk)-1)]) -
      (Lnk[,2:NCOL(Lnk)] - Ljk[,2:NCOL(Ljk)]))))
  ## Distances 2
  D2=rowSums(atan(abs(
    (u[,1:(NROW(u)-1)] - v[,1:(NROW(v)-1)]) -
      (u[,2:NROW(u)] - v[,2:NROW(v)]))))
  ## Weights based on the amount of variance explained by each axis
  w=sqrt(wsvd$d[1:nfreqs])/sum(sqrt(wsvd$d[1:nfreqs]))
  D=weighted.mean(D1+D2, w)
  return (D)
}
