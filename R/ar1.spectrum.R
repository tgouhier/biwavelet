ar1.spectrum <-
function (ar1, periods) {
  p=(1-ar1^2)/abs(1-ar1*exp(-2*1i*pi*(1/periods)))^2
  return (p)
}

