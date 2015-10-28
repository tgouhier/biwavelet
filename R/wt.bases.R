wt.bases <-
function (mother = c("morlet", "paul", "dog"), k, scale, param = -1) {

  mothers=c("morlet", "paul", "dog")
  mother=match.arg(tolower(mother), mothers)
  
  n=length(k)
  if (mother == "morlet") {
    if (param == -1) {
      param = 6
    }
    k0 = param
    expnt = -(scale * k - k0)^2/2 * (k > 0)
    norm = sqrt(scale * k[2]) * (pi^(-0.25)) * sqrt(n)
    daughter = norm * exp(expnt)
    daughter = daughter * (k > 0)
    fourier.factor = (4 * pi) / (k0 + sqrt(2 + k0^2))
    coi = fourier.factor / sqrt(2)
    dof = 2
  }
  else if (mother == "paul") {
    if (param == -1) {
      param = 4
    }
    m = param
    expnt = -(scale * k) * (k > 0)
    norm = sqrt(scale * k[2]) * (2^m / sqrt(m * prod( 2:(2*m - 1) ))) * sqrt(n)
    daughter = norm * ((scale * k)^m) * exp(expnt)
    daughter = daughter * (k > 0)
    fourier.factor = (4 * pi) / (2*m + 1)
    coi = fourier.factor * sqrt(2)
    dof = 2
  }
  else if (mother == "dog") {
    if (param == -1) {
      param = 2
    }
    m = param
    expnt = -(scale * k)^2 / 2
    norm = sqrt(scale * k[2] / gamma(m + 0.5)) * sqrt(n)
    daughter = -norm * (1i^m) * ((scale * k)^m) * exp(expnt)
    fourier.factor = 2 * pi * sqrt(2 / (2*m + 1))
    coi = fourier.factor / sqrt(2)
    dof = 1
  }
  else {
    stop("mother wavelet parameter must be 'morlet', 'paul', or 'dog'")
  }

  return (list(daughter=daughter, 
               fourier.factor=fourier.factor, 
               coi=coi, 
               dof=dof))
}
