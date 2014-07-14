.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    paste("biwavelet", 
          utils::packageDescription("biwavelet",
                                    field="Version"),
          "loaded."),
    appendLF = TRUE)
}
