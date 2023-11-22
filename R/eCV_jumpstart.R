.onAttach <- function(libname, pkgname) {
  pkg_message <- paste0("

 __________________________________________________________________
|\t                     ______     __                         |
|\t              ___   / ___\\ \\   / /                         |
|\t             / _ \\ | |    \\ \\ / /                          |
|\t            |  __/ | |___  \\ V /                           |
|\t             \\___|  \\____|  \\_/                            |
|__________________________________________________________________|

        Enhanced Coefficient of Variation and IDR Extensions 
                  for Reproducibility Assessment

This package provides extensions and alternative methods to  IDR  to 
measure the reproducibility of omic data with an arbitrary number of 
replicates. It introduces an enhanced Coefficient of Variation (eCV)
metric to assess the likelihood of omic features being reproducible. 


")
  packageStartupMessage(pkg_message)
}
