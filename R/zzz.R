# zzz.R
#
# Package startup and unload functions




.onLoad <- function(libname, pkgname) {
  # Install required package

  # a. For ngs.plot
  for (pkg in c("doMC", "caTools", "utils")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      utils::install.packages(pkg, dep=T)
    }
  }


  # b. For functions in this package generally
  if (!requireNamespace("BiocManager", quietly = TRUE)){
    utils::install.packages("BiocManager")
  }

  for (pkg in c("edgeR", "limma", "Rsubread", "BSgenome", "Rsamtools", "ShortRead", "rtracklayer")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      BiocManager::install(pkg, version = "3.8")
    }
  }

  if (! requireNamespace("biomaRt", quietly = TRUE)) {
    BiocManager::install("biomaRt")
  }

  invisible()
}


.onAttach <- function(libname, pkgname) {
  # Startup message
  m <- character()
  m[1] <- "\nWelcome to BCBProject.\n"

  packageStartupMessage(paste(m, collapse=""))
}


# .onUnload <- function(libname, pkgname) {
#
# }



# [END]
