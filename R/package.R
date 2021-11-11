.onUnload <- function (libpath) {
  library.dynam.unload("bayesianVARs", libpath)
}
