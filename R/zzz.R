.onLoad <- function(libname, pkgname) {

  print("Compiling sraplus support files")
  TMB::compile(system.file("tmb", "baranov.cpp", package = "sraplus"), "-O0")

  dyn.load(TMB::dynlib(
    gsub(".cpp","",system.file("tmb", "baranov.cpp", package = "sraplus"))))

}
