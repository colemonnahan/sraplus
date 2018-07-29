#' @param fit A returned list from run.SIR
#' @param ... A series of fits as returned by run.SIR
#' @param names An optional character string with length equal to the
#'   number of fits, used to identify fits on the plot.
#' @param plot Boolean for whether to print the plot.
#' @return Functions using base plots return nothing, but those using
#'   ggplot2 return an invisible ggplot object which can then be used to
#'   print or save as normally done. In all cases if 'plot=TRUE' then a
#'   plot is created.
