#' @importFrom bayestestR distribution_gamma distribution_beta hdi
#' @importFrom prodlim row.match
#' @importFrom tidyr expand_grid
#' @importFrom dplyr mutate if_else
#' @importFrom reshape2 melt
#' @importFrom ggpubr ggarrange
#' @import ggplot2
#' @import ggrepel
#' @import methods
#' @import ceg

NULL

#Below are what I think are the proper imports, but couldn't get to work without loading ceg (even importing ceg doesn't seem to work)
#@importClassesFrom ceg Stratified.event.tree Stratified.staged.tree Stratified.ceg.model
#@importMethodsFrom ceg plot
