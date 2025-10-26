# 00_load_libraries.R

load_libraries <- function() {
  if (!require("pacman")) install.packages("pacman")
  
  pacman::p_load(
    reshape, 
    ggplot2, 
    plotly, 
    sigmoid, 
    Rlab, 
    Rmpfr, 
    cowplot, 
    lavaan, 
    torch, 
    reticulate, 
    data.table, 
    dplyr, 
    tidyverse, 
    rstudioapi,
    here      
  )
}