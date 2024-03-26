load_libraries <- function(libs = c("survival", "plot")){
library(tidyverse)
library(purrr)
library(furrr)
library(tictoc)
library(janitor)
library(forcats)
  
if("survival" %in% libs){
  library(survival)
  library(rms)
  library(mgcv)
}
  
if("plot" %in% libs){
  library(patchwork)
}
}
