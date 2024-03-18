load_libraries <- function(libs = c("survival", "plot")){
library(tidyverse)
library(purrr)
library(furrr)
library(tictoc)
library(janitor)
  
if("survival" %in% libs){
  library(survival)
  library(rms)
}
  
if("plot" %in% libs){
  library(patchwork)
}
}