orig_met_from_spline_cox <- function(cox_model){
  original_data <- model.frame(cox_model)
  len_og <- nrow(original_data)
  odm <- original_data$`rcs(met, 4)`
  og_met <- map_dbl(odm, function(e){
    e[[1]]
  })
  og_met[1:len_og]
}