spline_minimum <- function(cox_model){
  met_pred <- 
    tibble(
      met = orig_met_from_spline_cox(cox_model),
      pred = predict(cox_model)
    )
  met_seq <- seq(min(met_pred$met), max(met_pred$met), length.out = 100)
  
  model <- gam(pred ~ s(met), data = met_pred)
  
  met_seq <- seq(min(met_pred$met), max(met_pred$met), length.out = 100)
  pred_seq <- predict(model, newdata = data.frame(met = met_seq))
  
  min_index <- which.min(pred_seq)
  if(min_index <= 5 | min_index >= 95){
    met_min <- as.numeric(NA)
  } else {
  met_min <- met_seq[min_index]
  }
  
  return(met_min)
}