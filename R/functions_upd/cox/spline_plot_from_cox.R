spline_plot_from_cox <- function(cox_model){
  original_data <- model.frame(cox_model)
  len_og <- nrow(original_data)
  odm <- original_data$`rcs(met, 4)`
  og_met <- map_dbl(odm, function(e){
    e[[1]]
  })
  og_met <- og_met[1:len_og]
    tibble(
      pred=predict(cox_model),
      met=og_met
    ) |> 
    ggplot(aes(x=met, y=pred))+
    geom_bin2d()+
    geom_smooth(color="red", size = 1.5) 
}