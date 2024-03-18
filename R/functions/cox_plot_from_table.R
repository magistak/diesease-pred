cox_plot_from_table <- function(list_table){
  if("sex" %in% names(list_table)){
    cox_df <- list_table[[1]]
    surv_obj <- with(cox_df, Surv(time = time, event = event))
    model_spline <- cph(surv_obj ~ spectrometer + rcs(age, 4) + rcs(met_log, 4) + sex + rcs(bmi), data = cox_df)
    p <- 
      tibble(
        pred=predict(model_spline),
        met=cox_df$met_log
      ) |> 
      ggplot(aes(x=met, y=pred))+
      geom_bin2d()+
      geom_smooth(color="red") 
      return(p)
  } else {
    cox_df <- list_table[[1]]
    surv_obj <- with(cox_df, Surv(time = time, event = event))
    model_spline <- cph(surv_obj ~ spectrometer + rcs(age, 4) + rcs(met_log, 4) + rcs(bmi), data = cox_df)
    p <- 
      tibble(
        pred=predict(model_spline),
        met=cox_df$met_log
      ) |> 
      ggplot(aes(x=met, y=pred))+
      geom_bin2d()+
      geom_smooth(color="red") 
    return(p)
    
  }
  
}
