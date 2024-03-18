cox_stats_from_table <- function(list_table){
  
  if("sex" %in% names(list_table)){
    
  cox_df <- list_table[[1]]
  surv_obj <- with(cox_df, Surv(time = time, event = event))
  model_linear <- cph(surv_obj ~ spectrometer + rcs(age, 4) + met_log + sex, data = cox_df)
  model_spline <- cph(surv_obj ~ spectrometer + rcs(age, 4) + rcs(met_log, 4) + sex, data = cox_df)
  lrtest_res <- lrtest(model_spline, model_linear)
  return(
  setNames(list(list(lrt=lrtest_res$stats[1][[1]],
       p=lrtest_res$stats[3][[1]])),
       names(list_table))
  )
  
  } else {
    cox_df <- list_table[[1]]
    surv_obj <- with(cox_df, Surv(time = time, event = event))
    model_linear <- cph(surv_obj ~ spectrometer + rcs(age, 4) + met_log, data = cox_df)
    model_spline <- cph(surv_obj ~ spectrometer + rcs(age, 4) + rcs(met_log, 4), data = cox_df)
    lrtest_res <- lrtest(model_spline, model_linear)
    return(
    setNames(list(list(lrt=lrtest_res$stats[1][[1]],
                       p=lrtest_res$stats[3][[1]])),
             names(list_table))
    )
  }

}