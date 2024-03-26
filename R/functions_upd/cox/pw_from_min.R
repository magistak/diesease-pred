pw_from_min <- function(cox_table, min, model_vars, var_main, spline_vars = character()){
  cox_table <- 
    cox_table |> 
    mutate(met_x = !!sym(var_main),
           met_y = pmax(!!sym(var_main) - min, 0))
  
  surv_obj <- with(cox_table, Surv(time = time, event = event))
  pw_vars <- c(model_vars[model_vars != var_main], "met_x", "met_y")
  formula_vars <- spline_string(pw_vars, spline_vars)
  formula <- as.formula(paste0("surv_obj ~ " , paste(formula_vars, collapse = " + ")))
  return(cph(formula, data = cox_table))
}