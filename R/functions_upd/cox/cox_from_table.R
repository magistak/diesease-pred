cox_from_table <- function(cox_table, model_vars, spline_vars=character(), spline_df=4){
  surv_obj <- with(cox_table, Surv(time = time, event = event))
  formula_vars <- spline_string(model_vars, spline_vars)
  formula <- as.formula(paste0("surv_obj ~ " , paste(formula_vars, collapse = " + ")))
  cph(formula, data = cox_table)
}