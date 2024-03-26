model_vars_update <- function(model_vars, disease, disease_sex){
  # if diseases sex:
  # 0: women
  # 1: men
  # 2: both
  sex_id <- disease_sex$sex[disease_sex$disease == disease]
  if(sex_id == 2){
    return(model_vars)
  } else {
    return(model_vars[model_vars != "sex"])
  }
  
}