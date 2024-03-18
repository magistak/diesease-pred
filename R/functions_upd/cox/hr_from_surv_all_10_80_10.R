hr_from_surv_all_10_80_10 <- function(surv_all, model_vars, eid_filter, metabolite, disease){
  surv_all <- surv_all |> 
    filter(eid %in% eid_filter)
  cwm <- as.character(unique(surv_all$sex))
  disease_sex <- ifelse(!"0" %in% cwm, "men", ifelse(!"1" %in% cwm, "women", "both"))
  if("sex" %in% model_vars & (disease_sex == "women" | disease_sex == "men")){
    model_vars <- model_vars[model_vars != "sex"]
  }
  
  cox_df <- 
    filter_surv_all(surv_all = surv_all,
                    metabolite = metabolite,
                    disease = disease,
                    eids = with(surv_all, eid),
                    model_vars = model_vars) |> 
    mutate(met = discret_by_people(met, c(0.1, 0.8, 0.1))) 
  
  # relevel
  levels <- levels(cox_df$met)
  cox_df <- cox_df |> 
    mutate(met = relevel(met, ref = levels[2]))
  
  cox_model <-cox_from_table(cox_table = cox_df,
                             model_vars = model_vars,
  )
  
  hr_with_conf(cox_model)
}