filter_surv_all <- function(surv_all, metabolite, disease, eids, model_vars, log_trans = TRUE){
  
  vars <- c("eid", "event", "time", model_vars)
  
  table <- 
    surv_all |>  
    mutate(met = ifelse(rep(log_trans,nrow(surv_all)), log10(surv_all[[metabolite]]), surv_all[[metabolite]]),
           event = surv_all[[paste("event", disease, sep = "_")]],
           time = surv_all[[paste("time", disease, sep = "_")]]) |> 
    filter(eid %in% eids) |> 
    select(all_of(vars)) |> 
    drop_na()
  
  return(
  table
  )
}