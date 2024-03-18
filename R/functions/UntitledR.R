filter_surv_all <- function(surv_all, metabolite, disease, eids){
  surv_all |>  
    filter(eid %in% eids) |> 
    mutate(met_log = log10(table_all[[metabolite]]),
           event = table_all[[paste("event", disease, sep = "_")]],
           time = table_all[[paste("time", disease, sep = "_")]])
}