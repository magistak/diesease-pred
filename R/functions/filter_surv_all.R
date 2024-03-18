filter_surv_all <- function(surv_all, metabolite, disease, eids, disease_sex){
  if(disease_sex == "both"){
  table <- 
  surv_all |>  
    mutate(met_log = log10(surv_all[[metabolite]]),
           event = surv_all[[paste("event", disease, sep = "_")]],
           time = surv_all[[paste("time", disease, sep = "_")]]) |> 
    filter(eid %in% eids[[1]]) |> 
    select(eid, sex, spectrometer, age, met_log, event, time, bmi) |> 
    drop_na()
  return(
    setNames(list(table),
             names(eids))
    # paste(metabolite, disease, names(eids), sep = "_"))
  )
  } 
  if(disease_sex == "men"){
    table <- 
      surv_all |>  
      mutate(met_log = log10(surv_all[[metabolite]]),
             event = surv_all[[paste("event", disease, sep = "_")]],
             time = surv_all[[paste("time", disease, sep = "_")]]) |> 
      filter(eid %in% eids[[1]],
             sex == 1) |> 
      select(eid, spectrometer, age, met_log, event, time, bmi) |> 
      drop_na()
    return(
      setNames(list(table),
               names(eids))
      # paste(metabolite, disease, names(eids), sep = "_"))
    )
  }
  if(disease_sex == "women"){
    table <- 
      surv_all |>  
      mutate(met_log = log10(surv_all[[metabolite]]),
             event = surv_all[[paste("event", disease, sep = "_")]],
             time = surv_all[[paste("time", disease, sep = "_")]]) |> 
      filter(eid %in% eids[[1]],
             sex == 0) |> 
      select(eid, spectrometer, age, met_log, event, time, bmi) |> 
      drop_na()
    return(
    setNames(list(table),
             names(eids))
   # paste(metabolite, disease, names(eids), sep = "_"))
    )
  }
}