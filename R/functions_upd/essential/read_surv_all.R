read_surv_all <- function(path, dis_list, nmr_list, clin_list, health_filter){
  cols_dis <- c(paste0("event_",dis_list),paste0("time_",dis_list))
  read_tsv(path) |> 
    mutate(eid = as.character(eid), sex=as.character(sex), spectrometer=as.character(spectrometer)) |> 
    select(eid, sex, spectrometer, age, bmi,
         all_of(c(cols_dis, nmr_list, clin_list, health_filter)))
}