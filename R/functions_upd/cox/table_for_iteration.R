table_for_iteration <- function(table_all, table_base, met_i, dis_i){
  table_base %>% 
    mutate(met_cox = table_all[[met_i]],
           event = table_all[[paste("event", dis_i, sep = "_")]],
           time = table_all[[paste("time", dis_i, sep = "_")]])
}