merge_diseases <- function(surv_table, disease_vector, merged_name){
  dis_merge_table <- 
    surv_table |> 
    select(any_of(contains(disease_vector))) 
  
  time_table <- 
    dis_merge_table |> 
    select(any_of(contains("time")))
  
  event_table <- 
    dis_merge_table |> 
    select(any_of(contains("event")))
  
  t <- 
    dis_merge_table |> 
    mutate(
      !!paste("time", merged_name, sep = "_") :=
        as.matrix(time_table) |> 
        matrixStats::rowMins(),
      !!paste("event", merged_name, sep = "_") :=
        as.matrix(event_table) |> 
        matrixStats::rowMaxs()
    ) 
  cbind(surv_table, t[ ,c(ncol(t)-1, ncol(t))]) |> 
    as_tibble()
}