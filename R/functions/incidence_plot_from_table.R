incidence_plot_from_table <- function(list_table){
  p <- list_table[[1]] |> 
    arrange(met_log) |> 
    mutate(group = cut(1:nrow(list_table[[1]]), 10, labels = FALSE)) |> 
    group_by(group) |> 
    summarise(events = sum(event)) |> 
    ggplot(aes(x=factor(group), y=events))+
    geom_col()
  return(p)
}
