lrt_stats <- function(cox1, cox2){
  test <- rms::lrtest(cox1, cox2)
  test$stats[c(1,3)] |> 
    rbind() |> 
    as_tibble() |> 
    rename(lrt_stat = `L.R. Chisq`, lrt_p = P)
}