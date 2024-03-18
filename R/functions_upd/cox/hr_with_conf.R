hr_with_conf <- function(cox_model){
  conf <- confint(cox_model)
  
  conf |> 
    as_tibble(rownames = NA) |> 
    rownames_to_column() |> 
    rename(lim_2_5 = `2.5 %`, lim_97_5 = `97.5 %`) |> 
    mutate(hr = (lim_2_5 + lim_97_5)/2) |> 
    mutate(across(
      c(lim_2_5, lim_97_5, hr),
      exp
    )) |> 
    filter(startsWith(rowname, "met"))
}