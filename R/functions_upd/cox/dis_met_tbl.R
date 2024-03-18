dis_met_tbl <- function(dis_list, nmr_list, clin_list){
  def_clinical <- c(nmr_list, clin_list)
  expand_grid(dis_list, def_clinical) |> 
    rename(disease = dis_list, met = def_clinical) |> 
    mutate(dm = paste(disease, met, sep = "_")) |> 
    mutate(nmr_clin = ifelse(met %in% nmr_list, "nmr", "clinical"))
}