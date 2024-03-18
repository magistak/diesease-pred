library(tidyverse)
library(DBI)
library(dbplyr)
library(duckdb)
library(janitor)



con <- DBI::dbConnect(duckdb::duckdb())
file_path <- "data/phenotype/ukbiobank_20240221_phenotypic_data_decoded/ukbiobank_20240221_phenotypic_data_decoded.tsv"

duckdb_read_csv(conn = con,files =  file_path, delim = "\t", name = "all")
table <- con |> tbl("all")
str(table)
colnames(table)
eids5000 <- 
  table |> 
  select(eid) |> 
  as_tibble() |> 
  distinct() |> 
  head(5000) |> 
  unlist(use.names = F)

fields <- 
table |> 
  filter(eid %in% eids5000) |> 
  select(data_type_flag, field_title,  value_type_name) |> 
  as_tibble() |> 
  distinct() 

pheno_long <- 
table |> 
  filter(instance_index == 0) |> 
  mutate(eid = as.character(eid)) |> 
  select(eid, field_title, value) |> 
  as_tibble()

pheno_long <- 
  pheno_long |> 
  pivot_wider(names_from = field_title,
              values_from = value)

pheno_long <- 
pheno_long |> 
  clean_names() 

# from pheno long I need: body_mass_index_bmi, clinical biomarkers
biochem_fields <- 
(fields |> 
  filter(data_type_flag == "blood_biochemistry"))$field_title 

make_numeric <- function(l){
  map_dbl(
    l, function(e){
      ifelse(length(e) == 0, as.numeric(NA), as.numeric(e))
    }
  )
}

pheno_out <- 
pheno_long |> 
  select(eid, all_of(make_clean_names(biochem_fields)), body_mass_index_bmi) |> 
  mutate(across(where(is.list), ~make_numeric(.x)))

pheno_out <- 
  pheno_out |> 
  rename(bmi = body_mass_index_bmi) |> 
  rename_with(~str_c(.,"_clinical"), where(is.numeric) & !any_of("bmi"))
dbDisconnect(con)
rm(con)

biochem_tibble <- 
  con |> 
  tbl("all") |> 
  select(eid, value, field_title, instance_index, data_type_flag) |> 
  filter((data_type_flag == "blood_biochemistry")) |> 
  mutate(value = as.numeric(value)) |> 
  select(eid, value, field_title) |> 
  pivot_wider(names_from = field_title, values_from = value) |> 
  as_tibble()


# surv_all
surv_all <- read_tsv("data/surv_all_log.tsv") 
health_filters <- read_tsv("data/ukb_healthy_subpopul_full_table.tsv") |> 
  mutate_all(as.character)
surv_all <- 
surv_all |> 
  clean_names() |> 
  mutate(eid = as.character(eid), sex=as.character(sex), spectrometer=as.character(spectrometer)) |> 
  left_join(pheno_out) |> 
  merge_diseases(dis_to_merge, "cvd") |> 
  mutate(non_hdl_clinical = cholesterol_clinical- hdl_cholesterol_clinical,
         non_hdl = total_cholesterol - hdl_cholesterol,
         across(
           all_of(clin_list),
           iqr4_to_na
         )
         ) |> 
  mutate(eid = as.character(eid)) |> 
  left_join(
    select(health_filters, eid, takes_medication_cholesterol, starts_with("is_healthy"))
  ) 

# add cvd, add non_hdl

dis_to_merge <- c("AAA", "Atrial Fibrillation", "Cerebral Stroke", "Heart Failure", "MACE", "PAD", "Venous Thrombosis") |> make_clean_names()


# iqr4 clinical

write_tsv(surv_all, file = "data/surv_all_bmi_clin.tsv")
rm(surv_all)








