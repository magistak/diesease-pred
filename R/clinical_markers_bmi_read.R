library(tidyverse)
library(data.table)
library(DBI)
library(dbplyr)
library(duckdb)

con <- DBI::dbConnect(duckdb::duckdb())
file_path <- "data/phenotype/ukbiobank_20240221_phenotypic_data_decoded/ukbiobank_20240221_phenotypic_data_decoded.tsv"

duckdb_read_csv(conn = con,files =  file_path, delim = "\t", name = "all")
colnames(con |> tbl("all"))
table <- con |> tbl("all")
table |> 
  distinct(data_type_flag) |> 
  as_tibble() |> 
  View()

table |> 
  select(data_type_flag, field_title) |> 
  distinct() |> 
  as_tibble() |> 
  View()

biochem_tibble <- 
con |> 
  tbl("all") |> 
  select(eid, value, field_title, instance_index, data_type_flag) |> 
  filter((data_type_flag == "blood_biochemistry")) |> 
  mutate(value = as.numeric(value)) |> 
  select(eid, value, field_title) |> 
  pivot_wider(names_from = field_title, values_from = value) |> 
  as_tibble()

bmi_tibble <- 
  con |> 
  tbl("all") |> 
  select(eid, value, field_title, instance_index, data_type_flag) |> 
  filter((field_title == "Body mass index (BMI)")) |> 
  mutate(value = as.numeric(value)) |> 
  select(eid, value) |>
  as_tibble() |> 
  rename(bmi = value)

# write biochem tibble
write_tsv(biochem_tibble, file = "data/raw/biochem_tibble.tsv")
write_tsv(bmi_tibble, file = "data/raw/bmi_tibble.tsv")


