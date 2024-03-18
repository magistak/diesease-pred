library(tidyverse)
# read tables
met_table_read <- read_tsv("data/ukbiobank_metabolome_base_data_metabolomics_20231120.tsv") 
colnames(met_table_read)

met_table_wide <- 
met_table_read |> 
  filter(instance_index == 0 & iqr4_log == 0) |> 
  select(eid, met_name, met_value) |> 
  pivot_wider(names_from = met_name, values_from = met_value)

write_tsv(met_table_wide,"data/part_met_wide_log.tsv")

rm(met_table_read)
# 
# # create smaller data with filtering eids
# eids_sample <- sample(met_table_read$eid, 50000)
# met_table_sample <- met_table_read |> 
#   filter(eid %in% eids_sample)
# rm(met_table_read)
# write_tsv(met_table_sample, "data/sample_part_met.tsv")
# 
# # read other table
# other_table_read <- read_tsv("data/part_other.tsv")
# other_table_sample <- other_table_read |> 
#   filter(eid %in% eids_sample)
# rm(other_table_read)
# write_tsv(other_table_sample, "data/sample_part_other.tsv")

# read sample data
# other_table_sample <- read_tsv("data/sample_part_other.tsv")
# met_table_sample <- read_tsv("data/sample_part_met.tsv")
# 
# # create table for survival analysis
# met_table_sample
# # create ids for metabolite names
# met_vec <- met_table_sample |> 
#   distinct(met_name) |> 
#   unlist(use.names = F) |> 
#   sort()
# met_key <- tibble(met_name=met_vec) |> 
#   mutate(met_id=row_number())
#   
# # update table with met ids
# met_table_sample_id <- met_table_sample |> 
#   left_join(met_key) |> 
#   select(-met_name)
# 
# print(object.size(met_table_wide), units="auto")
# 
# met_table_wide <- 
# met_table_sample |> 
#   filter(instance_index==0 & iqr4 == 0) |> 
#   select(-instance_index, -iqr4, met_value_log) |> 
#   pivot_wider(id_cols = eid, names_from=met_name, values_from=met_value)

# pivot metabolite data wider
# met_table_wide <- 
# read_tsv("data/part_metabolomics.tsv") |> 
#   filter(instance_index==0 & iqr4 == 0) |> 
#   select(-instance_index, -iqr4, met_value_log) |> 
#   pivot_wider(id_cols = eid, names_from=met_name, values_from=met_value)
# write_tsv(met_table_wide, "data/part_met_wide.tsv")

met_table_wide <- read_tsv("data/part_met_wide.tsv")
other_table_sample <- read_tsv("data/sample_part_other.tsv")
met_table_sample <- read_tsv("data/sample_part_met.tsv")

met_table_sample
other_table_read <- read_tsv("data/part_other.tsv")

# filter other data
fields <- 
other_table_sample |> 
  distinct(field_title)

field_titles_other <- c("Date of attending assessment centre",
                        "Spectrometer",
                        "Sex",
                        "Year of birth",
                        "Month of birth")
other_shorter <- 
other_table_read |>
  select(eid, field_title, value, instance_index) |> 
  filter(instance_index == 0) |> 
  select(-instance_index) |> 
  filter(field_title %in% field_titles_other)
write_tsv(other_shorter, "data/other_shorter.tsv")

rm(other_table_read)

other_table_sample |> 
  filter(eid=="1000021") |> 
           View()

# diagnosis table
disease <- read_tsv("data/ukbiobank_diagnosis_table_20230922_with_mortality.txt")
disease |> 
  distinct(disease_group) |> 
  View()
disease |> 
#  filter(disease_group!="all-cause mortality") |> 
  arrange(desc(diagnosis_date)) |> 
  View()

dis_groups <- unique(disease$disease_group)
length(dis_groups)

disease |> 
  filter(disease_group %in% dis_groups[21:30]) |> 
  ggplot(aes(x=diagnosis_date)) +
  geom_histogram() +
  facet_grid(rows = vars(disease_group))
