library(tidyverse)
library(rms)
library(Hmisc)
library(survival)


surv_all <- read_tsv("data/surv_all_log.tsv")
health_data <- read_tsv("data/ukbiobank_20231215_health_data_decoded.tsv") 

health_data_met <- 
health_data |> 
  filter(eid %in% surv_all$eid) |> 
  left_join(select(surv_all, eid, sex)) 

table(health_data_met$sex)

health_data_met |> 
  distinct(field_title, eid, .keep_all = TRUE) |> 
  group_by(sex, field_title) |> 
  summarise(n = n()) |> 
  pivot_wider(names_from = sex, values_from = n) |> 
  rename("s0" = `0`,"s1" = `1`) |> 
  View()
  

met_v <- colnames(surv_all)[5:151]
dis_v <- colnames(surv_all)[str_detect(colnames(surv_all), "event")]
dis_v <- sub("^event_", "", dis_v)

health_data |> 
  select(value_type_name) |> 
  distinct() |> 
  View()
health_data |> 
  #filter(value_type_name %in% c("Categorical (single)", "Continuous")) |> 
  group_by(field_title) |> 
  distinct(eid, .keep_all = TRUE) |> 
  summarise(n = n()) |> 
  View()
unique(med$eid)
med <- 
health_data |> 
  filter(field_title == "Medication for cholesterol, blood pressure or diabetes")
med <- 
med |>   
  filter(instance_index == 0) |> 
  group_by(eid) |> 
  mutate(med_n = n()) |> 
  ungroup()
med |> 
  filter(med_n == 3) |> 
  View()
unique(med$value_meaning)
table(med$value_meaning)
table((health_data |> 
         filter(field_title == "Smoking status"))$value_meaning)
chol_t <- 
med |> 
  filter(med_n == 1) |> 
  filter(value_meaning %in% c("None of the above", 
                            "Cholesterol lowering medication")) |> 
  mutate(chol_med = ifelse(value_meaning == "None of the above", "no", "yes")) |> 
  select(eid, chol_med) 

save(chol_t,file =  "data/filters/eids/chol_t.RData")

healthy_w <- 
health_data |>
  filter(instance_index == 0) |> 
  filter(field_title %in% c("Medication for cholesterol, blood pressure or diabetes",
                            "Smoking status",
                            "Overall health rating")) |> 
  select(eid, field_title, value_meaning) |> 
  pivot_wider(names_from = field_title,
              values_from = value_meaning)

# filter for healthy people
health_data |>
  left_join(select(surv_all, eid, sex)) |> 
  filter(instance_index == 0) |> 
  filter(field_title %in% c("Medication for cholesterol, blood pressure or diabetes",
                            "Smoking status",
                            "Overall health rating")) |> 
  select(eid, field_title, value_meaning) |> 
  filter((field_title == "Medication for cholesterol, blood pressure or diabetes" &
  pivot_wider(names_from = field_title,
              values_from = value_meaning) |> 
  filter(`Smoking status` == "Never") |> 
  mutate(have_met = eid %in% surv_all$eid) |> 
  filter(have_met,
         `Overall health rating` %in% 
           c("Excellent", "Fair", "Good")) |> 
  drop_na()


no_drugs_eids <- no_drugs$eid
save(no_drugs_eids,
     file = "data/filters/eids/no_drugs_eids.RData")
# survival with no drugs
