library(tidyverse)
library(foreach)

# read data
met_table <- read_tsv("data/part_met_wide_log.tsv") |> 
  mutate(eid=as.character(eid))
met_table[met_table == 0] <- NA
met_table <- na.omit(met_table)
# mets <- colnames(met_table)
# mets_v <- mets[mets!="eid"]
other_table <- read_tsv("data/other_shorter.tsv") |> 
  mutate(eid=as.character(eid))
disease_table <- read_tsv("data/ukbiobank_diagnosis_table_20230922_with_mortality.txt") |> 
  mutate(eid=as.character(eid))
# disease_v <- unique(disease_table$disease_group)
no_followup_table <- read_tsv("data/ukbiobank_20230926_mortalityCensoring_noFollowup.tsv") |> 
  mutate(eid=as.character(eid))

disease_eids <- unique(disease_table$eid)
no_disease_eids <- setdiff(met_table$eid, disease_eids)
save(no_disease_eids, file = "data/filters/eids/no_disease_eids.RData")
# save(mets_v, disease_v, file="data/rdata/met_dis.RData")

other_table |> 
  distinct(field_title)
# keep first diagnosis
disease_table <- 
disease_table |> 
  group_by(eid, disease_group) |> 
  mutate(first_date=min(diagnosis_date)) |> 
  ungroup() |> 
  filter(diagnosis_date==first_date) |> 
  distinct(eid, disease_group, .keep_all = TRUE) |> 
  select(-first_date)

# filter for eids with all met data

eids <- sort(unique(met_table$eid))
# create survival data with right censoring
eid_meta <- 
tibble(eid=eids) |> 
  full_join(disease_table |> 
              filter(eid %in% eids)) |> 
  left_join(select(no_followup_table, eid, nfo=censoring_date)) |> 
  left_join(select(other_table |> 
                     filter(field_title=="Date of attending assessment centre"), eid, blood_date=value)) |> 
  filter(diagnosis_date>blood_date) |> 
  select(-blood_date)

  
# birth
surv_covariates <- 
tibble(eid=eids) |> 
left_join(select(other_table |> 
                   filter(field_title=="Year of birth"), eid, year=value)) |> 
  left_join(select(other_table |> 
                     filter(field_title=="Month of birth"), eid, month=value)) |> 
  mutate(birth=as.Date(paste0(year, "-", month, "-16"))) |> 
  select(-year, -month)  |> 
  left_join(select(other_table |> 
                     filter(field_title=="Sex"), eid, sex=value)) |> 
  left_join(select(other_table |> 
                     filter(field_title=="Spectrometer"), eid, spectrometer=value)) |> 
  left_join(select(other_table |> 
                     filter(field_title=="Date of attending assessment centre"), eid, blood_date=value)) |> 
  mutate(blood_date=as.Date(blood_date)) |>
  mutate(age=(blood_date-birth)/365.25,
         age=as.numeric(age)) |> 
  select(-birth, -blood_date)

# create survival data table
disease_groups <- sort(unique(disease_table$disease_group))
eid <- "5298556"
event <- eid_meta

eid_disease <- expand.grid(eids, disease_groups) |> 
  rename(eid=Var1, disease_group=Var2) |> 
  left_join(eid_meta) |> 
  left_join(select(other_table |> 
                     filter(field_title=="Date of attending assessment centre"), eid, blood_date=value)) |> 
  mutate(time=ifelse(!is.na(diagnosis_date), diagnosis_date-as.Date(blood_date), ifelse(!is.na(nfo), nfo-as.Date(blood_date), ifelse(disease_group=="all-cause mortality", as.Date("2022-12-19")-as.Date(blood_date), as.Date("2022-10-31")-as.Date(blood_date))))) |> 
  mutate(event=ifelse(is.na(diagnosis_date),0,1))


summary(eid_disease |> 
            select(time) |> 
            mutate(time=time/365.25) |> 
            unlist(use.names = FALSE)
            )
  
surv_data_all <- 
eid_disease |> 
  select(eid, disease_group, time, event) |> 
  mutate(time=time/365.25) |> 
  pivot_wider(names_from = disease_group, values_from = c(event, time))

eid_disease |> 
  mutate()
print(object.size(eid_meta), units = "auto")
# add censoring
disease_table |> 
  left_join(birth_table) |> 
  left_join(no_followup_table) |> 
  mutate(event=ifelse(!is.na(censoring_date),ifelse(censoring_date<diagnosis_date,0,1),1)) |> 
  View()

  mutate(censoring_date=if(!is.na(censoring_date)){
    
  })
disease_table |> 
  filter(eid=="4578564")

disease_table |> 
  left_join(birth_table) |> 
  left_join(no_followup_table) |> 
  distinct(eid)
# mortality table

# mt <- 
# disease_table |> 
#   filter(disease_group=="all-cause mortality") |> 
#   rename(death_date=diagnosis_date) |> 
#   select(eid, death_date)
# 
# disease_table |> 
#   filter(disease_group!="all-cause mortality") |> 
#   left_join(mt) |> 
#   filter(!is.na(death_date)) |> 
#   filter(diagnosis_date>death_date)
# 
# disease_table |> 
#   filter(disease_group=="all-cause mortality") |> 
#   filter(value_meaning=="lost to follow-up, death reported") |> 
#   View()
# 


# survival data all
surv_all <- 
tibble(eid=eids) |> 
  left_join(surv_covariates) |> 
  left_join(met_table) |> 
  left_join(surv_data_all)
write_tsv(surv_all, "data/surv_all_log.tsv")












