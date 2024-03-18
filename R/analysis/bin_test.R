# read surv_all table with needed variables
source("R/functions_upd/source_functions.R")
source_functions()
load_libraries()
library(forcats)

# source("other_functions")

met_list <- c("Albumin","Creatinine","HDL Cholesterol","Clinical LDL Cholesterol","Apolipoprotein A1","Apolipoprotein B", "Total Cholesterol", "Total Triglycerides", "non_HDL") |> make_clean_names()
clin_list <- c("Albumin_clinical", "Creatinine_clinical", "HDL cholesterol_clinical", "LDL direct_clinical", "Apolipoprotein A_clinical", "Apolipoprotein B_clinical", "Cholesterol_clinical", "Triglycerides_clinical", "non_HDL_clinical") |> make_clean_names()
biom_list <- c(met_list, clin_list)
disease_list <- c("all-cause mortality", "T2 Diabetes", "COPD", "Liver Disease", "cvd") |> make_clean_names()
health_filters <- c("is_healthy_1",  "is_healthy_3", "is_healthy_7")

# run with: both sex, health_filter_1, no bmi filter
# read survival data
surv_all <- 
  read_surv_all(path = "data/surv_all_bmi_clin.tsv",
                dis_list = disease_list,
                nmr_list = met_list,
                clin_list = clin_list,
                health_filter = health_filters) 

eid_filter <- with(surv_all, eid[is_healthy_1 == 1])

model_vars <- c("sex", "spectrometer", "age", "met")

dis_met <- dis_met_tbl(disease_list, met_list, clin_list)

cox_variations <- 
  expand_grid(
    sex_filter=c("men", "women", "both"),
    health_filter=c("is_healthy_1", "is_healthy_3", "is_healthy_7", "is_healthy_all"),
    bmi_filter=c("yes", "no"),
    disease = unique(dis_met$disease),
    met = unique(dis_met$met)) |> 
  left_join(dis_met)


filters_tbl <- 
  surv_all |> 
  select(eid, starts_with("is_healthy"), sex, bmi) |> 
  mutate(bmi_filter = ifelse(bmi >= 18.5 & bmi <= 24.9, 1, 0)) |> 
  mutate(is_healthy_all = 1) 

# test
metabolite <- "albumin"
disease <- "all_cause_mortality"
eids <- surv_all$eid
model_vars <- c("sex", "spectrometer", "age", "met")


table <- filter_surv_all(surv_all = surv_all,
                         metabolite = metabolite,
                         disease = disease,
                         eids = eids,
                         model_vars = model_vars)


# cox from table

cox_model <- 
cox_from_table(cox_table = table,
               model_vars = model_vars,
               spline_vars = "met")
met_pred <- 
tibble(
met = orig_met_from_spline_cox(cox_model ),
pred = predict(cox_model)
)
smooth_fit <- loess(pred ~ met, data = met_pred)

p <- 
tibble(
  pred=predict(cox_model),
  met=og_met
) |> 
  ggplot(aes(x=met, y=pred))+
  geom_bin2d()+
  geom_smooth(color="red", size = 1.5) 
print(p)
