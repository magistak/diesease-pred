# read surv_all table with needed variables
source("R/functions_upd/source_functions.R")
source_functions()
load_libraries()


# source("other_functions")

met_list <- c("Albumin","Creatinine","HDL Cholesterol","Clinical LDL Cholesterol","Apolipoprotein A1","Apolipoprotein B", "Total Cholesterol", "Total Triglycerides", "non_HDL") |> make_clean_names()
clin_list <- c("Albumin_clinical", "Creatinine_clinical", "HDL cholesterol_clinical", "LDL direct_clinical", "Apolipoprotein A_clinical", "Apolipoprotein B_clinical", "Cholesterol_clinical", "Triglycerides_clinical", "non_HDL_clinical") |> make_clean_names()
biom_list <- c(met_list, clin_list)
disease_list <- c("all-cause mortality", "T2 Diabetes", "COPD", "Liver Disease", "cvd") |> make_clean_names()
health_filters <- c("is_healthy_1",  "is_healthy_3", "is_healthy_7")


surv_all <- 
  read_surv_all(path = "data/surv_all_bmi_clin.tsv",
                dis_list = disease_list,
                nmr_list = met_list,
                clin_list = clin_list,
                health_filter = health_filters) 

metabolite <- met_list[2]
disease <- disease_list[1]
eids <- with(surv_all, eid[is_healthy_1 == 1])

disease_sex <- tibble(
  disease = disease_list,
  sex = 2
)

model_vars <- c("sex", "spectrometer", "age", "met")
model_vars <- model_vars_update(model_vars, disease, disease_sex)


cox_table <- filter_surv_all(surv_all = surv_all,
                         metabolite = metabolite,
                         disease = disease,
                         eids = eids,
                         model_vars = model_vars)

minimum_cox <- cox_from_table(cox_table, model_vars, 
                              spline_vars = c("met", "age"))

met_min <- spline_minimum(minimum_cox)

pw_from_min(cox_table, met_min, model_vars, "met", c("age"))


