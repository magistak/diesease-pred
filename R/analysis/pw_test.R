# read surv_all table with needed variables
source("R/functions_upd/source_functions.R")
source_functions()
load_libraries()


load("data/rdata/met_dis.RData")

met_filtered <- mets_v[c(3:10,14,18,22,26,30:63, 64,66, 70, 72, 76, 78, 82, 84, 88, 90, 94,96,100,102,106,108,112,114,118,120,124,126,130,132,136,138,142,144)] |> make_clean_names()
met_all <- c(met_filtered, "non_hdl")
dis_men <- disease_v[disease_v == "Prostate Cancer"] |> make_clean_names()
dis_women <- disease_v[disease_v %in% c("Breast Cancer", "Ovarian Cancer")] |> make_clean_names()
dis_both_sex <- c(disease_v[! disease_v %in% dis_men & ! disease_v %in% dis_women & disease_v != "T1 Diabetes"] |> make_clean_names(), "cvd")

disease_sex <- tibble(
  disease = c(dis_men, dis_women, dis_both_sex),
  sex = c(rep(1, length(dis_men)), rep(0, length(dis_women)), rep(2, length(dis_both_sex)))
)


# 
# met_list <- c("Albumin","Creatinine","HDL Cholesterol","Clinical LDL Cholesterol","Apolipoprotein A1","Apolipoprotein B", "Total Cholesterol", "Total Triglycerides", "non_HDL") |> make_clean_names()
# clin_list <- c("Albumin_clinical", "Creatinine_clinical", "HDL cholesterol_clinical", "LDL direct_clinical", "Apolipoprotein A_clinical", "Apolipoprotein B_clinical", "Cholesterol_clinical", "Triglycerides_clinical", "non_HDL_clinical") |> make_clean_names()
# biom_list <- c(met_list, clin_list)
# disease_list <- c("all-cause mortality", "T2 Diabetes", "COPD", "Liver Disease", "cvd") |> make_clean_names()

health_filters <- c("is_healthy_1",  "is_healthy_3", "is_healthy_7")


surv_all <- 
  read_surv_all(path = "data/surv_all_bmi_clin.tsv",
                dis_list = disease_sex$disease,
                nmr_list = met_all,
                clin_list = c(),
                health_filter = health_filters) 


eids <- with(surv_all, eid[is_healthy_1 == 1])

dis_met <- 
dis_met_tbl(disease_sex$disease, met_all, c())
# test
dis_met_test <- dis_met |> head(5)
met <- "total_cholesterol"
disease <- "all_cause_mortality"





model_vars <- c("sex", "spectrometer", "age", "met")
model_vars <- model_vars_update(model_vars, disease, disease_sex)

cox_table <- filter_surv_all(surv_all = surv_all,
                         met = met,
                         disease = disease,
                         eids = eids,
                         model_vars = model_vars)

# standardize metabolite level
cox_table <- cox_table |> 
  mutate(met = as.numeric(scale(met)))


minimum_cox <- cox_from_table(cox_table, model_vars, 
                              spline_vars = c("met", "age"))

met_min <- spline_minimum(minimum_cox)
met_min

pw_cox <- pw_from_min(cox_table, met_min, model_vars, "met", c("age"))
lin_cox <- cox_from_table(cox_table, model_vars, "age")

hrs <- hrs_from_pw(pw_cox)
lrt_p <- lrt_stats(pw_cox, lin_cox)

tibble(
  disease = disease,
  met = met,
  hr_l = hrs[1],
  hr_r = hrs[2],
  lrt_p = lrt_p
)







