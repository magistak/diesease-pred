# read surv_all table with needed variables
source("R/functions_upd/source_functions.R")
source_functions()
load_libraries()
library(mgcv)
met_list <- c("Albumin","Creatinine","HDL Cholesterol","Clinical LDL Cholesterol","Apolipoprotein A1","Apolipoprotein B", "Total Cholesterol", "Total Triglycerides", "non_HDL") |> make_clean_names()
clin_list <- c("Albumin_clinical", "Creatinine_clinical", "HDL cholesterol_clinical", "LDL direct_clinical", "Apolipoprotein A_clinical", "Apolipoprotein B_clinical", "Cholesterol_clinical", "Triglycerides_clinical", "non_HDL_clinical") |> make_clean_names()
biom_list <- c(met_list, clin_list)
disease_list <- c("all-cause mortality", "T2 Diabetes", "COPD", "Liver Disease", "cvd") |> make_clean_names()
health_filters <- c("is_healthy_1",  "is_healthy_3", "is_healthy_7")

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

hr_from_met_lh <- function(met_lh){
  met <- met_lh$met[1:2]
  lh <- met_lh$lh[1:2]
  first <- ifelse(met[1]<met[2],1,2)
  sec <- ifelse(first == 1, 2,1)
  coef <- (lh[sec]-lh[first])/(met[sec]-met[first])
  exp(coef)
}


eids <- surv_all$eid
model_vars <- c("sex", "spectrometer", "age", "met")
cox_variations <- 
  expand_grid(
    sex_filter=c("men", "women", "both"),
    health_filter=c("is_healthy_1", "is_healthy_3", "is_healthy_7", "is_healthy_all"),
    disease = unique(dis_met$disease),
    met = unique(dis_met$met)) |> 
  left_join(dis_met) |> 
  select(metabolite = met, 
         disease,
         nmr_clin,
         sex_filter,
         health_filter)

filters_tbl <- 
  surv_all |> 
  select(eid, starts_with("is_healthy"), sex, bmi) |> 
  mutate(bmi_filter = ifelse(bmi >= 18.5 & bmi <= 24.9, 1, 0)) |> 
  mutate(is_healthy_all = 1) 

plan(multisession, workers = 6)
cox_variations_head <- cox_variations |> head(3)
tic()
met_min_tbl <- 
future_pmap(cox_variations, function(metabolite, disease, nmr_clin, sex_filter, health_filter){

  a_s <- if (sex_filter == "both") c("0", "1") else if (sex_filter == "men") 1 else 0

  eid_filter <- 
    filters_tbl |> 
    filter(!!sym(health_filter) == 1) |> 
    filter(sex %in% a_s)
  
  eids <- eid_filter$eid
  model_vars <- c("sex", "spectrometer", "age", "met")
  
if(sex_filter!="both"){
  model_vars <- model_vars[model_vars != "sex"]
}

table <- filter_surv_all(surv_all = surv_all,
                         metabolite = metabolite,
                         disease = disease,
                         eids = eids,
                         model_vars = model_vars)


cox_model <- 
  cox_from_table(cox_table = table,
                 model_vars = model_vars,
                 spline_vars = c("met", "age"))
met_pred <- 
  tibble(
    met = orig_met_from_spline_cox(cox_model),
    pred = predict(cox_model)
  )

met_seq <- seq(min(met_pred$met), max(met_pred$met), length.out = 100)


model <- gam(pred ~ s(met), data = met_pred)

met_seq <- seq(min(met_pred$met), max(met_pred$met), length.out = 100)
pred_seq <- predict(model, newdata = data.frame(met = met_seq))


min_index <- which.min(pred_seq)

if(min_index < 95 & min_index > 5){
met_min <- met_seq[min_index]
} else {
met_min <-  as.numeric(NA)
}

tibble(
  disease = disease,
  metabolite = metabolite,
  bmi_filter = "no",
  health_filter = health_filter,
  sex_filter = sex_filter,
  met_min = met_min,
  nmr_clin = nmr_clin
)
},.options = furrr_options(seed = TRUE), .progress = TRUE
) |> bind_rows()
toc()
met_min_tbl

# save met_min_tbl
write_rds(
  met_min_tbl,
  "data/v_shape_piecewise/met_min_tbl_v2.rds"
)


# Create the V-shaped spline term

# Fit the Cox model using coxph() from the survival package
met_min_tbl <- read_rds("data/v_shape_piecewise/met_min_tbl_v2.rds")
disease <- met_min_tbl$disease[2]
metabolite <- met_min_tbl$metabolite[2]
met_min <- met_min_tbl$met_min[2]
met_min_tbl_head <- met_min_tbl |> 
  drop_na(met_min) |> 
  head(3)

hr_from_met_lh <- function(met_lh){
  met <- met_lh$met[1:2]
  lh <- met_lh$lh[1:2]
  first <- ifelse(met[1]<met[2],1,2)
  sec <- ifelse(first == 1, 2,1)
  coef <- (lh[sec]-lh[first])/(met[sec]-met[first])
  exp(coef)
}
met_min_tbl_head <- met_min_tbl |> 
  drop_na(met_min) |> 
  head(10)
# test
test_list <- as.list(met_min_tbl_head[1,])
disease = test_list[[1]]
metabolite = test_list[[2]]
bmi_filter = test_list[[3]]
health_filter = test_list[[4]]
sex_filter = test_list[[5]]
met_min = test_list[[6]]
nmr_clin = test_list[[7]]
met_min_tbl_v <- met_min_tbl |> drop_na(met_min)

tic()
v_hr_tbl <- 
future_pmap(met_min_tbl_v, function(disease, metabolite, bmi_filter, health_filter, sex_filter, met_min, nmr_clin){
  a_s <- if (sex_filter == "both") c("0", "1") else if (sex_filter == "men") 1 else 0
  
  eid_filter <- 
    filters_tbl |> 
    filter(!!sym(health_filter) == 1) |> 
    filter(sex %in% a_s)
  eids <- eid_filter$eid
  model_vars <- c("sex", "spectrometer", "age", "met")
  
  if(sex_filter!="both"){
    model_vars <- model_vars[model_vars != "sex"]
  }
  
table <- filter_surv_all(surv_all = surv_all,
                         metabolite = metabolite,
                         disease = disease,
                         eids = eids,
                         model_vars = model_vars)

table_v <- 
  table |> 
  mutate(met_left = pmax(met_min - table$met, 0 ),
         met_right = pmax(met - met_min, 0))
table_v$spectrometer

if(sex_filter == "both"){
cox_v_covars <- cph(Surv(time = time, event = event) ~ met_left + met_right + sex + spectrometer + age, data = table_v)
cox_lin_covars <- cph(Surv(time = time, event = event) ~ met + sex + spectrometer + age, data = table_v)
} else {
  cox_v_covars <- cph(Surv(time = time, event = event) ~ met_left + met_right +  age + spectrometer, data = table_v)
  cox_lin_covars <- cph(Surv(time = time, event = event) ~ met + spectrometer + age, data = table_v)
}

cox_v <- cph(Surv(time = time, event = event) ~ met_left + met_right, data = table_v)
cox_lin <- cph(Surv(time = time, event = event) ~ met, data = table_v)


met_lh_left <- 
  tibble(
    met=table$met[table$met < met_min],
    lh = predict(cox_v)[table$met < met_min]
  ) 

met_lh_right <- 
  tibble(
    met=table$met[table$met >= met_min],
    lh = predict(cox_v)[table$met >= met_min]
  ) 


covars_tbl <- 
tibble(
  lh = predict(cox_v_covars),
  met = table_v$met,
) |> 
  arrange(met)

nocovars_tbl <- 
  tibble(
    lh = predict(cox_v),
    met = table_v$met,
  ) |> 
  arrange(met)


model1 <- lm(lh ~ met, data = subset(nocovars_tbl, met <= met_min))
model2 <- lm(lh ~ met, data = subset(nocovars_tbl, met > met_min))

p <- 
covars_tbl |> 
  arrange(met) |> 
  mutate(lms = c(predict(model1), predict(model2))) |> 
  ggplot(aes(x = met, y = lh))+
  geom_bin2d() +
  geom_line(aes(x = met, y= lms),linewidth = 1)+
  xlab("Log biomarker concentration")+
  ylab("Log hazard")+
  labs(
    title = paste(disease, metabolite),
    subtitle = sprintf("Sex: %s\nHealth filter: %s", sex_filter, health_filter)
  ) 
  

file_name <- paste(disease, metabolite, sex_filter, health_filter, bmi_filter,nmr_clin, sep = "#")
file_name <- paste0(file_name, ".png")
ggsave(file_name, p, path = "plots/pics/v_shape_v1/",
       width = 8, height = 6, units = "in", dpi = 100)

hrs <- map_dbl(list(met_lh_left, met_lh_right), hr_from_met_lh)

# lrtest
lrt <- lrtest(cox_v_covars,cox_lin_covars)
lrt_p <- lrt$stats["P"]

tibble(
  disease = disease,
  metabolite = metabolite,
  bmi_filter = bmi_filter,
  health_filter = health_filter,
  sex_filter = sex_filter,
  met_min = met_min,
  nmr_clin = nmr_clin,
  lrt_p = lrt_p,
  hr_left = hrs[1],
  hr_right = hrs[2]
)
}
) |> bind_rows()
toc()

write_rds(
  v_hr_tbl,
  "data/v_shape_piecewise/v_hr_tbl.rds"
)

met_min_tbl <- read_rds("data/v_shape_piecewise/met_min_tbl_v2.rds")
v_hr_tbl <- read_rds("data/v_shape_piecewise/v_hr_tbl.rds")

met_min_tbl |> 
  left_join(v_hr_tbl) |> 
  mutate(p_signif = lrt_p<0.05,
         v_shape = hr_left < 1 & hr_right > 1) |> 
  filter(p_signif & v_shape) |> 
  group_by(health_filter) |> 
  summarise(n = n())





