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

write_rds(filters_tbl, "apps/app_v_v1/data/filters_tbl.rds")

# test
metabolite <- met_list[2]
disease <- disease_list[1]
eids <- surv_all$eid
model_vars <- c("sex", "spectrometer", "age", "met")


table <- filter_surv_all(surv_all = surv_all,
                         metabolite = metabolite,
                         disease = disease,
                         eids = eids,
                         model_vars = model_vars)

# piecewise cox


# cox from table
# calculate minimum


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

tibble(met_seq, pred_seq) |> 
  ggplot(aes(x = met_seq, y = pred_seq)) +
  geom_point()


min_index <- which.min(pred_seq)
met_min <- met_seq[min_index]

# Create the V-shaped spline term

# Fit the Cox model using coxph() from the survival package


table_v <- 
  table |> 
  mutate(met_left = pmax(met_min - table$met, 0 ),
         met_right = pmax(met - met_min, 0))

table_v2 <- 
  table |> 
  mutate(met_x = met,
         met_y = pmax(met - met_min, 0))


cox_v <- cph(Surv(time = time, event = event) ~ met_left + met_right, data = table_v)
cox_v2 <- cph(Surv(time = time, event = event) ~ met_x + met_y, data = table_v2)
coefs <- coefficients(cox_v2)
exp(coefs[1])
exp(coefs[1]+coefs[2])

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

hr_from_met_lh <- function(met_lh){
  met <- met_lh$met[1:2]
  lh <- met_lh$lh[1:2]
  first <- ifelse(met[1]<met[2],1,2)
  sec <- ifelse(first == 1, 2,1)
  coef <- (lh[sec]-lh[first])/(met[sec]-met[first])
  exp(coef)
}

tibble(
  lh = predict(cox_v),
  met = table_v$met,
) |> 
  ggplot(aes(x = met, y = lh))+
  geom_bin2d()

tibble(
  lh = predict(cox_v2),
  met = table_v$met,
) |> 
  ggplot(aes(x = met, y = lh))+
  geom_bin2d()

identical(predict(cox_v), predict(cox_v2))
predict(cox_v) == predict(cox_v2)
hrs <- map_dbl(list(met_lh_left, met_lh_right), hr_from_met_lh)
predict(cox_v) |> tail(5) == predict(cox_v2) |> tail(5)

# lrtest
lrtest(cox_v,cox_lin)

cox_test

tibble(
  pred = predict(cox_v),
  met = table_v$met,
) |> 
  mutate(index = row_number()) |> 
  arrange(met) |> 
  filter(met > met_min)




met_pred_v$met[met_pred_v$met == met_min]

coef <- (met_pred_v$pred[12]-met_pred_v$pred[2])/(met_pred_v$met[12]-met_pred_v$met[2]) 
coef |> exp()

cox_v$coefficients |> exp()

confint(cox_v) |> exp()
predict(cox_test, type = "lp")

tibble(
log_hazard = predict(cox_test, type = "lp"),
met = table_v$met
) |> 
  ggplot(aes(x = met, y = log_hazard))+
  geom_bin2d()
log_hazard = predict(cox_test, type = "lp")
met = table_v$met
coef <- (log_hazard[2]-log_hazard[1])/(met[2]-met[1]) 
coef |> exp()
coefficients(cox_test) |> exp()


table

table_l <- 
  table |> 
  filter(met < met_min)

table_r <- 
  table |> 
  filter(met >= met_min)

cox_model_l <- 
  cox_from_table(cox_table = table_l,
                 model_vars = model_vars,
                 spline_vars = c("age"))

str(cox_model_l)

cox_model_l$scale.pred
coef_values <- cox_model_l$coefficients
vcov_matrix <- vcov(cox_model_l)

# Calculate the standard errors
std_errors <- sqrt(diag(vcov_matrix))

# Calculate the z-values
z_values <- coef_values / std_errors

# Calculate the p-values using the standard normal distribution
p_values <- 2 * pnorm(-abs(z_values))

# Print the p-values
p_values["met"]

lrtest()

confint(cox_model_l) |> 
  as_tibble(rownames = NA) |> 
  rownames_to_column() |> 
  rename(lim_2_5 = `2.5 %`, lim_97_5 = `97.5 %`) |> 
  mutate(hr = (lim_2_5 + lim_97_5)/2) |> 
  mutate(across(
    c(lim_2_5, lim_97_5, hr),
    exp
  )) |> 
  filter(startsWith(rowname, "met"))

met_pred |> 
  ggplot(aes(x = met, y = pred))+
  geom_bin2d() +
  geom_smooth()

# v-shape term



# Fit the Cox model using cph()
cox_model <- cph(Surv(time, status) ~ met_spline, data = data)


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



# segmented
