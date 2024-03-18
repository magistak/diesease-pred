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
  left_join(dis_met) |> 
  head(1)


filters_tbl <- 
  surv_all |> 
  select(eid, starts_with("is_healthy"), sex, bmi) |> 
    mutate(bmi_filter = ifelse(bmi >= 18.5 & bmi <= 24.9, 1, 0)) |> 
    mutate(is_healthy_all = 1) 

plan(multisession, workers = 6)
tic()
out_conf <- future_pmap(select(cox_variations, disease, met, sex_filter, health_filter, bmi_filter), function(disease, metabolite, sex_filter, health_filter, bmi_filter){
  a_s <- if (sex_filter == "both") c("0", "1") else if (sex_filter == "men") 1 else 0
  a_bmi <- if (bmi_filter == "yes") 1 else c(0, 1)
  
eid_filter <- 
filters_tbl |> 
    filter(!!sym(health_filter) == 1) |> 
    filter(sex %in% a_s) |> 
    filter(bmi_filter %in% a_bmi)
eid_filter <- eid_filter$eid
    
hr_from_surv_all_10_80_10(surv_all = surv_all,
                          model_vars = model_vars,
                          eid_filter = eid_filter,
                          metabolite = metabolite,
                          disease = disease) |> 
    mutate(disease = disease,
           metabolite = metabolite,
           edge_name = c("lower", "upper"),
           sex_filter =sex_filter,
           health_filter = health_filter,
           bmi_filter = bmi_filter)
}
)
out_conf <- list_rbind(out_conf)
toc()



heat_tbl <- 
  out_conf |> 
  left_join(distinct(select(dis_met, disease,metabolite= met, nmr_clin)), by = c("disease", "metabolite")) |> 
  mutate(hr_conf = ifelse(lim_2_5 > 1 & lim_97_5 > 1, "above", ifelse(lim_2_5 < 1 & lim_97_5 < 1, "below", "is1"))) |> 
  select(disease:hr_conf, ) |> 
  pivot_wider(names_from = edge_name,
              values_from = hr_conf) |> 
  mutate(heat_group = ifelse(lower == "above" & upper == "above", "above", ifelse(lower == "below" & upper == "below", "below", "other"))) |> 
  arrange(match(disease, disease_list), match(metabolite, biom_list)) 

write_rds(heat_tbl, "apps/app_10_80_10/data/heat_tbl.rds")

out_conf |> 
  left_join(distinct(select(dis_met, disease,metabolite= met, nmr_clin)), by = c("disease", "metabolite")) |> 
  arrange(match(disease, disease_list), match(metabolite, biom_list))

plot_vars <- 
expand_grid(
  sex_filter_var=c("men", "women", "both"),
  health_filter_var=c("is_healthy_1", "is_healthy_3", "is_healthy_7", "is_healthy_all"),
  bmi_filter_var=c("yes", "no")) |> 
  arrange(sex_filter_var, health_filter_var, bmi_filter_var)

pdf(width = 20,
    height = 10,
    file = "plots/bin_10_80_10v2.pdf")
pwalk(plot_vars, function(sex_filter_var, health_filter_var, bmi_filter_var){
p <- 
heat_tbl |> 
    filter(sex_filter == sex_filter_var,
         health_filter == health_filter_var,
         bmi_filter == bmi_filter_var) |> 
  left_join(distinct(select(dis_met, disease,metabolite= met, nmr_clin)), by = c("disease", "metabolite")) |> 
  arrange(match(disease, disease_list), match(metabolite, biom_list)) |> 
  ggplot(aes(x = fct_inorder(metabolite), y = fct_inorder(disease),  fill = heat_group))+
  geom_tile(color = "black") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  facet_wrap(~nmr_clin, scales = "free_x") +
  labs(
    title = sprintf("Sex: %s\nHealth filter: %s\nBmi filter: %s", sex_filter_var, health_filter_var, bmi_filter_var),
    x = "",
    y = ""
  ) +
  scale_fill_manual(values = c("above" = "green", "below" = "red", "other" = "blueviolet")) +
  theme(text = element_text(size = 16))
print(p)
}
)
dev.off()

sex_filter_var = "men"
health_filter_var = "is_healthy_1"
bmi_filter_var = "yes"

heat_tbl |> 
  left_join(distinct(select(dis_met, disease,metabolite= met, nmr_clin)), by = c("disease", "metabolite")) |> 
  filter(heat_group == "above",
         nmr_clin == "nmr") |> 
  group_by(disease, metabolite) |> 
  summarise(n = n()) |> 
  arrange(desc(n)) |> 
  ungroup() |> 
  gt()

# create spline plots
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


plan(multisession, workers = 6)
tic()
future_pwalk(select(cox_variations, disease, met, sex_filter, health_filter, bmi_filter), function(disease, metabolite, sex_filter, health_filter, bmi_filter){
  a_s <- if (sex_filter == "both") c("0", "1") else if (sex_filter == "men") 1 else 0
  a_bmi <- if (bmi_filter == "yes") 1 else c(0, 1)
  
  eid_filter <- 
    filters_tbl |> 
    filter(!!sym(health_filter) == 1) |> 
    filter(sex %in% a_s) |> 
    filter(bmi_filter %in% a_bmi)
  eid_filter <- eid_filter$eid
  
  table <- filter_surv_all(surv_all = surv_all,
                           metabolite = metabolite,
                           disease = disease,
                           eids = eid_filter,
                           model_vars = model_vars)
  model_vars <- if (sex_filter == "men" | sex_filter == "women") model_vars[model_vars!="sex"] else model_vars
  
  cox_model <- 
    cox_from_table(cox_table = table,
                   model_vars = model_vars,
                   spline_vars = "met")
  
  p <- spline_plot_from_cox(cox_model) +
    xlab(metabolite)+
    labs(
      title = paste(disease, metabolite),
      subtitle = sprintf("Sex: %s\nHealth filter: %s\nBmi filter: %s", sex_filter, health_filter, bmi_filter)
    ) 
  
  clin_nmr <- ifelse(metabolite %in% met_list, "nmr", "clin")
  file_name <- paste(disease, metabolite, sex_filter, health_filter, bmi_filter,clin_nmr, sep = "#")
  file_name <- paste0(file_name, ".png")
  ggsave(file_name, p, path = "apps/app_10_80_10/data/spline_pics/",
         width = 8, height = 6, units = "in", dpi = 100)

}
)
out_conf <- list_rbind(out_conf)
toc()
