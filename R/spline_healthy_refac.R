library(tidyverse)
library(purrr)
library(furrr)
library(survival)
library(rms)
library(Hmisc)
library(patchwork)
library(cowplot)
library(tictoc)
library(forcats)

surv_all <- read_tsv("data/surv_all_log.tsv") |> 
  mutate(eid = as.character(eid), sex=as.character(sex), spectrometer=as.character(spectrometer))
health_filters <- read_tsv("data/ukb_healthy_subpopul_full_table.tsv") |> 
  mutate_all(as.character)
load("data/rdata/met_dis.RData")
source("R/functions/source_functions.R")

surv_all <- surv_all |> 
  left_join(
    select(health_filters, eid, takes_medication_cholesterol, is_healthy_3, is_healthy_7)
  ) 

# plan: 
# for each metabolite-disease pairs, cohorts, models
# 1. select metabolites and diseases
# 2. have the table
# I have table for metabolite-disease, cohort
#   first create numbers, statistics, after that plots- 
#     statistics list structure: one 
#   when plotting, can refer back to the numbers
#   create incidence with table
#   create Cox model with table
# 3. goals: spline plots, incidence plots, model stats, like lrt, fdr p-value
# 
# functions: 
# filter surv_all
# 


#### create metabolite, disease list ####
mets_filtered <- mets_v[c(2,3:10,14,18,22,26,30:63, 64,66, 70, 72, 76, 78, 82, 84, 88, 90, 94,96,100,102,106,108,112,114,118,120,124,126,130,132,136,138,142,144)]
dis_men <- disease_v[disease_v == "Prostate Cancer"]
dis_women <- disease_v[disease_v %in% c("Breast Cancer", "Ovarian Cancer")]
dis_both_sex <- disease_v[! disease_v %in% dis_men & ! disease_v %in% dis_women & disease_v != "T1 Diabetes"]
dis_all <- c(dis_men, dis_women, dis_both_sex)

#### create cohorts ####
eid_cohorts <- list(
  all = surv_all$eid,
  is_healthy_3 = with(surv_all, eid[is_healthy_3 == "1"]),
  is_healthy_7 = with(surv_all, eid[is_healthy_7 == "1"])
)
# create if with sex-specific diseases
# done
# run analysis
dis_met <- 
  expand_grid(dis_all, mets_filtered) |> 
  rename(disease = dis_all, met = mets_filtered) |> 
  mutate(dm = paste(disease, met, sep = "_"))

plan(multisession, workers = 6)

tic()
test <- 
future_map2(dis_met$disease, dis_met$met, function(disease, metabolite){
  cohort_level <- list()
  for(i in seq_along(eid_cohorts)){
    table <- filter_surv_all(surv_all, metabolite, disease, eid_cohorts[i], ifelse(disease %in% dis_men, "men", ifelse(disease %in% dis_women, "women", "both")))
    cohort_level[c(names(eid_cohorts))[i]] <- cox_stats_from_table(table)
  }
  setNames(list(cohort_level),
           paste(disease, metabolite, sep = "_"))
}
   )
toc()
t <- unlist(test, recursive = F)
cohort <- eid_cohorts[1]
i=2
save(t, file = "data/spline_stats/dis_met_spline_stats_list.RData")

# fdr correct
# create table
extract_values <- function(lst) {
  lst %>%
    map_dfr(~as_tibble(flatten(.), .name_repair = "universal"))
}
your_tibble <- map_dfr(t, ~as_tibble(extract_values(.), .name_repair = "universal"), .id = "Variable")
as_tibble(do.call(rbind, t))[[1]]
str(t)
t[1]
as_tibble(do.call(rbind, head(t, 10))) |> 
  unnest(col = everything()) |> 
  View()
stats_table <- 
enframe(t) |> 
  unnest_longer(value) |> 
  rename(cohort=value_id) |> 
  unnest_longer(value) |> 
  rename(
    dm=name,
    stat_id=value_id,
    stat=value
  )
write_tsv(stats_table, file="data/spline_stats/stats_table_raw.tsv")

# fdr correct
stats_table_fdr <- 
stats_table |> 
  pivot_wider(names_from = stat_id,
              values_from = stat) |> 
  group_by(cohort) |> 
  mutate(across(p,\(x) p.adjust(x, method = "fdr"))) |> 
  ungroup() |> 
  separate(dm, into = c("disease", "metabolite"), sep = "_", remove = F)
  
# create ordering for metabolites and diseases
dis_n_sig <- 
stats_table_fdr |> 
filter(cohort == "is_healthy_7") |> 
  mutate(p_bin = ifelse(p < 0.1, "significant", "non-significant")) |> 
  group_by(disease) |> 
  summarise(n = sum(p_bin == "significant")) |> 
  arrange(n)
met_n_sig <- 
  stats_table_fdr |> 
  filter(cohort == "is_healthy_7") |> 
  mutate(p_bin = ifelse(p < 0.1, "significant", "non-significant")) |> 
  group_by(metabolite) |> 
  summarise(n = sum(p_bin == "significant"))  |> 
  arrange(desc(n))
  
  
stats_table_fdr |> 
  filter(cohort == "is_healthy_7") |> 
  mutate(p_bin = ifelse(p < 0.1, "significant", "non-significant")) |> 
  ggplot(aes(x =factor(metabolite, levels = met_n_sig$metabolite), 
             y = factor(disease, levels = dis_n_sig$disease), 
             fill = p_bin)) +
  geom_tile(color = "black") +
  theme_minimal() +
  scale_fill_manual(values = c("grey", "blue")) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    title = "In healthy_7 cohort significance of LR-test of spline Cox over linear",
    x = "Metabolites",
    y = "Diseases"
  ) +
  guides(fill=guide_legend(title="FDR p < 0.1"))
  
stats_table_fdr |> 
  filter(cohort == "is_healthy_7") |> 
  mutate(p_bin = ifelse(p < 0.1, "significant", "non-significant")) |> 
  group_by(disease) |> 
  summarise(n = sum(p_bin == "significant")) |> 
  ggplot(aes(x=factor(disease, rev(dis_n_sig$disease)), 
             y = n/nrow(met_n_sig))) +
  geom_col() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    title = "Frequency of significant spline metabolite LR tests",
    x = "Diseases",
    y = "Frequency"
  )

stats_table_fdr |> 
  filter(cohort == "is_healthy_7") |> 
  mutate(p_bin = ifelse(p < 0.1, "significant", "non-significant")) |> 
  group_by(metabolite) |> 
  summarise(n = sum(p_bin == "significant")) |> 
  ggplot(aes(x=factor(metabolite, met_n_sig$metabolite), 
             y = n/nrow(met_n_sig))) +
  geom_col() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    title = "Frequency of significant spline disease LR tests",
    x = "Metabolites",
    y = "Frequency"
  ) +
  theme(plot.margin = margin(1,1,1,2.2, "cm"))

plot_table <- 
stats_table_fdr |> 
  filter(cohort == "is_healthy_7") |> 
  mutate(p_bin = ifelse(p < 0.1, "significant", "non-significant")) |> 
  filter(p_bin == "significant") |> 
  select(dm, disease, metabolite, lrt) |> 
  arrange(desc(lrt))

# test
metabolite
disease
cox_df <- table[[1]]
surv_obj <- with(cox_df, Surv(time = time, event = event))
model_spline <- cph(surv_obj ~ spectrometer + rcs(age, 4) + rcs(met_log, 4) , data = cox_df)
  tibble(
    pred=predict(model_spline),
    met=cox_df$met_log
  ) |> 
  ggplot(aes(x=met, y=pred))+
  geom_bin2d()+
  geom_smooth(color="red") 


table <- filter_surv_all(surv_all, metabolite, disease, eid_cohorts[3], ifelse(disease %in% dis_men, "men", ifelse(disease %in% dis_women, "women", "both")))
cox_stats_from_table(table)
cox_plot_from_table(table)


# incidence
cox_df <- table[[1]]
quantile_breaks <- quantile(cox_df$met_log, probs = seq(from = 0, to = 1, length.out = 11))
cox_df$value_cut <- cut(cox_df$met_log, breaks=quantile_breaks)
incidence_table <- 
  with(cox_df, table(value_cut, event)) |> 
  as_tibble() |> 
  pivot_wider(names_from = event, values_from = n) |> 
  mutate(freq=`1`/sum(`0`, `1`))

freqs <- incidence_table$freq

incidence_table |> 
  ggplot(aes(x=value_cut, y=freq))+
  geom_col()
incidence_plot_from_table(table)

plot_table <- 
  stats_table_fdr |> 
  filter(cohort == "is_healthy_7") |> 
  mutate(p_bin = ifelse(p < 0.1, "significant", "non-significant")) |> 
  filter(p_bin == "significant") |> 
  select(dm, disease, metabolite, lrt) |> 
  arrange(desc(lrt))

pdf(file = "plots/spline_incidence4.pdf",
    width = 12,
    height = 5)
walk2(plot_table$disease, plot_table$metabolite,function(disease, metabolite){
  table <- filter_surv_all(surv_all, metabolite, disease, eid_cohorts[3], ifelse(disease %in% dis_men, "men", ifelse(disease %in% dis_women, "women", "both")))
  p1 <- cox_plot_from_table(table)
  p2 <- incidence_plot_from_table(table)
  print(p1 + p2 + plot_annotation(paste(disease, metabolite, sep = "_")))
}
)
dev.off()

cox_df
quantile(cox_df$met_log, probs = seq(from = 0, to = 1, length.out = 11))




cox_df |> 
  arrange(met_log) |> 
  mutate(group = cut(1:nrow(cox_df), 10, labels = FALSE)) |> 
  group_by(group) |> 
  summarise(events = sum(event)) |> 
  ggplot(aes(x=factor(group), y=events))+
  geom_col()

tic()
test <- 
  future_map2(dis_met$disease, dis_met$met, function(disease, metabolite){
    cohort_level <- list()
    for(i in seq_along(eid_cohorts)){
      table <- filter_surv_all(surv_all, metabolite, disease, eid_cohorts[i], ifelse(disease %in% dis_men, "men", ifelse(disease %in% dis_women, "women", "both")))
      cohort_level[c(names(eid_cohorts))[i]] <- cox_stats_from_table(table)
    }
    setNames(list(cohort_level),
             paste(disease, metabolite, sep = "_"))
  }
  )
toc()






map(head(t, 10), rbind)
as_tibble(flatten(t))
as_tibble(t, .name_repair = "universal")
as_tibble(t) |> 
  unnest_longer(col = everything()) |> 
  View()
colnames(as_tibble(t))
as_tibble(t(as_tibble(t))) |> 
  unnest_longer(col = everything()) |> 
  mutate(id = colnames(as_tibble(t)))
#### test variables ####
surv_all
# metabolite: 
metabolite <- "Glucose"
disease <- "Prostate Cancer"
cohort <- eid_cohorts[3]

eids <- list(all=surv_all$eid)
disease_sex <- ifelse(disease %in% dis_women, "women", "both")


# test cox_stats_from_table
test_table <- filter_surv_all(surv_all, metabolite, disease, eids, ifelse(disease %in% dis_men, "men", ifelse(disease %in% dis_women, "women", "both")))
cox_stats_from_table(test_table)


names(list_out)=names(list_table)
rm(list_out)
list_out

cox_stats_from_table(test_table)



plan(multisession, workers = 6)


model_cox <-  function(surv_all, eids, diseases, metabolites, comment){
  
  surv_all <- surv_all |> 
    filter(eid %in% eids)
  
  dis_met <- 
    expand_grid(diseases, metabolites) |> 
    rename(disease = diseases, met = metabolites) |> 
    mutate(dm = paste(disease, met, sep = "_"))
  
  
  cox_df_base <- tibble(
    spec = as.factor(surv_all$spectrometer),
    age = surv_all$age,
    sex = as.factor(surv_all$sex))
  
  table_my <- 
    future_map2(.x = dis_met$met, .y = dis_met$disease, .f = function(met_i, dis_i){
      
      cox_df <- table_for_iteration(surv_all, cox_df_base, met_i, dis_i)
      surv_obj <- with(cox_df, Surv(time = time, event = event))
      model <- cph(surv_obj ~ spec + rcs(age, 4) + rcs(met_cox, 4) + sex, data = cox_df)
      model_base <- cph(surv_obj ~ spec + rcs(age, 4) + sex, data = cox_df)
      model_linear <- cph(surv_obj ~ spec + rcs(age, 4) + met_cox + sex, data = cox_df)
      
      
      pc_tibble <- tibble(
        psl=lrtest(model, model_linear)$stats[3],
        psb=lrtest(model, model_base)$stats[3],
        csl=lrtest(model, model_linear)$stats[1],
        csb=lrtest(model, model_base)$stats[1]
      )
      
      
      tibble(disease=dis_i, metabolite=met_i) |> 
        bind_cols(pc_tibble) 
      
    }
    ) |> 
    list_rbind()
  
  return(table_my)
}