library(tidyverse)
library(furrr)
library(tictoc)
library(rms)
library(patchwork)
library(broom)

surv_all <- read_tsv("data/surv_all_log.tsv") |> 
  mutate(eid = as.character(eid), sex=as.character(sex), spectrometer=as.character(spectrometer))
health_filters <- read_tsv("data/ukb_healthy_subpopul_full_table.tsv") |> 
  mutate_all(as.character)
clinical_read <- read_tsv("data/raw/biochem_tibble.tsv") |> 
  mutate(eid = as.character(eid)) 
bmi_read <- read_tsv("data/raw/bmi_tibble.tsv") |> 
  mutate(eid = as.character(eid)) 
colnames(clinical_read)[2:ncol(clinical_read)] <- paste(colnames(clinical_read)[2:ncol(clinical_read)], "clinical", sep = "_")
def_clinical <- c("Albumin_clinical", "Creatinine_clinical", "HDL cholesterol_clinical", "LDL direct_clinical", "Apolipoprotein A_clinical", "Apolipoprotein B_clinical", "Cholesterol_clinical", "Triglycerides_clinical", "non_HDL_clinical")


# iqr filter for clinical biomarkers

iqr4_filter <- function(v){
  l <- quantile(v, probs = c(0,0.25,0.5,0.75,1), na.rm = TRUE)[c(2,3,4)]
  v>l[2]-4*(l[2]-l[1]) & v<l[2]+4*(l[3]-l[2])
}

log_add <- function(v1, v2, do="add"){
  switch (do,
    "add" = log10((10^v1 + 10^v2)),
    "sub" = log10((10^v1 - 10^v2)),
    "mult" = log10((10^v1 * 10^v2)),
    "div" = log10((10^v1 / 10^v2))
  )
}

clinical_read <- 
clinical_read |>
  pivot_longer(!eid) |> 
  filter(name %in% def_clinical) |> 
  group_by(name) |> 
  filter(iqr4_filter(value)) |> 
  pivot_wider(names_from = name,
              values_from = value)

load("data/rdata/met_dis.RData")
source("R/functions/source_functions.R")

surv_all <- surv_all |> 
  left_join(
    select(health_filters, eid, takes_medication_cholesterol, is_healthy_3, is_healthy_7, is_healthy_1)
  ) |> 
  left_join(clinical_read) |> 
  left_join(bmi_read)

dis_to_merge <- c("AAA", "Atrial Fibrillation", "Cerebral Stroke", "Heart Failure", "MACE", "PAD", "Venous Thrombosis")


surv_all <- 
  merge_diseases(surv_all, dis_to_merge, "cvd") |> 
  filter(`Cholesterol_clinical`>`HDL cholesterol_clinical` & `Total Cholesterol` > `HDL Cholesterol`) |> 
  mutate(non_HDL_clinical = `Cholesterol_clinical`- `HDL cholesterol_clinical`,
         non_HDL = `Total Cholesterol`- `HDL Cholesterol` ) 

# define metabolites, clinical markers, diseases

def_mets <- c("Albumin","Creatinine","HDL Cholesterol","Clinical LDL Cholesterol","Apolipoprotein A1","Apolipoprotein B", "Total Cholesterol", "Total Triglycerides", "non_HDL")
clin_met <- tibble(def_clinical, def_mets)
dis_men <- disease_v[disease_v == "Prostate Cancer"]
dis_women <- c("Breast Cancer", "Ovarian Cancer")
mets_clin <- tibble(clin=def_clinical,met= def_mets)




colnames(surv_all)[grep("*time*", colnames(surv_all))]
# define diseases:
def_dis <- c("all-cause mortality", "T2 Diabetes", "COPD", "Liver Disease", "cvd")

eid_cohorts <- list(
  all = surv_all$eid,
  is_healthy_7 = with(surv_all, eid[is_healthy_7 == "1"]),
  is_healthy_3 = with(surv_all, eid[is_healthy_3 == "1"]),
  is_healthy_1 = with(surv_all, eid[is_healthy_1 == "1"])
)

dis_met <- 
  expand_grid(def_dis, def_clinical) |> 
  rename(disease = def_dis, met = def_clinical) |> 
  mutate(dm = paste(disease, met, sep = "_"))

plan(multisession, workers = 6)

t <- unlist(test, recursive = F)

#test


surv_all <- 
  drop_na(surv_all)
tic()
lrt_stats <- 
  map2(dis_met$disease, dis_met$met, function(disease, metabolite){
      table <- filter_surv_all(surv_all, metabolite, disease, eid_cohorts["is_healthy_1"], ifelse(disease %in% dis_men, "men", ifelse(disease %in% dis_women, "women", "both")))
      cox_stats_from_table(table)
      stats <- cox_stats_from_table(table)
      t <- rbind(as.numeric(unlist(stats))) |> as_tibble()
      colnames(t) <- c("lrt", "p")
      t |> 
        mutate(disease = disease,
               metabolite = metabolite,
               type = "clinical")

  }
  )
toc()
colnames(lrt_stats)
lrt_stats <- 
  list_rbind(lrt_stats)
lrt_stats_clinical <- 
  lrt_stats |> 
  mutate(type = "clinical") |> 
  rename(clin = metabolite) |> 
  left_join(mets_clin)

dis_met_met <- 
  expand_grid(def_dis, def_mets) |> 
  rename(disease = def_dis, met = def_mets) |> 
  mutate(dm = paste(disease, met, sep = "_"))

tic()
lrt_stats <- 
  map2(dis_met_met$disease, dis_met_met$met, function(disease, metabolite){
    table <- filter_surv_all(surv_all, metabolite, disease, eid_cohorts["is_healthy_1"], ifelse(disease %in% dis_men, "men", ifelse(disease %in% dis_women, "women", "both")))
    cox_stats_from_table(table)
    stats <- cox_stats_from_table(table)
    t <- rbind(as.numeric(unlist(stats))) |> as_tibble()
    colnames(t) <- c("lrt", "p")
    t |> 
      mutate(disease = disease,
             metabolite = metabolite,
             type = "met")
    
  }
  )
toc()
lrt_stats_met <- 
  lrt_stats |> 
  list_rbind() |> 
  rename(met = metabolite) |> 
  left_join(mets_clin)

stats_clin_met <- 
rbind(lrt_stats_clinical, lrt_stats_met) |> 
  group_by(type) |> 
  mutate(p = p.adjust(p, method = "fdr")) |> 
  ungroup() |> 
  pivot_wider(names_from = type,
              values_from = c(lrt, p)) |> 
  head(2)

with(stats_clin_met, cor.test(lrt_clinical, lrt_met, method = "spearman"))
stats_clin_met |> 
  ggplot(aes(x = lrt_clinical, y = lrt_met))+
  geom_point()

plot_table <- 
  stats_clin_met |> 
  arrange(desc(lrt_clinical))

plan(multisession, workers = 6)

# plot tables:
metabolite_met <- plot_table$met[9]
metabolite_clin <- plot_table$clin[9]
disease <- plot_table$disease[9]

plots_bmi <- 
  future_pmap(list(plot_table$disease, plot_table$met, plot_table$clin),function(disease, metabolite_met, metabolite_clin){
    metabolite = metabolite_met
    table <- filter_surv_all(surv_all, metabolite, disease, eid_cohorts[3], ifelse(disease %in% dis_men, "men", ifelse(disease %in% dis_women, "women", "both")))
    p1 <- cox_plot_from_table(table) + theme(legend.position="none") + labs(title = paste(disease, metabolite,sep =  "_"), subtitle = paste("FDR_p:", signif(with(plot_table, p_met[disease==disease & met==metabolite]),4)))
    metabolite = metabolite_clin
    table <- filter_surv_all(surv_all, metabolite, disease, eid_cohorts[3], ifelse(disease %in% dis_men, "men", ifelse(disease %in% dis_women, "women", "both")))
    p2 <- cox_plot_from_table(table) + theme(legend.position="none") + labs(title = paste(disease, metabolite,sep =  "_"), subtitle = paste("FDR_p:", signif(with(plot_table, p_clinical[disease==disease & clin==metabolite]),4)))

    p <- p1 + p2
    p
  }
  )

plots[[2]]
p <- wrap_plots(plots)
pdf(file = "plots/clin_met_op.pdf",
    width = 12,
    height = 5)
for(i in seq_along(plots)){
  print(plots[[i]])
}
dev.off()

stats_table <- 
enframe(test) |> 
  unnest_longer(value) |> 
  rename(cohort=value_id) |> 
  unnest_longer(value) |> 
  rename(
    dm=name,
    cohort=value_id,
    stat=value,
    dm_real=cohort
  ) |> 
  unnest_longer(stat) |> 
  mutate(dm = dm_real) |> 
  select(!dm_real)

stats_table_fdr <- 
  stats_table |> 
  pivot_wider(names_from = stat_id,
              values_from = stat) |> 
  group_by(cohort) |> 
  mutate(across(p,\(x) p.adjust(x, method = "fdr"))) |> 
  ungroup() |> 
  separate(dm, into = c("disease", "metabolite"), sep = "_", remove = F, extra = "merge")

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

plot_table <- 
  stats_table_fdr |> 
  filter(cohort == "is_healthy_7") |> 
  mutate(p_bin = ifelse(p < 0.1, "significant", "non-significant")) |> 
  select(dm, disease, metabolite, lrt) |> 
  arrange(desc(lrt))

pdf(file = "plots/spline_incidence4.pdf",
    width = 12,
    height = 5)
# test
disease <- plot_table$disease[1]
metabolite <- plot_table$metabolite[1]
plots <- 
future_map2(plot_table$disease, plot_table$metabolite,function(disease, metabolite){
  table <- filter_surv_all(surv_all, metabolite, disease, eid_cohorts[3], ifelse(disease %in% dis_men, "men", ifelse(disease %in% dis_women, "women", "both")))
  p1 <- cox_plot_from_table(table) + theme(legend.position="none")
  p2 <- incidence_plot_from_table(table) + theme(legend.position="none")
  p <- p1 + p2 + plot_annotation(paste(disease, metabolite, sep = "_"))
  p
}
)
dev.off()
png(width = 2400, height = 1200, filename = "plots/all_clinical.png", units = "px")
p <- wrap_plots(plots)
pp <- plots[1:2]
names(pp) <- c("aa","bb")
wrap_plots(plots[1:2], tag_level = 'keep')
print(p)
dev.off()
plots[1]

# non_hdl in different cohorts for all-cause mortality
i = 1
non_hdl <- c("non_HDL_clinical", "non_HDL")
disease="all-cause mortality"
ec=1

metabolite = non_hdl[1]
plots_clin <- list()
for(ec in seq_along(eid_cohorts)){
table <- filter_surv_all(surv_all, metabolite, disease, eid_cohorts[ec], ifelse(disease %in% dis_men, "men", ifelse(disease %in% dis_women, "women", "both")))
cox_df <- table[[1]]
cox_df <- 
  cox_df |> 
  mutate(met_quant = quantize_continous(met_log, 5))
surv_obj <- with(cox_df, Surv(time = time, event = event))
model_spline <- cph(surv_obj ~ spectrometer + age + rcs(met_log, 4) + sex , data = cox_df)
model_lin <- cph(surv_obj ~ spectrometer + age + met_log + sex , data = cox_df)
model_quant <- cph(surv_obj ~ met_quant +spectrometer + age +  sex , data = cox_df)
split_vector <- seq(0, 1, length.out = 5 + 1)
breaks <- quantile(cox_df$met_log,probs=split_vector)
hrs <- 
  exp(model_quant$coefficients)[1:4] |> 
  unname() |> 
  signif(3)
lrt_p <- lrtest(model_spline, model_lin)$stats[[3]]
p <- 
  tibble(
    pred=predict(model_spline),
    met=cox_df$met_log
  ) |> 
  ggplot(aes(x=met, y=pred))+
  geom_bin2d()+
  geom_smooth(color="red", size = 1.5) +
  geom_vline(xintercept = breaks, color = "red", size=1.5)+
  labs(title = paste(disease,"_", metabolite, "\n", "p: ", signif(lrt_p,4),sep =  ""), 
       subtitle = paste("cohort:", names(eid_cohorts)[ec], "#people:", length(eid_cohorts[[ec]]), "\n","HR for quantiles:", hrs[1], hrs[2], hrs[3], hrs[4]))
plots_clin[[ec]] <- p
}
plots_clin[[1]]

interleave <- function(list1, list2){
  list_all <- list()
  for(i in seq_along(list1)){
    list_all[[2*i-1]] <- list1[[i]] 
    list_all[[2*i]] <- list2[[i]]
  }
  return(list_all)
}
combined_list <- interleave(plots_clin, plots_nmr)
wrap_plots(combined_list, ncol = 2)
v=cox_df$met_log
n=5
quantize_continous <- function(v, n){
  split_vector <- seq(0, 1, length.out = n + 1)
  breaks <- quantile(v,probs=split_vector)
  cuts <- cut(v, 
              breaks = breaks,
              include.lowest = TRUE,
              right = FALSE)
  as.character(cuts)
}
# quantize continous variable
table <- filter_surv_all(surv_all, metabolite, disease, eid_cohorts[ec], ifelse(disease %in% dis_men, "men", ifelse(disease %in% dis_women, "women", "both")))
cox_df <- table[[1]]
cox_df <- 
cox_df |> 
  mutate(met_quant = quantize_continous(met_log, 5))

surv_obj <- with(cox_df, Surv(time = time, event = event))
model_spline <- cph(surv_obj ~ spectrometer + age + rcs(met_log, 4) + sex , data = cox_df)
model_lin <- cph(surv_obj ~ spectrometer + age + met_log + sex , data = cox_df)
model_quant <- cph(surv_obj ~ met_quant +spectrometer + age +  sex , data = cox_df)
hrs <- 
exp(model_quant$coefficients)[1:4] |> 
  unname()

lrt_p <- lrtest(model_spline, model_lin)$stats[[3]]

pdf(
  height = 10,
  width = 20,
  file = "plots/all_cause_nhdl.pdf"
)
for(i in seq_along(plots_clin)){
  p <- plots_clin[[i]] + plots_nmr[[i]]
  print(p)
}
dev.off()

