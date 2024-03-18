library(tidyverse)
library(broom)
library(rlang)
library(purrr)
library(patchwork)

surv_all <- read_tsv("data/surv_all_log.tsv") |> 
  mutate(eid = as.character(eid), sex=as.character(sex), spectrometer=as.character(spectrometer))
health_filters <- read_tsv("data/ukb_healthy_subpopul_full_table.tsv") |> 
  mutate_all(as.character)
clinical_read <- read_tsv("data/raw/biochem_tibble.tsv") |> 
  mutate(eid = as.character(eid)) 
colnames(clinical_read)[2:ncol(clinical_read)] <- paste(colnames(clinical_read)[2:ncol(clinical_read)], "clinical", sep = "_")
load("data/rdata/met_dis.RData")
source("R/functions/source_functions.R")

surv_all <- surv_all |> 
  left_join(
    select(health_filters, eid, takes_medication_cholesterol, is_healthy_3, is_healthy_7)
  ) |> 
  left_join(clinical_read)

# define metabolites, clinical markers, diseases
# end-point: all-cause mortality, cardiovascular diseases
# create cardiovascular disease endpoint

dis_to_merge <- c("AAA", "Atrial Fibrillation", "Cerebral Stroke", "Heart Failure", "MACE", "PAD", "Venous Thrombosis")

surv_all <- 
merge_diseases(surv_all, dis_to_merge, "cvd") |> 
  mutate(non_HDL_clinical = `Cholesterol_clinical` - `HDL cholesterol_clinical`,
         non_HDL = `Total Cholesterol` - `HDL Cholesterol`)

spearman_stat <- (with(surv_all, cor.test(non_HDL_clinical, non_HDL, method = "spearman")) |> 
  broom::tidy())[[1]]

surv_all |> 
  ggplot(aes(x= non_HDL_clinical, y = non_HDL)) +
  geom_bin_2d() +
  geom_smooth() +
  labs(subtitle = paste("Spearman rho:", signif(spearman_stat, 3)))

#### plot correlations between nmr and clinical biomarkers ####
def_clinical <- c("Albumin_clinical", "Creatinine_clinical", "HDL cholesterol_clinical", "LDL direct_clinical", "Apolipoprotein A_clinical", "Apolipoprotein B_clinical", "Cholesterol_clinical", "Triglycerides_clinical", "non_HDL_clinical")
def_mets <- c("Albumin","Creatinine","HDL Cholesterol","Clinical LDL Cholesterol","Apolipoprotein A1","Apolipoprotein B", "Total Cholesterol", "Total Triglycerides", "non_HDL")
clin_met <- tibble(def_clinical, def_mets)
clin_mark <- def_clinical[9]
met_mark <- def_mets[9]

p <- plot_cor(surv_all, clin_mark, clin_met)

plots <- 
  map2(clin_met$def_clinical, clin_met$def_mets, function(clin_mark, met_mark){
    spearman_stat <- (cor.test(surv_all[[clin_mark]], surv_all[[met_mark]], method = "spearman") |> 
                        broom::tidy())[[1]]
    
    surv_all |> 
      rename(clin = !!clin_mark, met = !!met_mark) |> 
      ggplot(aes(x= clin, y = met)) +
      geom_bin_2d() +
      geom_smooth() +
      labs(title = paste(clin_mark, met_mark, sep = "_"),
          subtitle = paste("Spearman rho:", signif(spearman_stat, 3)))
    
  })
p <- wrap_plots(plots)






