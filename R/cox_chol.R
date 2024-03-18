library(tidyverse)
library(rms)
library(Hmisc)
library(survival)
library(purrr)
library(furrr)
library(patchwork)

load("data/filters/eids/no_drugs_eids.RData")
load("data/filters/eids/chol_t.RData")
# read table for cox models
surv_all <- read_tsv("data/surv_all_log.tsv")
surv_c <- surv_all |> 
  filter(eid %in% chol_t$eid,
         sex == 1) |> 
  left_join(chol_t)

# surv_c, disease_v, met_list
# example cox
# cox table
cox_df_base <- tibble(
  spec = as.factor(surv_c$spectrometer),
  age = surv_c$age,
  chol = as.factor(surv_c$chol_med))

mets <- c("Glucose", "LDL Cholesterol", "HDL Cholesterol")
diseases <- c("T2 Diabetes", "MACE", "Heart Failure")
dis_met <- expand_grid(mets, diseases) |> 
  mutate(dm = paste(mets, diseases, by = "_"))
i = 1
res_tibble <-   tibble(
  pbl = NA,
  pls = NA,
  cbl = NA,
  cls = NA,
  pred0 = NA,
  pred1 = NA ,
  ap = NA
)
res_tibble <- res_tibble[1,]

pdf("plots/diab_heart_chol.pdf",
    width = 15,
    height = 4)

for(i in 1:nrow(dis_met)){
  met_i  <- dis_met$mets[i]
  dis_i <- dis_met$diseases[i]
  
  cox_df <- 
    cox_df_base %>% 
    mutate(met_cox = surv_c[[met_i]],
           event = surv_c[[paste("event", dis_i, sep = "_")]],
           time = surv_c[[paste("time", dis_i, sep = "_")]])
  
  cox_dfs <- group_split(cox_df, by= chol)
  cox_dfs[[3]] <- bind_rows(cox_dfs[[1]], cox_dfs[[2]])

plots <- list()
  for(j in 1:3) {
 
    
    cox_df <- cox_dfs[[j]]
    title_g <- paste(unique(cox_df$by), collapse = "_")
    
    surv_obj <- with(cox_df, Surv(time = time, event = event))
    
    csp <- cph(surv_obj ~ spec + age + rcs(met_cox, 4), data = cox_df)
    
    plots[[j]] <- tibble(
    pred=predict(csp),
    met=cox_df$met_cox) |> 
    ggplot(aes(x=met, y=pred))+
    geom_bin2d()+
    geom_smooth(color="red") +
    labs(x = met_i,
         y = "Hazard",
         title = paste(met_i, dis_i),
         subtitle = paste("cholesterol med", title_g)
    )
  
  }
print(plots[[1]] + plots[[2]] + plots[[3]])
}
dev.off()

