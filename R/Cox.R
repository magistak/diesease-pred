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
surv_all <- surv_all |> 
  filter(eid %in% chol_t$eid,
         sex == 1) |> 
  left_join(chol_t)
  
# create surv_all_c
met_list <- read_csv("data/biom_meta.csv") |> 
  filter(in_strict == 1) |> 
  select(met_name) |> 
  unlist(use.names = FALSE)
met_list <- colnames(surv_all)[5:151]
load("data/rdata/met_dis.RData")
# mets no need
nn_mets <- setdiff(mets_v, met_list)
csa <- colnames(surv_all)
surv_c_cols <- csa[!csa %in% nn_mets]
surv_c <- surv_all[ , surv_c_cols]

# surv_c, disease_v, met_list
# example cox
# cox table
surv_c <- surv_all
cox_df_base <- tibble(
                 spec = as.factor(surv_c$spectrometer),
                 age = surv_c$age,
                 chol = as.factor(surv_c$chol_med))

table(cox_df_base$sex)
disease_v
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
pdf("plots/diab_heart.pdf")
for(i in 1:nrow(dis_met)){
  met_i  <- dis_met$mets[i]
  dis_i <- dis_met$diseases[i]
  
cox_df <- 
  cox_df_base %>% 
  mutate(met_cox = surv_c[[met_i]],
         event = surv_c[[paste("event", dis_i, sep = "_")]],
         time = surv_c[[paste("time", dis_i, sep = "_")]])

cox_dfs <- group_split(cox_df, by= chol)

cox_df <- cox_dfs[[2]]

surv_obj <- with(cox_dfs[[2]], Surv(time = time, event = event))

cb <- cph(surv_obj ~ spec + age, data = cox_df)
clin <- cph(surv_obj ~ spec + age + met_cox, data = cox_df)
csp <- cph(surv_obj ~ spec + age + rcs(met_cox, 4), data = cox_dfs[[2]])

lt <- lrtest(clin, cb)
st <- lrtest(csp, clin)
pred00 <- Predict(csp, sex = 0, age = mean(cox_df$age), spec = 3, met_cox = seq(min(cox_df$met_cox), max(cox_df$met_cox), length = 200))
pred11 <- Predict(csp, sex = 1, age = mean(cox_df$age), spec = 3, met_cox = seq(min(cox_df$met_cox), max(cox_df$met_cox), length = 200))

tib <- tibble(
  pbl = lt$stats[3] %>% unname(),
  pls = st$stats[3] %>% unname(),
  cbl = lt$stats[1] %>% unname(),
  cls = st$stats[1] %>% unname()
)

all_plot <- tibble(
  pred=predict(csp),
  met=cox_df$met_cox) |> 
  ggplot(aes(x=met, y=pred))+
  geom_bin2d()+
  geom_smooth(color="red") +
  labs(x = met_i,
       y = "Hazard",
       title = paste(met_i, dis_i),
       subtitle = paste(colnames(tib), signif(as.numeric(tib[1,]), 4), collapse = " ")
       )
  
print(all_plot)
}
dev.off()
cox_df <- 

dis_met <- 
  expand_grid(disease_v, met_list) |> 
  rename(disease = disease_v, met = met_list) |> 
  mutate(dm = paste(disease, met, sep = "_")) |> 
  head(10)

# map
surv_all <- surv_all |> 
  mutate(sex=as.factor(sex),
         spectrometer=as.factor(spectrometer)) |> 
  drop_na()

# analyze results
# fdr
dis_met2 <- 
dis_met |> 
  mutate(pbl=p.adjust(pbl, method="fdr"),
         pls=p.adjust(pls, method="fdr")) |> 
  mutate(sig_s = pls<0.05,
         sig_l = pbl<0.05)

mean(dis_met2$sig_s, na.rm=T)
dis_met2 |> 
  group_by(disease) |> 
  summarise(n=n(), nsig=sum(sig_s)) |> 
  ungroup() |> 
  mutate(percent=nsig/n*100) |> 
  View()

dis_met2 |> 
  group_by(met) |> 
  summarise(n=n(), nsig=sum(sig_s)) |> 
  ungroup() |> 
  mutate(percent=nsig/n*100) |> 
  View()

dis_met2 |> 
  ggplot(aes(x=pls))+
  geom_histogram()
# number 

# statistics for significant spline


write_tsv(dis_met, "data/dis_met1.tsv")
# working with pairs



plan(multisession, workers = 7)

res <- dis_met %>%
  mutate(
    result = future_pmap(list(met, disease), ~ {
      cox_df <- 
        cox_df_base %>% 
        mutate(met_cox = surv_c[[.x]],
               event = surv_c[[paste("event", .y, sep = "_")]],
               time = surv_c[[paste("time", .y, sep = "_")]])

      surv_obj <- with(cox_df, Surv(time = time, event = event))
      
      cb <- cph(surv_obj ~ sex + spec + age, data = cox_df)
      clin <- cph(surv_obj ~ sex + spec + age + met_cox, data = cox_df)
      csp <- cph(surv_obj ~ sex + spec + age + rcs(met_cox, 4), data = cox_df)
      
      lt <- lrtest(clin, cb)
      st <- lrtest(csp, clin)
      pred00 <- Predict(csp, sex = 0, age = mean(cox_df$age), spec = 3, met_cox = seq(min(cox_df$met_cox), max(cox_df$met_cox), length = 200))
      pred11 <- Predict(csp, sex = 1, age = mean(cox_df$age), spec = 3, met_cox = seq(min(cox_df$met_cox), max(cox_df$met_cox), length = 200))
      
      all_plot <- tibble(
        pred=predict(csp),
        met=cox_df$met_cox) |> 
        ggplot(aes(x=met, y=pred))+
        geom_bin2d()+
        geom_smooth(color="red")
      
      tibble(
        pbl = lt$stats[3] %>% unname(),
        pls = st$stats[3] %>% unname(),
        cbl = lt$stats[1] %>% unname(),
        cls = st$stats[1] %>% unname(),
        pred0 = list(pred00),
        pred1 = list(pred11),
        ap <- list(all_plot)
      )
    })
  ) 
res
sex_dis <- c("Prostate Cancer", "Breast Cancer", "Ovarian Cancer")
resw <- res|> 
  filter(!disease %in% sex_dis) |> 
  select(disease, met, dm, result) |> 
  unnest_wider(result)
save(resw,file =  "data/resw.RData")
resw |> 
  mutate(pls = p.adjust(pls, method = "fdr"),
         pbl = p.adjust(pbl, method = "fdr"),
         sigls = pls < 0.05,
         sigbl = pbl < 0.05
         ) |>
  group_by(met) |> 
  mutate(m_sigls = sum(sigls),
         m_all = n()) |> 
  View()
i=1
i
pdf(width=10, height = 6,
    file="plots/termplots2.pdf")
resw <- resw |> 
  arrange(desc(cls))
for(i in 1:nrow(resw)){
pred0 <- resw$pred0[[i]][[1]]
pred1 <- resw$pred1[[i]][[1]]
p <- ggplot(pred0)+
  xlab("")+
  ggtitle(resw$met[i]) +
  labs(subtitle=paste0("Linear lrt: ",signif(resw$cbl[i],3)))+
  ggplot(pred1)+
  xlab("") +
  ggtitle(resw$disease[i]) +
  labs(subtitle=paste0("Spline lrt: ",signif(resw$cls[i],3)),
       caption = paste0("P(fdr): ", signif(resw$pls[i],3)))

print(p)
}
dev.off()
# heatplot

resw |> 
  select(-c("pred0", "pred1")) |> 
  mutate(pls = p.adjust(pls, method = "fdr"),
         pbl = p.adjust(pbl, method = "fdr"),
         sigls = pls < 0.05,
         sigbl = pbl < 0.05
  ) |> 
  ggplot(aes(met, reorder(disease, as.numeric(sigls)), fill = sigls)) +
  geom_tile() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
resw <- 
resw |> 
  mutate(slr=cls/cbl) |> 
  arrange(desc(slr))
  
resw$pred0[1]

res |> 
  unnest_wider(result)

res$result[1]
#%>%
  unnest_wider(result) %>%
  select(-result)
dis_met <- 
    expand_grid(disease_v, met_list) |> 
    rename(disease = disease_v, met = met_list) |> 
    mutate(dm = paste(disease, met, sep = "_"))
i=1
pdf(file = "plots/term_sss78.pdf")
for(i in 601:720){
met <- dis_met$met[i]
dis <- dis_met$disease[i] 

cox_df <- 
    cox_df_base %>% 
    mutate(met_cox = surv_c[[met]],
           event = surv_c[[paste("event", dis, sep = "_")]],
           time = surv_c[[paste("time", dis, sep = "_")]])
  

 
#cb <- cph(surv_obj ~ sex + spec + age, data = cox_df)
#clin <- cph(surv_obj ~ sex + spec + age + met_cox, data = cox_df)

cox_df_0 <- cox_df |> 
  filter(sex == 0)
surv_obj0 <- with(cox_df_0, Surv(time = time, event = event))
cox_df_1 <- cox_df |> 
  filter(sex == 1)
surv_obj1 <- with(cox_df_1, Surv(time = time, event = event))
dd <- datadist(cox_df_0)
options(datadist="dd")
csp0 <- cph(surv_obj0 ~ spec + rcs(age, 4) + rcs(met_cox, 4), data = cox_df_0)
dd <- datadist(cox_df_1)
options(datadist="dd")
csp1 <- cph(surv_obj1 ~ spec +rcs(age, 4) + rcs(met_cox, 4), data = cox_df_1)
pred0 <- Predict(csp0)
pred1 <- Predict(csp1)
p0 <- ggplot(pred0)
p1 <- ggplot(pred1)
p <- patchwork::wrap_plots(p0, p1) + 
  plot_annotation(
  title = paste(met, dis)
  )

print(p)
}
surv_c |> 
  group_by(sex) |> 
  summarise(prostate = sum(`event_Rectal Cancer`))

dev.off()
csp
data.class(csp)
Predict(csp) 
?predict.cph
all_plot <- tibble(
pred=predict(csp),
met=cox_df$met_cox
)
p_all <- 
all_plot |> 
  ggplot(aes(x=met, y=pred))+
  geom_bin2d()+
  geom_smooth(color="red")
pred00 <- Predict(csp, sex = 0, age = mean(cox_df$age), spec = 3, met_cox = seq(min(cox_df$met_cox), max(cox_df$met_cox), length = 200))
ggplot(pred00) + p_all

plot(pr)  
pred0 <- Predict(csp, sex = 0, age = mean(cox_df$age), spec = 3, met_cox = seq(min(cox_df$met_cox), max(cox_df$met_cox), length = 200))
  metterm<- pred0$met_cox
  center<-  
   
pred1 <- Predict(csp, sex = 1, age = mean(cox_df$age), spec = 3, met_cox = seq(min(cox_df$met_cox), max(cox_df$met_cox), length = 200))
  
pred0 <- Predict(csp, sex = 0, age = mean(cox_df$age), spec = 3, met_cox = seq(min(cox_df$met_cox), max(cox_df$met_cox), length = 200))
pred1 <- Predict(csp, sex = 1, age = mean(cox_df$age), spec = 3, met_cox = seq(min(cox_df$met_cox), max(cox_df$met_cox), length = 200))

ggplot(pred0) + ggtitle("0") + ggplot(pred1)+ ggtitle("1")
age_quartiles <- quantile(cox_df_base$age, c( 0.25, 0.5, 0.75))


pred_list <- lapply(age_quartiles, function(age_q) {
  pred <- Predict(csp, sex = 0, age = age_q, spec = 3, met_cox = seq(min(cox_df$met_cox), max(cox_df$met_cox), length = 200))
  return(pred)
})

# Plot the predictions in three rows
plot_list <- lapply(seq_along(pred_list), function(i) {
  ggplot(pred_list[[i]]) + ggtitle(paste("Age Quartile ", i - 1))
})

# Combine the plots
gridExtra::grid.arrange(grobs = plot_list, ncol = 1)
plot(pred)
object.size(pred)

set_datadist <- function(met, disease) {
  cox_df <- 
    cox_df_base %>% 
    mutate(met_cox = surv_c[[met]],
           event = surv_c[[paste("event", disease, sep = "_")]],
           time = surv_c[[paste("time", disease, sep = "_")]])
  
  surv_obj <- with(cox_df, Surv(time = time, event = event))
  
  datadist_obj <- datadist(cox_df)
  options(datadist = "datadist_obj")
  
  return(list(cox_df, surv_obj))
}

  unnest_wider(result)
