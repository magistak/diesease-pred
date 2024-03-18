library(Hmisc)
library(tidyverse)
library(purrr)
library(furrr)
library(survival)
library(rms)
library(patchwork)
library(ggsurvfit)
library(tidycmprsk)


surv_all <- read_tsv("data/surv_all_log.tsv")
health_filters <- read_tsv("data/ukb_healthy_subpopul_full_table.tsv")
load("data/rdata/met_dis.RData")
source("R/source_functions.R")

surv_all <- surv_all |> 
  left_join(
    select(health_filters, eid, takes_medication_cholesterol, is_healthy_3, is_healthy_7)
  ) 
diseases <- disease_v
metabolites <- mets_v

dis_met <- 
  expand_grid(diseases, metabolites) |> 
  rename(disease = diseases, met = metabolites) |> 
  mutate(dm = paste(disease, met, sep = "_"))


cox_df_base <- tibble(
  spec = as.factor(surv_all$spectrometer),
  age = surv_all$age,
  sex = as.factor(surv_all$sex))

met_i <- dis_met$met[1]
dis_i <- dis_met$disease[1]
cox_df <- table_for_iteration(surv_all, cox_df_base, met_i, dis_i) |> 
  mutate(met_cox=log10(met_cox))
qs <- quantile(cox_df$met_cox, c(0.2, 0.8))

cox_df <- 
cox_df |> 
  mutate(lh=ifelse(met_cox<=quantile(cox_df$met_cox, c(0.33, 0.67))[1],"low",ifelse(met_cox>=quantile(cox_df$met_cox, c(0.33, 0.67))[2],"high","mid")))|> 
  mutate(event=as.factor(event))
cox_df_lh <- cox_df |> 
  filter(lh %in% c("low", "high")) |> 
  mutate(event = as.factor(event))

cuminc(Surv(time, event) ~ lh, data = cox_df) %>% 
  ggcuminc() + 
  labs(
    x = "Years"
  ) + 
  add_confidence_interval() +
  add_risktable()

surv_obj <- with(cox_df, Surv(time = time, event = event))


cox_df |> 
  ggplot(aes(x=met_cox))+
  geom_histogram()
quantile_breaks <- quantile(cox_df$met_cox, probs = seq(from = 0, to = 1, length.out = 11))
cox_df$value_cut <- cut(cox_df$met_cox, breaks=quantile_breaks)
incidence_table <- 
with(cox_df, table(value_cut, event)) |> 
  as_tibble() |> 
  pivot_wider(names_from = event, values_from = n) |> 
  mutate(freq=`1`/sum(`0`, `1`))

freqs <- incidence_table$freq

incidence_table |> 
  ggplot(aes(x=value_cut, y=freq))+
  geom_col()

#### freqs list ####

diseases <- diseases[diseases!="T1 Diabetes"]
dis_met <- 
  expand_grid(diseases, metabolites) |> 
  rename(disease = diseases, met = metabolites) |> 
  mutate(dm = paste(disease, met, sep = "_")) #|> 
#  head(10)

plan(multisession, workers = 6)

# test
dis_i <- dis_met$disease[1]
met_i <- dis_met$met[1]

freqs_list <- future_map2(dis_met$disease, dis_met$met, 
                          #safely(
                            function(dis_i, met_i){
  
cox_df_in <- 
tibble(met_cox = surv_all[[met_i]],
       event = surv_all[[paste("event", dis_i, sep = "_")]]) |> 
  mutate(met_cox=log10(met_cox))

cox_df_in |> 
  arrange(met_cox) |> 
  mutate(group = cut(1:nrow(cox_df_in), 10, labels = FALSE)) |> 
  group_by(group) |> 
  summarise(events = sum(event))



# with(cox_df_in, table(value_cut, event)) |> 
#   as_tibble() |> 
#   pivot_wider(names_from = event, values_from = n) |> 
#   mutate(freq=`1`/sum(`0`, `1`)) |> 
#   select(freq) |> 
#   unlist(use.names = FALSE)
}
#)
)
names(freqs_list) <- dis_met$dm

n_people <- future_map2(dis_met$disease, dis_met$met, 
          function(dis_i, met_i){
            cox_df_in <- 
            tibble(met_cox = surv_all[[met_i]],
            event = surv_all[[paste("event", dis_i, sep = "_")]])
            round(nrow(cox_df_in)/10)  
}
)

# cvs
cvs <- 
imap(
  freqs_list, function(table, name){
  with(table, sd(events)/mean(events))
  }
)
# order fl_filtered
fl_ordered <- freqs_list[order(unlist(cvs), decreasing = TRUE)]
# filter for U or anti-U shaped
u_filter <- 
imap(
  fl_ordered, function(table, name){
    v <- table$events
    ifelse(!min(v)%in%c(v[1],v[length(v)]) | !max(v)%in%c(v[1],v[length(v)]),
           c(TRUE),
           c(FALSE))
  }
)

fl_filtered <- fl_ordered[unlist(u_filter, use.names = FALSE)]

# pdf plot
i=1
pdf("plots/incidence_v3.pdf")
for(i in seq_along(fl_filtered)[1:300]){
p <-   fl_filtered[[i]] |> 
    ggplot(aes(x=as.factor(group), y=events))+
    geom_col()+
    geom_label(aes(label = signif(events/23426,2)), vjust = 1.5)+
    labs(title=names(fl_filtered)[i],
         x=paste("Number of people in each: 23 426"),
         y="Number of events",
         caption = "Frequencies on bars")
print(p)
}
dev.off()
names(freqs_list[1])



?quantile
