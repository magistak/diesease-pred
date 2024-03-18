library(tidyverse)
library(purrr)
library(furrr)
library(survival)
library(rms)
library(Hmisc)
library(patchwork)
library(cowplot)

surv_all <- read_tsv("data/surv_all_log.tsv")
health_filters <- read_tsv("data/ukb_healthy_subpopul_full_table.tsv")
load("data/rdata/met_dis.RData")

surv_all <- surv_all |> 
  left_join(
  select(health_filters, eid, takes_medication_cholesterol, is_healthy_3, is_healthy_7)
) 

# two types of filters:
# healthy: no.3
# cholesterol lowering medication

# model building
# use map
# map goes through metabolite, disease pairs, for each creates a list





# functions for cox

# creates data for iteration
table_for_iteration <- function(table_all, table_base, met_i, dis_i){
    table_base %>% 
    mutate(met_cox = table_all[[met_i]],
           event = table_all[[paste("event", dis_i, sep = "_")]],
           time = table_all[[paste("time", dis_i, sep = "_")]])
}

.grobify_ggplot <- function(p) {
  if ("ggplot" %in% class(p)) {
    p <- ggplot2::ggplotGrob(p)
    p <- ggpubr::as_ggplot(p)
  }
  return(p)
}


# goal
# to plot the disease-metabolite relationships



# map2 in function with: eids, diseases, metabolites, 

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

#### run ####

disease_v1 <- disease_v[c(2,3,5,11,15,16,20,23,24,26,29,30)]
disease_v2 <- disease_v1[c(1,2,7,12)]
mets_v1 <- mets_v[c(2,3,6,7,8,9,10,18,22,26,30:63)]
mets_v2 <- mets_v[c(9,10)]

filters <- c("is_healthy_7", "takes_medication_cholesterol", "all")

plot_tables <- 
future_map(.x = filters, .f = function(healthy){
if(healthy %in% c("is_healthy_7")){
  eids <- with(surv_all, eid[get(healthy)==1])
} else if(healthy == "takes_medication_cholesterol"){
  eids <- with(surv_all, eid[get(healthy)==0])
} else {
  eids <- surv_all$eid
}
model_cox(surv_all = surv_all, eids = eids, diseases = disease_v1, metabolites = mets_v1, comment = healthy)
}
)


names(plot_tables2) <- filters
plot_tables_all <- 
map2(plot_tables, plot_tables2, function(t1,t2){
  rbind(t1,t2)
}
)

save(plot_tables_all, file = "data/plot_tables_v1.RData")

plot_tables  

# fdr-correct
plot_tables_fdr <- map(plot_tables_all, function(table){
  table |> 
    mutate(across(psl:psb, ~ p.adjust(.x,method = "fdr"))) |> 
    mutate(dis_met=paste(disease, metabolite, sep = "_"))
})

# get significant ps in order
surv_all <- read_tsv("data/surv_all_log.tsv")
surv_all <- surv_all |> 
  left_join(
    select(health_filters, eid, takes_medication_cholesterol, is_healthy_3, is_healthy_7)
  ) 


dis_met_ordered <- 
plot_tables_fdr$is_healthy_7 |> 
  filter(psl<0.1) |> 
  arrange(psl) |> 
  select(disease,met= metabolite, psl) 

#test
healthy <- "is_healthy_7"
eids <- with(surv_all, eid[get(healthy)==1])
surv_all_in = surv_all
dis_met = dis_met_ordered_test
comment = healthy
dis_i <- dis_met[[1,1]]
met_i<- dis_met[[1,2]]
psl_i<- dis_met[[1,3]]

# model
dis_met_ordered_test <- head(dis_met_ordered, 3)
preplot_model_cox <-  function(surv_all_in, eids, dis_met, comment){
  
  surv_all_pre <- surv_all_in |> 
    filter(eid %in% eids)
  
  cox_df_base <- tibble(
    spec = as.factor(surv_all_pre$spectrometer),
    age = surv_all_pre$age,
    sex = as.factor(surv_all_pre$sex))
  
  table_my <- 
    pmap(dis_met, .f = function(disease, met, psl){
      dis_i <- disease
      met_i <- met
      psl_i <- psl
      cox_df <- table_for_iteration(surv_all_pre, cox_df_base, met_i, dis_i)
      surv_obj <- with(cox_df, Surv(time = time, event = event))
      model <- cph(surv_obj ~ spec + rcs(age, 4) + rcs(met_cox, 4) + sex, data = cox_df)
  
  all_plot <- 
        tibble(
          pred=predict(model),
          met=cox_df$met_cox
        ) |> 
        ggplot(aes(x=met, y=pred))+
        geom_bin2d()+
        geom_smooth(color="red") +
        labs(title = paste(dis_i, met_i, "p:",signif(psl_i,3)),
             subtitle = paste(comment, nrow(surv_all_pre)),
             x= met_i,
             y= "risk")+
    scale_x_continuous(trans = "log10")
      
  all_plot <- .grobify_ggplot(all_plot)
      
      tibble(disease=dis_i, metabolite=met_i, plot=list(all_plot)) 
      
    }
    ) |> 
   list_rbind()
  
  return(table_my)
}
filters <- c("is_healthy_7", "takes_medication_cholesterol", "all")
load("data/plot_tables_v1.RData")
unique(plot_tables_all$all$metabolite)
preplot_tables2 <- 
  future_map(.x = filters, .f = function(healthy){
    if(healthy %in% c("is_healthy_7")){
      eids <- with(surv_all, eid[get(healthy)==1])
    } else if(healthy == "takes_medication_cholesterol"){
      eids <- with(surv_all, eid[get(healthy)==0])
    } else {
      eids <- surv_all$eid
    }
    preplot_model_cox(surv_all_in = surv_all, eids = eids, dis_met = dis_met_ordered2, comment = healthy)
  }
  )




names(preplot_tables2) <- filters

i=1
pdf(file = "plots/healthy_rev_v02.pdf",
    width = 10,
    height = 15)
for(i in 1:nrow(dis_met_ordered2)){
  p1 <- preplot_tables2$is_healthy_7$plot[[i]]
  p2 <- preplot_tables2$takes_medication_cholesterol$plot[[i]]
  p3 <- preplot_tables2$all$plot[[i]]

  p <-  p1 / p2 / p3
  print(p)
}
dev.off()



preplot_tables[[1]]$plot[1]


# print
map(dis_met_ordered,function(dis_met){
  
  #model
  
  #plot
  
  #add plots
  
  
  
})





p_all <- 
  all_plot |> 
  ggplot(aes(x=met, y=pred))+
  geom_bin2d()+
  geom_smooth(color="red") +
  labs(title = paste(dis_i, met_i),
       x= met_i,
       y= "risk")

print(p_all)



all_plot <- 
  tibble(
    pred=predict(model),
    met=cox_df$met_cox
  ) |> 
  ggplot(aes(x=met, y=pred))+
  geom_bin2d()+
  geom_smooth(color="red") +
  labs(title = paste(dis_i, met_i),
       subtitle = comment,
       x= met_i,
       y= "risk") 
all_plot <- .grobify_ggplot(all_plot)
# viridis options=mako scale fill

# analysis
table <- plot_tables_fdr[[1]]
signif_table <- 
map(plot_tables_fdr,function(table){
table |> 
  mutate(across(psl:psb,~.x<0.1)) 
}
)

dis_met2 <- setdiff(with(signif_table$all, dis_met[psl]),with(signif_table$is_healthy_7, dis_met[psl]))
dis_met_ordered2 <- 
plot_tables_fdr$all |> 
  filter(dis_met %in% dis_met2) |> 
  arrange(desc(csl)) |> 
  select(disease,met= metabolite, psl)
