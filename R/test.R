library(tidycmprsk)
library(tidyverse)
library(ggsurvfit)

data(Melanoma, package = "MASS")
Melanoma <- 
  Melanoma %>% 
  mutate(
    status = as.factor(recode(status, `2` = 0, `1` = 1, `3` = 2))
  )

cuminc(Surv(time, status) ~ 1, data = Melanoma) %>% 
  ggcuminc() + 
  labs(
    x = "Days"
  ) + 
  add_confidence_interval() +
  add_risktable()

cuminc(Surv(time, status) ~ ulcer, data = Melanoma) %>% 
  ggcuminc() + 
  labs(
    x = "Days"
  ) + 
  add_confidence_interval() +
  add_risktable()
