# Insulin Sensitivity
#
# Summary of relationship between M-value and various insulin sensitivity indices.
#
# Summarizes population, mean, SD, and effect size for baseline, endpoint and change values.
# Creates results for normal and logged data; for TRT_A only and for all treatments combined. 
# Above done for baseline, endpoint, and change values.
# 
# Data was created using rnorm(nrow(), mean=1, sd=.3) on real data and filtered down to variables of interest for this particular project.
#


source("H:/R/InsulinSensitivity/R/Indices_Data.R")


### Summary of Variables
detach("package:plyr", unload=TRUE) 
summary_trt_a <- full_set %>% group_by(Parameter) %>% summarise( n=n(),mean_base= mean(Baseline), sd_base =sd(Baseline), effect_size_base = mean(Baseline)/sd(Baseline),
                                                           mean_end= mean(Endpoint), sd_end =sd(Endpoint),effect_size_end = mean(Endpoint)/sd(Endpoint),
                                                           mean_change= mean(Change), sd_change =sd(Change),effect_size_change = mean(Change)/sd(Change))
summary <- all_trt %>% group_by(Parameter) %>% summarise( n=n(),mean_base= mean(Baseline), sd_base =sd(Baseline), effect_size_base = mean(Baseline)/sd(Baseline),
                                                           mean_end= mean(Endpoint), sd_end =sd(Endpoint),effect_size_end = mean(Endpoint)/sd(Endpoint),
                                                           mean_change= mean(Change), sd_change =sd(Change),effect_size_change = mean(Change)/sd(Change))
summary_by_trt <- all_trt %>% group_by(Parameter,Treatment) %>% summarise( n=n(),mean_base= mean(Baseline), sd_base =sd(Baseline), effect_size_base = mean(Baseline)/sd(Baseline),
                                                                           mean_end= mean(Endpoint), sd_end =sd(Endpoint),effect_size_end = mean(Endpoint)/sd(Endpoint),
                                                                           mean_change= mean(Change), sd_change =sd(Change),effect_size_change = mean(Change)/sd(Change))


write.csv(summary_trt_a, "InsulinSensitivity_trt_a_summary.csv")
write.csv(summary, "InsulinSensitivity_summary.csv")
write.csv(summary__by_trt, "InsulinSensitivity_summary_byTRT.csv")


### Split Summaries (trt_a only)
summary_baseline <- full_set %>% group_by(Parameter) %>% summarise(n=n(), Baseline_Means= mean(Baseline), sd_base =sd(Baseline), effect_size_base = mean(Baseline)/sd(Baseline),)
summary_change <-  full_set %>% dplyr::group_by(Parameter) %>%dplyr::summarize(n=n(), Change_Means= mean(Change), SD_change =sd(Change), effect_size_change = mean(Change)/sd(Change))
summary_endpoint <-  full_set %>% dplyr::group_by(Parameter) %>%dplyr::summarize(n=n(), Endpoint_Means= mean(Endpoint), SD_endpoint =sd(Endpoint), effect_size_endpoint = mean(Endpoint)/sd(Endpoint))

summary_baseline_comp <-  baseline_compins_trt_a %>%dplyr::summarize(n=n(), Baseline_Means= mean(Composite_Insulin), sd_base =sd(Composite_Insulin), effect_size_base = mean(Composite_Insulin)/sd(Composite_Insulin))
summary_baseline_comp$Parameter <- "Mastuda Index"
summary_baseline1 <- merge(summary_baseline, summary_baseline_comp,all=TRUE)

summary_change_comp <-  change_compins_trt_a %>%dplyr::summarize(n=n(), Change_Means= mean(Composite_Insulin), SD_change =sd(Composite_Insulin), effect_size_change = mean(Composite_Insulin)/sd(Composite_Insulin))
summary_change_comp$Parameter <- "Mastuda Index"
summary_change1 <- merge(summary_change, summary_change_comp,all=TRUE)

summary_endpoint_comp <-  endpoint_compins_trt_a %>%dplyr::summarize(n=n(), Endpoint_Means= mean(Composite_Insulin), SD_endpoint =sd(Composite_Insulin), effect_size_endpoint = mean(Composite_Insulin)/sd(Composite_Insulin))
summary_endpoint_comp$Parameter <- "Mastuda Index"
summary_endpoint1 <- merge(summary_endpoint, summary_endpoint_comp,all=TRUE)


### Split Summaries (All Treatments)
summary_baseline <- all_trt %>% group_by(Parameter) %>% summarise(n=n(), Baseline_Means= mean(Baseline), sd_base =sd(Baseline), effect_size_base = mean(Baseline)/sd(Baseline),)
summary_change <-  all_trt %>% dplyr::group_by(Parameter) %>%dplyr::summarize(n=n(), Change_Means= mean(Change), SD_change =sd(Change), effect_size_change = mean(Change)/sd(Change))
summary_endpoint <-  all_trt %>% dplyr::group_by(Parameter) %>%dplyr::summarize(n=n(), Endpoint_Means= mean(Endpoint), SD_endpoint =sd(Endpoint), effect_size_endpoint = mean(Endpoint)/sd(Endpoint))

summary_baseline_comp <-  baseline_compins_all %>%dplyr::summarize(n=n(), Baseline_Means= mean(Composite_Insulin), sd_base =sd(Composite_Insulin), effect_size_base = mean(Composite_Insulin)/sd(Composite_Insulin))
summary_baseline_comp$Parameter <- "Mastuda Index"
summary_baseline1 <- merge(summary_baseline, summary_baseline_comp,all=TRUE)

summary_change_comp <-  change_compins_all %>%dplyr::summarize(n=n(), Change_Means= mean(Composite_Insulin), SD_change =sd(Composite_Insulin), effect_size_change = mean(Composite_Insulin)/sd(Composite_Insulin))
summary_change_comp$Parameter <- "Mastuda Index"
summary_change1 <- merge(summary_change, summary_change_comp,all=TRUE)

summary_endpoint_comp <-  endpoint_compins_all %>%dplyr::summarize(n=n(), Endpoint_Means= mean(Composite_Insulin), SD_endpoint =sd(Composite_Insulin), effect_size_endpoint = mean(Composite_Insulin)/sd(Composite_Insulin))
summary_endpoint_comp$Parameter <- "Mastuda Index"
summary_endpoint1 <- merge(summary_endpoint, summary_endpoint_comp,all=TRUE)



#### Log Data
source("H:/R/InsulinSensitivity/R/Indices_Ln_Data.R")

## Log Summary
summary_ln_trt_a <- lnfull_set %>% group_by(Parameter) %>% summarise( n=n(),mean_base= mean(Baseline), sd_base =sd(Baseline), effect_size_base = mean(Baseline)/sd(Baseline),
                                                                mean_end= mean(Endpoint), sd_end =sd(Endpoint),effect_size_end = mean(Endpoint)/sd(Endpoint),
                                                                mean_change= mean(Change), sd_change =sd(Change),effect_size_change = mean(Change)/sd(Change))
summary_ln <- all_trt_ln %>% group_by(Parameter) %>% summarise( n=n(),mean_base= mean(Baseline), sd_base =sd(Baseline), effect_size_base = mean(Baseline)/sd(Baseline),
                                                                mean_end= mean(Endpoint), sd_end =sd(Endpoint),effect_size_end = mean(Endpoint)/sd(Endpoint),
                                                                mean_change= mean(Change), sd_change =sd(Change),effect_size_change = mean(Change)/sd(Change))
summary_by_trt_ln <- all_trt_ln %>% group_by(Parameter,Treatment) %>% summarise( n=n(),mean_base= mean(Baseline), sd_base =sd(Baseline), effect_size_base = mean(Baseline)/sd(Baseline),
                                                                                 mean_end= mean(Endpoint), sd_end =sd(Endpoint),effect_size_end = mean(Endpoint)/sd(Endpoint),
                                                                                 mean_change= mean(Change), sd_change =sd(Change),effect_size_change = mean(Change)/sd(Change))


write.csv(summary_trt_a, "ln_InsulinSensitivity_trt_a_summary.csv")
write.csv(summary, "ln_InsulinSensitivity_summary.csv")
write.csv(summary__by_trt, "ln_InsulinSensitivity_summary_byTRT.csv")


### Split Summaries (trt_a only)
ln_summary_baseline <- lnfull_set %>% group_by(Parameter) %>% summarise(n=n(), Baseline_Means= mean(Baseline), sd_base =sd(Baseline), effect_size_base = mean(Baseline)/sd(Baseline),)
ln_summary_change <-  lnfull_set %>% dplyr::group_by(Parameter) %>%dplyr::summarize(n=n(), Change_Means= mean(Change), SD_change =sd(Change), effect_size_change = mean(Change)/sd(Change))
ln_summary_endpoint <-  lnfull_set %>% dplyr::group_by(Parameter) %>%dplyr::summarize(n=n(), Endpoint_Means= mean(Endpoint), SD_endpoint =sd(Endpoint), effect_size_endpoint = mean(Endpoint)/sd(Endpoint))

ln_summary_baseline_comp <-  baseline_compins_trt_a %>%dplyr::summarize(n=n(), Baseline_Means= mean(Composite_Insulin), sd_base =sd(Composite_Insulin), effect_size_base = mean(Composite_Insulin)/sd(Composite_Insulin))
ln_summary_baseline_comp$Parameter <- "log(Mastuda Index)"
ln_summary_baseline1 <- merge(ln_summary_baseline, ln_summary_baseline_comp,all=TRUE)

ln_summary_change_comp <-  change_compins_trt_a %>%dplyr::summarize(n=n(), Change_Means= mean(Composite_Insulin), SD_change =sd(Composite_Insulin), effect_size_change = mean(Composite_Insulin)/sd(Composite_Insulin))
ln_summary_change_comp$Parameter <- "log(Mastuda Index)"
ln_summary_change1 <- merge(ln_summary_change, ln_summary_change_comp,all=TRUE)

ln_summary_endpoint_comp <-  endpoint_compins_trt_a %>%dplyr::summarize(n=n(), Endpoint_Means= mean(Composite_Insulin), SD_endpoint =sd(Composite_Insulin), effect_size_endpoint = mean(Composite_Insulin)/sd(Composite_Insulin))
ln_summary_endpoint_comp$Parameter <- "log(Mastuda Index)"
ln_summary_endpoint1 <- merge(ln_summary_endpoint, ln_summary_endpoint_comp,all=TRUE)


### Split Summaries (All Treatments)
ln_summary_baseline <- all_trt_ln %>% group_by(Parameter) %>% summarise(n=n(), Baseline_Means= mean(Baseline), sd_base =sd(Baseline), effect_size_base = mean(Baseline)/sd(Baseline),)
ln_summary_change <-  all_trt_ln %>% dplyr::group_by(Parameter) %>%dplyr::summarize(n=n(), Change_Means= mean(Change), SD_change =sd(Change), effect_size_change = mean(Change)/sd(Change))
ln_summary_endpoint <-  all_trt_ln %>% dplyr::group_by(Parameter) %>%dplyr::summarize(n=n(), Endpoint_Means= mean(Endpoint), SD_endpoint =sd(Endpoint), effect_size_endpoint = mean(Endpoint)/sd(Endpoint))

ln_summary_baseline_comp <-  baseline_compins_all %>%dplyr::summarize(n=n(), Baseline_Means= mean(Composite_Insulin), sd_base =sd(Composite_Insulin), effect_size_base = mean(Composite_Insulin)/sd(Composite_Insulin))
ln_summary_baseline_comp$Parameter <- "log(Mastuda Index)"
ln_summary_baseline1 <- merge(ln_summary_baseline, ln_summary_baseline_comp,all=TRUE)

ln_summary_change_comp <-  change_compins_all %>%dplyr::summarize(n=n(), Change_Means= mean(Composite_Insulin), SD_change =sd(Composite_Insulin), effect_size_change = mean(Composite_Insulin)/sd(Composite_Insulin))
ln_summary_change_comp$Parameter <- "log(Mastuda Index)"
ln_summary_change1 <- merge(ln_summary_change, ln_summary_change_comp,all=TRUE)

ln_summary_endpoint_comp <-  endpoint_compins_all %>%dplyr::summarize(n=n(), Endpoint_Means= mean(Composite_Insulin), SD_endpoint =sd(Composite_Insulin), effect_size_endpoint = mean(Composite_Insulin)/sd(Composite_Insulin))
ln_summary_endpoint_comp$Parameter <- "log(Mastuda Index)"
ln_summary_endpoint1 <- merge(ln_summary_endpoint, ln_summary_endpoint_comp,all=TRUE)










