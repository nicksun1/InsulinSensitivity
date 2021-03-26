# Insulin Sensitivity
#
# Evaluates the correlation between commmonly used insulin sensitivity indices with the M-value (gold standard) from euglycemic clamp test ()
#
# Calculates Mastuda Index (composite insulin)
# Outputs correlation and corresponding p-value between M-value and various insulin sensitivity indices
# Outputs Regression graphs of M-value against various insulin sensitivity indices with significance and variability factors (p-value and R^2)
# Summary of analysis data with N, Mean, SD, and effect size for by treatment and combined.
# 
# Data was created using rnorm(nrow(), mean=1, sd=.3) on real data and filtered down to variables of interest for this particular project.
# *Library coastr is an internal package with the sole purpose of increasing the efficiency and speed of pulling data from the internal server
#



#suppressMessages(library(coastr))
suppressMessages(library(readr))
#library(plyr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(psych)
library(gridExtra)

#####################

new_params <- read.csv("gway_IS_parameters.csv")
clamp <- read.csv("clamp.csv")
efficacy <- read.csv("efficacy.csv")
testmeal <- read.csv("testmeal.csv")


# Filter out non complete data from csv (SI, Quicki, OGIS, logPREDIM)
non_dup <-  data.frame(new_params[!duplicated(new_params$PATIENT),])
sort_dup <- new_params[order(new_params$PATIENT, -new_params$VISIT),]
dup <- data.frame(sort_dup[duplicated(new_params$PATIENT),])
diff <- setdiff(non_dup, dup)
final_new_params <- setdiff(new_params,diff)
final_new_params$THERAPY <- as.character(final_new_params$THERAPY)
final_new_params$THERAPY[final_new_params$THERAPY=="ROS"]<- "Rosiglitazone"
final_new_params$THERAPY[final_new_params$THERAPY=="EXENROSI"]<- "Exenatide+Rosiglt" 
final_new_params$THERAPY[final_new_params$THERAPY=="EXENATI1"]<- "Exenatide" 
long_newparams <- final_new_params


# SI
si_long <- long_newparams %>% select(PATIENT, VISIT, THERAPY, SI)
si <- reshape(si_long, idvar="PATIENT", timevar = "VISIT", direction ="wide")
si$Change <- si$SI.7-si$SI.2
si$Parameter <- "SI"
si<- si[,-4]
names(si) <- c("Subject ID","Treatment", "Baseline", "Endpoint", "Change", "Parameter")

# Quicki
quicki_long <- long_newparams %>% select(PATIENT, VISIT, THERAPY, QUICKI)
quicki <- reshape(quicki_long, idvar="PATIENT", timevar = "VISIT", direction ="wide")
quicki$Change <- quicki$QUICKI.7-quicki$QUICKI.2
quicki$Parameter <- "QUICKI"
quicki<- quicki[,-4]
names(quicki) <- c("Subject ID","Treatment", "Baseline", "Endpoint", "Change", "Parameter")

# OGIS
ogis_long <- long_newparams %>% select(PATIENT, VISIT, THERAPY, OGIS)
ogis <- reshape(ogis_long, idvar="PATIENT", timevar = "VISIT", direction ="wide")
ogis$Change <- ogis$OGIS.7-ogis$OGIS.2
ogis$Parameter <- "OGIS"
ogis<- ogis[,-4]
names(ogis) <- c("Subject ID","Treatment", "Baseline", "Endpoint", "Change", "Parameter")

# logPREDIM
predim_long <- long_newparams%>% select(PATIENT, VISIT, THERAPY, logPREDIM)
predim <- reshape(predim_long, idvar="PATIENT", timevar = "VISIT", direction ="wide")
predim$Change <- predim$logPREDIM.7-predim$logPREDIM.2
predim$Parameter <- "PREDIM"
predim<- predim[,-4]
names(predim) <- c("Subject ID","Treatment", "Baseline", "Endpoint", "Change", "Parameter")
predim$Baseline <- exp(predim$Baseline)
predim$Endpoint <- exp(predim$Endpoint)
predim$Change <- exp(predim$Endpoint)-exp(predim$Baseline)

# Organize and filter data from sas datasets
### M-value
m_value_full <- clamp %>% select (PATIENT, TRT, PARACODE, VAL_B, VAL_E, VAL_C)
m_value_full <- m_value_full %>% filter(PARACODE == "M_VALUE")
m_value <- m_value_full[!duplicated(m_value_full[,c('PATIENT','TRT','PARACODE')]),]
m_value <- m_value[complete.cases(m_value),]
names(m_value) <- c("Subject ID", "Treatment", "Parameter", "Baseline","Endpoint", "Change")

unit_check <- efficacy %>% select (PATIENT, VISIT, TRT, PARACODE,  VAL_B,VAL, VAL_C, UNIT)
unit_check <- unit_check %>% filter(PARACODE == "HOMASINS"|PARACODE == "TINSUL"|PARACODE == "FGLU")

### Homa_S
homa_s_full <- efficacy %>% select (PATIENT, VISIT, TRT, PARACODE,  VAL_B,VAL, VAL_C)
homa_s_full <- homa_s_full %>% filter(PARACODE == "HOMASINS")
homa_s <- homa_s_full[complete.cases(homa_s_full),]
homa_s <- homa_s %>% group_by(PATIENT) %>%
  top_n(n=1, wt= VISIT)
names(homa_s) <- c("Subject ID","Visit", "Treatment", "Parameter", "Baseline","Endpoint", "Change")

### Fasting Insulin
fasting_insulin_full <- efficacy %>% select (PATIENT, VISIT, TRT, PARACODE, VAL, VAL_B, VAL_C)
fasting_insulin_full <- fasting_insulin_full %>% filter(PARACODE == "TINSUL") 
fasting_insulin <- fasting_insulin_full[complete.cases(fasting_insulin_full),]
fasting_insulin <- fasting_insulin %>% group_by(PATIENT) %>%
  top_n(n=1, wt= VISIT)
names(fasting_insulin) <- c("Subject ID","Visit", "Treatment", "Parameter", "Endpoint", "Baseline", "Change")

### Fasting Glucose
fast_glucose_full <-  efficacy %>% select (PATIENT, VISIT, TRT, PARACODE, VAL, VAL_B, VAL_C)
fast_glucose_full <- fast_glucose_full %>% filter(PARACODE == "FGLU")
fast_glucose <- fast_glucose_full[complete.cases(fast_glucose_full),]
fast_glucose <- fast_glucose %>% group_by(PATIENT) %>%
  top_n(n=1, wt= VISIT)
names(fast_glucose) <- c("Subject ID","Visit", "Treatment", "Parameter", "Endpoint", "Baseline", "Change")

### Insulin AUC
insulin_auc_full <-  testmeal %>% select (PATIENT, VISIT, TRT, ANALYTE, PARACODE, VAL, VAL_B, VAL_E, VAL_C)
insulin_auc_full <-  insulin_auc_full %>% filter(PARACODE == "AUC") 
insulin_auc_full <-  insulin_auc_full %>% filter(ANALYTE =="TINSUL")
insulin_auc_full$PARACODE <- "TINSUL-AUC"
insulin_auc_full$ANALYTE <- NULL
insulin_auc_full$VAL <- NULL
insulin_auc <- insulin_auc_full[complete.cases(insulin_auc_full),]
insulin_auc <- insulin_auc %>% group_by(PATIENT) %>%
  top_n(n=1, wt= VISIT)
names(insulin_auc) <- c("Subject ID","Visit", "Treatment", "Parameter", "Baseline", "Endpoint", "Change")

### Glucose AUC
glucose_auc_full <-  testmeal %>% select (PATIENT, VISIT, TRT, ANALYTE, PARACODE, VAL, VAL_B, VAL_E, VAL_C)
glucose_auc_full <-  glucose_auc_full %>% filter(PARACODE == "AUC") 
glucose_auc_full <-  glucose_auc_full %>% filter(ANALYTE =="FGLU")
glucose_auc_full$PARACODE <- "FGLU-AUC"
glucose_auc_full$ANALYTE <- NULL
glucose_auc_full$VAL <- NULL
glucose_auc <- glucose_auc_full[complete.cases(glucose_auc_full),]
glucose_auc <- glucose_auc %>% group_by(PATIENT) %>%
  top_n(n=1, wt= VISIT)
names(glucose_auc) <- c("Subject ID","Visit", "Treatment", "Parameter", "Baseline", "Endpoint", "Change")

#### Merge full sets into one
full_set_testmeal <- merge(insulin_auc, glucose_auc, all=TRUE)
full_set_testmeal <- merge(full_set_testmeal, m_value, all=TRUE)
full_set_efficacy <- merge(fast_glucose, fasting_insulin, all=TRUE)
full_set <- merge(full_set_testmeal, full_set_efficacy, all=TRUE)
full_set <- merge(full_set, homa_s, all=TRUE)
full_set <- merge(full_set, si, all=TRUE)
full_set <- merge(full_set, ogis, all=TRUE)
full_set <- merge(full_set, predim, all=TRUE)
full_set <- merge(full_set, quicki, all=TRUE)
full_set <- subset(full_set, select=-c(Visit))
all_trt <- full_set
full_set <- full_set %>% filter (Treatment=="Rosiglitazone")
full_set_complete <- full_set[complete.cases(full_set),]


# Calculate Composite Insulin 
baseline_compins <- full_set[,c(1,3:4)]
baseline_compins <- reshape(baseline_compins, idvar="Subject ID", timevar="Parameter", direction="wide")
baseline_compins$`Subject ID` <- NULL
baseline_compins$Composite_Insulin <- 10000/(sqrt(baseline_compins[,1]*baseline_compins[,9]*baseline_compins[,2]/165*baseline_compins[,10]/165))
baseline_compins <- baseline_compins[complete.cases(baseline_compins[,11]),]

endpoint_compins <- full_set[,c(1,3,5)]
endpoint_compins <- reshape(endpoint_compins, idvar="Subject ID", timevar="Parameter", direction="wide")
endpoint_compins$`Subject ID` <- NULL
endpoint_compins$Composite_Insulin <- 10000/(sqrt(endpoint_compins[,1]*endpoint_compins[,9]*endpoint_compins[,2]/165*endpoint_compins[,10]/165))
endpoint_compins <- endpoint_compins[complete.cases(endpoint_compins[,11]),]

change_compins <- full_set[,c(1,3,6)]
change_compins <- reshape(change_compins, idvar="Subject ID", timevar="Parameter", direction="wide")
change_compins$`Subject ID` <- NULL
change_compins <- change_compins[complete.cases(change_compins[,c(1,2,9,10)]),]
change_compins$Composite_Insulin <- endpoint_compins$Composite_Insulin-baseline_compins$Composite_Insulin

# Output summary values
summary_change <-  full_set_complete %>% dplyr::group_by(Parameter) %>%dplyr::summarize(n=n(), Change_Means= mean(Change), SD_change =sd(Change), effect_size_change = mean(Change)/sd(Change))
write.csv(summary_change, "Log_GWAY_insulin_sensitivity_summary_change.csv")

summary_change_comp <-  change_compins %>%dplyr::summarize(n=33, Change_Means= mean(Composite_Insulin), SD_change =sd(Composite_Insulin), effect_size_change = mean(Composite_Insulin)/sd(Composite_Insulin))
summary_change_comp$Parameter <- "log(Composite_Insulin)"
summary_change1 <- merge(summary_change, summary_change_comp,all=TRUE)
write.csv(summary_change1, "GWAY_ROS_summary_change.csv")


library(plyr)
## Correlation Matrices
# Create data to do correlation on
baseline_corr <- full_set[,c(1,3:4)]
baseline_corr <- reshape(baseline_corr, idvar="Subject ID", timevar="Parameter", direction="wide")
baseline_corr$`Subject ID` <- NULL
baseline_corr <- baseline_corr[complete.cases(baseline_corr),]
baseline_corr$Composite_Insulin <- 10000/(sqrt(baseline_corr[,1]*baseline_corr[,9]*baseline_corr[,2]/165*baseline_corr[,10]/165))

endpoint_corr <- full_set[,c(1,3,5)]
endpoint_corr <- reshape(endpoint_corr, idvar="Subject ID", timevar="Parameter", direction="wide")
endpoint_corr$`Subject ID` <- NULL
endpoint_corr <- endpoint_corr[complete.cases(endpoint_corr),]
endpoint_corr$Composite_Insulin <- 10000/(sqrt(endpoint_corr[,1]*endpoint_corr[,9]*endpoint_corr[,2]/165*endpoint_corr[,10]/165))

change_corr <- full_set[,c(1,3,6)]
change_corr <- reshape(change_corr, idvar="Subject ID", timevar="Parameter", direction="wide")
change_corr$`Subject ID` <- NULL
change_corr <- change_corr[complete.cases(change_corr),]
change_corr$Composite_Insulin <- baseline_corr$Composite_Insulin-endpoint_corr$Composite_Insulin

#####  Add Composite Insulin
names(baseline_corr) <- c("FGLU", "FGLU_AUC", "HOMA_S", "M_VALUE","OGIS","PREDIM","QUICKI","SI", "TINSUL","TINSUL_AUC","Composite_Insulin")
names(endpoint_corr) <- c("FGLU", "FGLU_AUC", "HOMA_S", "M_VALUE","OGIS","PREDIM","QUICKI","SI", "TINSUL","TINSUL_AUC","Composite_Insulin")
names(change_corr) <- c("FGLU", "FGLU_AUC", "HOMA_S", "M_VALUE","OGIS","PREDIM","QUICKI","SI", "TINSUL","TINSUL_AUC","Composite_Insulin")
  

###################### Make Correlation Matrix
x <- data.frame(baseline_corr$M_VALUE)
y <- subset(baseline_corr, select=-c(M_VALUE))
temp_pvalue_base <- corr.test(x,y, adjust="none")
baseline_pvalue <- data.frame(temp_pvalue_base[[4]])
correlation_for_baseline <- data.frame(cor(x,y))

w <- data.frame(endpoint_corr$M_VALUE)
z <- subset(endpoint_corr, select=-c(M_VALUE))
temp_pvalue <- corr.test(w,z, adjust="none")
endpoint_pvalue <- data.frame(temp_pvalue[[4]])
correlation_for_endpoint <- data.frame(cor(w,z))

t <- data.frame(change_corr$M_VALUE)
v <- subset(change_corr, select=-c(M_VALUE))
temp_pvalue_change <- corr.test(t,v, adjust="none")
change_pvalue <- data.frame(temp_pvalue_change[[4]])
correlation_for_change <- data.frame(cor(t,v))

correlation_final <- merge(correlation_for_baseline, correlation_for_endpoint, all= TRUE)
correlation_final <- merge(correlation_final, correlation_for_change, all= TRUE)
row.names(correlation_final) <- c("Change","Baseline","Endpoint")
correlation_final <- correlation_final[c(2,3,1),]

pvalue_final <- merge(baseline_pvalue,endpoint_pvalue,all=TRUE)
pvalue_final <- merge(pvalue_final,change_pvalue,all=TRUE)
row.names(pvalue_final) <- c("Endpoint","Change","Baseline")
pvalue_final <- pvalue_final[c(3,1,2),]

write.csv(pvalue_final, "GWAY_ROS_pvalues.csv")
write.csv(correlation_final, "GWAY_ROS_correlation.csv")




# Labels creation for pvalue/r
## Baseline
baseline_pvalue1 <- cbind(q=0, baseline_pvalue)
correlation_for_baseline1 <- cbind(u=0, correlation_for_baseline)
lm_eqn_base <- function(a){
  eq <- substitute(italic(p)~"="~x *","~italic(r)~"="~y, 
                   list(x = format(baseline_pvalue1[,a], digits = 3),
                        y= format(correlation_for_baseline1[,a], digits = 3)))
  as.character(as.expression(eq));
}
labels = data.frame(x = 3, y = 3, label = lm_eqn_base(1))

## Endpoint
endpoint_pvalue1 <- cbind(q=0, endpoint_pvalue)
correlation_for_endpoint1 <- cbind(u=0, correlation_for_endpoint)
lm_eqn_end <- function(a){
  eq <- substitute(italic(p)~"="~x *","~italic(r)~"="~y, 
                   list(x = format(endpoint_pvalue1[,a], digits = 3),
                        y= format(correlation_for_endpoint1[,a], digits = 3)))
  as.character(as.expression(eq));
}
labels = data.frame(x = 3, y = 3, label = lm_eqn_end(1))

## Change
change_pvalue1 <- cbind(q=0, change_pvalue)
correlation_for_change1 <- cbind(u=0, correlation_for_change)
lm_eqn_change <- function(a){
  eq <- substitute(italic(p)~"="~x *","~italic(r)~"="~y, 
                   list(x = format(change_pvalue1[,a], digits = 3),
                        y= format(correlation_for_change1[,a], digits = 3)))
  as.character(as.expression(eq));
}
labels = data.frame(x = 3, y = 3, label = lm_eqn_change(1))


### Make Graphs
# Baseline graphs
b<- list()
baseline_corr1 <- baseline_corr[,c(4,1,2,3,5,6,7,8,9,10)]
for(i in 2:length(baseline_corr1)){
  b1= ggplot(baseline_corr1, aes(x=M_VALUE, y=baseline_corr1[,i])) + geom_point()+
    geom_smooth(method='lm') + geom_text(data = labels, aes(x = max(baseline_corr1$M_VALUE)-2, y = max(baseline_corr1[,i]), label = lm_eqn_base(i)), parse = TRUE) +
    ggtitle(paste("Baseline Comparison for M-value v",colnames(baseline_corr1)[i])) + ylab(colnames(baseline_corr1)[i])
  b[[i]] <- b1
}
pdf("Baseline_plots.pdf")
for(i in 2:length(baseline_corr1)){
  print(b[[i]])
}
dev.off()


#################### Endpoint graphs

e<- list()
endpoint_corr1 <- endpoint_corr[,c(4,1,2,3,5,6,7,8,9,10)]
for(i in 2:length(endpoint_corr1)){
  e1= ggplot(endpoint_corr1, aes(x=M_VALUE, y=endpoint_corr1[,i])) + geom_point()+
    geom_smooth(method='lm') + geom_text(data = labels, aes(x = max(endpoint_corr1$M_VALUE)-2, y = max(endpoint_corr1[,i]), label = lm_eqn_end(i)), parse = TRUE) +
    ggtitle(paste("Baseline Comparison for M-value v",colnames(endpoint_corr1)[i])) + ylab(colnames(endpoint_corr1)[i])
  e[[i]] <- e1
}
pdf("Endpoint_plots.pdf")
for(i in 2:length(endpoint_corr1)){
  print(e[[i]])
}
dev.off()


#################### Change graphs

c<- list()
change_corr1 <- change_corr[,c(4,1,2,3,5,6,7,8,9,10)]
for(i in 2:length(change_corr1)){
  c1= ggplot(change_corr1, aes(x=M_VALUE, y=change_corr1[,i])) + geom_point()+
    geom_smooth(method='lm') + geom_text(data = labels, aes(x = max(change_corr1$M_VALUE)-2, y = max(change_corr1[,i]), label = lm_eqn_change(i)), parse = TRUE) +
    ggtitle(paste("Baseline Comparison for M-value v",colnames(change_corr1)[i])) + ylab(colnames(change_corr1)[i])
  c[[i]] <- c1
}
pdf("Change_plots.pdf")
for(i in 2:length(change_corr1)){
  print(c[[i]])
}
dev.off()


### Summary of Variables
detach("package:plyr", unload=TRUE) 
full_set_baseline <- full_set[,1:4]
summary_baseline <- full_set %>% group_by(Parameter) %>% summarise(N=n(), Baseline_Means= mean(Baseline), sd_base =sd(Baseline), effect_size_base = mean(Baseline)/sd(Baseline),)

baseline_compins <- full_set[,c(1,3:4)]
baseline_compins <- reshape(baseline_compins, idvar="Subject ID", timevar="Parameter", direction="wide")
baseline_compins$`Subject ID` <- NULL
names(baseline_compins) <- c("FGLU", "FGLU_AUC", "HOMA_S", "M_VALUE","OGIS","PREDIM","QUICKI","SI", "TINSUL","TINSUL_AUC")
baseline_compins$Composite_Insulin <- 10000/(sqrt(baseline_compins$FGLU*baseline_compins$TINSUL*baseline_compins$FGLU_AUC/165*baseline_compins$TINSUL_AUC/165))
baseline_compins <- baseline_compins[complete.cases(baseline_compins[,10]),]
summary_baseline_comp <-  baseline_compins %>%dplyr::summarize(n=33, Baseline_Means= mean(Composite_Insulin), SD_base =sd(Composite_Insulin), effect_size_base = mean(Composite_Insulin)/sd(Composite_Insulin))


summary <- full_set %>% group_by(Parameter) %>% summarise( n=n(),mean_base= mean(Baseline), sd_base =sd(Baseline), effect_size_base = mean(Baseline)/sd(Baseline),
                                                           mean_end= mean(Endpoint), sd_end =sd(Endpoint),effect_size_end = mean(Endpoint)/sd(Endpoint),
                                                           mean_change= mean(Change), sd_change =sd(Change),effect_size_change = mean(Change)/sd(Change))
summary_by_trt <- all_trt %>% group_by(Parameter,Treatment) %>% summarise( n=n(),mean_base= mean(Baseline), sd_base =sd(Baseline), effect_size_base = mean(Baseline)/sd(Baseline),
                                                                            mean_end= mean(Endpoint), sd_end =sd(Endpoint),effect_size_end = mean(Endpoint)/sd(Endpoint),
                                                                            mean_change= mean(Change), sd_change =sd(Change),effect_size_change = mean(Change)/sd(Change))
write.csv(summary_baseline, "GWAY_insulin_sensitivity_baseline_summary.csv")
write.csv(summary_by_trt, "GWAY_insulin_sensitivity_summary_by_treatment.csv")



###  Summary of Variables for all Treatments

full_set_baseline <- all_trt[,1:4]
summary_baseline <- all_trt %>% group_by(Parameter) %>% summarise(N=n(), Baseline_Means= mean(Baseline), sd_base =sd(Baseline), effect_size_base = mean(Baseline)/sd(Baseline),)

baseline_compins <- all_trt[,c(1,3:4)]
baseline_compins <- reshape(baseline_compins, idvar="Subject ID", timevar="Parameter", direction="wide")
baseline_compins$`Subject ID` <- NULL
names(baseline_compins) <- c("FGLU", "FGLU_AUC", "HOMA_S", "M_VALUE","OGIS","PREDIM","QUICKI","SI", "TINSUL","TINSUL_AUC")
baseline_compins$Composite_Insulin <- 10000/(sqrt(baseline_compins$FGLU*baseline_compins$TINSUL*baseline_compins$FGLU_AUC/165*baseline_compins$TINSUL_AUC/165))
baseline_compins <- baseline_compins[complete.cases(baseline_compins[,10]),]
summary_baseline_comp <-  baseline_compins %>%dplyr::summarize(n=101, Baseline_Means= mean(Composite_Insulin), SD_base =sd(Composite_Insulin), effect_size_base = mean(Composite_Insulin)/sd(Composite_Insulin))

write.csv(summary_baseline, "GWAY_insulin_sensitivity_baseline_summary_alltrt.csv")



## Correlation and P-Value for all treatments
all_trt <- all_trt[complete.cases(all_trt),]

baseline_corr <- all_trt[,c(1,3:4)]
baseline_corr <- reshape(baseline_corr, idvar="Subject ID", timevar="Parameter", direction="wide")
baseline_corr$`Subject ID` <- NULL
baseline_corr <- baseline_corr[complete.cases(baseline_corr),]
baseline_corr$Composite_Insulin <- 10000/(sqrt(baseline_corr[,1]*baseline_corr[,9]*baseline_corr[,2]/165*baseline_corr[,10]/165))

endpoint_corr <- all_trt[,c(1,3,5)]
endpoint_corr <- reshape(endpoint_corr, idvar="Subject ID", timevar="Parameter", direction="wide")
endpoint_corr$`Subject ID` <- NULL
endpoint_corr <- endpoint_corr[complete.cases(endpoint_corr),]
endpoint_corr$Composite_Insulin <- 10000/(sqrt(endpoint_corr[,1]*endpoint_corr[,9]*endpoint_corr[,2]/165*endpoint_corr[,10]/165))

change_corr <- all_trt[,c(1,3,6)]
change_corr <- reshape(change_corr, idvar="Subject ID", timevar="Parameter", direction="wide")
change_corr$`Subject ID` <- NULL
change_corr <- change_corr[complete.cases(change_corr),]
change_corr$Composite_Insulin <- baseline_corr$Composite_Insulin-endpoint_corr$Composite_Insulin

#####  Add Composite Insulin
names(baseline_corr) <- c("FGLU", "FGLU_AUC", "HOMA_S", "M_VALUE","OGIS","PREDIM","QUICKI","SI", "TINSUL","TINSUL_AUC","Composite_Insulin")
names(endpoint_corr) <- c("FGLU", "FGLU_AUC", "HOMA_S", "M_VALUE","OGIS","PREDIM","QUICKI","SI", "TINSUL","TINSUL_AUC","Composite_Insulin")
names(change_corr) <- c("FGLU", "FGLU_AUC", "HOMA_S", "M_VALUE","OGIS","PREDIM","QUICKI","SI", "TINSUL","TINSUL_AUC","Composite_Insulin")

x <- data.frame(baseline_corr$M_VALUE)
y <- subset(baseline_corr, select=-c(M_VALUE))
temp_pvalue_base <- corr.test(x,y, adjust="none")
baseline_pvalue <- data.frame(temp_pvalue_base[[4]])
correlation_for_baseline <- data.frame(cor(x,y))

w <- data.frame(endpoint_corr$M_VALUE)
z <- subset(endpoint_corr, select=-c(M_VALUE))
temp_pvalue <- corr.test(w,z, adjust="none")
endpoint_pvalue <- data.frame(temp_pvalue[[4]])
correlation_for_endpoint <- data.frame(cor(w,z))

t <- data.frame(change_corr$M_VALUE)
v <- subset(change_corr, select=-c(M_VALUE))
temp_pvalue_change <- corr.test(t,v, adjust="none")
change_pvalue <- data.frame(temp_pvalue_change[[4]])
correlation_for_change <- data.frame(cor(t,v))

correlation_final <- merge(correlation_for_baseline, correlation_for_endpoint, all= TRUE)
correlation_final <- merge(correlation_final, correlation_for_change, all= TRUE)
row.names(correlation_final) <- c("Change","Baseline","Endpoint")
correlation_final_all <- correlation_final[c(2,3,1),]

pvalue_final <- merge(baseline_pvalue,endpoint_pvalue,all=TRUE)
pvalue_final <- merge(pvalue_final,change_pvalue,all=TRUE)
row.names(pvalue_final) <- c("Endpoint","Change","Baseline")
pvalue_final_all <- pvalue_final[c(3,1,2),]

write.csv(pvalue_final_all, "GWAY_ALL_pvalues.csv")
write.csv(correlation_final_all, "GWAY_ALL_correlation.csv")











