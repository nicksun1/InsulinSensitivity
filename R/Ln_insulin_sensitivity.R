suppressMessages(library(coastr))
suppressMessages(library(readr))
suppressMessages(library(dplyr))
library(tidyr)
library(psych)
library(gridExtra)

#Pull in datasets
new_params <- read.csv("/lrlhps/users/c269016/gway_IS_parameters_20190326.csv")
clamp <-import_cluwe_data(source_path="/lillyce/prd/ly2148568/h8o_us_gway/final/data/analysis/shared",data_file="clamp.sas7bdat")
efficacy <-import_cluwe_data(source_path="/lillyce/prd/ly2148568/h8o_us_gway/final/data/analysis/shared",data_file="efficacy.sas7bdat")
testmeal <- import_cluwe_data(source_path="/lillyce/prd/ly2148568/h8o_us_gway/final/data/analysis/shared",data_file="testmeal.sas7bdat")


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
si$Parameter <- "log(SI)"
si<- si[,-4]
names(si) <- c("Subject ID","Treatment", "Baseline", "Endpoint", "Change", "Parameter")
si$ln_Baseline <- log(si$Baseline)
si$ln_Endpoint <- log(si$Endpoint)
si$ln_Change <- log(si$Endpoint)-log(si$Baseline)

# Quicki
quicki_long <- long_newparams %>% select(PATIENT, VISIT, THERAPY, QUICKI)
quicki <- reshape(quicki_long, idvar="PATIENT", timevar = "VISIT", direction ="wide")
quicki$Change <- quicki$QUICKI.7-quicki$QUICKI.2
quicki$Parameter <- "log(QUICKI)"
quicki<- quicki[,-4]
names(quicki) <- c("Subject ID","Treatment", "Baseline", "Endpoint", "Change", "Parameter")
quicki$ln_Baseline <- log(quicki$Baseline)
quicki$ln_Endpoint <- log(quicki$Endpoint)
quicki$ln_Change <- log(quicki$Endpoint)-log(quicki$Baseline)

# OGIS
ogis_long <- long_newparams %>% select(PATIENT, VISIT, THERAPY, OGIS)
ogis <- reshape(ogis_long, idvar="PATIENT", timevar = "VISIT", direction ="wide")
ogis$Change <- ogis$OGIS.7-ogis$OGIS.2
ogis$Parameter <- "log(OGIS)"
ogis<- ogis[,-4]
names(ogis) <- c("Subject ID","Treatment", "Baseline", "Endpoint", "Change", "Parameter")
ogis$ln_Baseline <- log(ogis$Baseline)
ogis$ln_Endpoint <- log(ogis$Endpoint)
ogis$ln_Change <- log(ogis$Endpoint)-log(ogis$Baseline)

# logPREDIM
logpredim_long <- long_newparams%>% select(PATIENT, VISIT, THERAPY, logPREDIM)
logpredim <- reshape(logpredim_long, idvar="PATIENT", timevar = "VISIT", direction ="wide")
logpredim$Change <- logpredim$logPREDIM.7-logpredim$logPREDIM.2
logpredim$Parameter <- "log(PREDIM)"
logpredim<- logpredim[,-4]
names(logpredim) <- c("Subject ID","Treatment", "Baseline", "Endpoint", "Change", "Parameter")
logpredim$ln_Baseline <- logpredim$Baseline
logpredim$ln_Endpoint <- logpredim$Endpoint
logpredim$ln_Change <- logpredim$Endpoint-m_value$Baseline
logpredim$Parameter <- "log(PREDIM)"

### M-value
m_value_full <- clamp %>% select (PATIENT, TRT, PARACODE, VAL_B, VAL_E, VAL_C)
m_value_full <- m_value_full %>% filter(PARACODE == "M_VALUE")
m_value <- m_value_full[!duplicated(m_value_full[,c('PATIENT','TRT','PARACODE')]),]
m_value <- m_value[complete.cases(m_value),]
names(m_value) <- c("Subject ID", "Treatment", "Parameter", "Baseline","Endpoint", "Change")
m_value$ln_Baseline <- log(m_value$Baseline)
m_value$ln_Endpoint <- log(m_value$Endpoint)
m_value$ln_Change <- log(m_value$Endpoint)-log(m_value$Baseline)
m_value$Parameter <- "log(M_VALUE)"

### Log of Homa_S
homa_s_full <- efficacy %>% select (PATIENT, VISIT, TRT, PARACODE, VAL_B, VAL, VAL_C)
homa_s_full <- homa_s_full %>% filter(PARACODE == "HOMASINS")
homa_s <- homa_s_full[complete.cases(homa_s_full),]
homa_s <- homa_s %>% group_by(PATIENT) %>%
  top_n(n=1, wt= VISIT)
names(homa_s) <- c("Subject ID","Visit", "Treatment", "Parameter", "Baseline","Endpoint", "Change")
ln_homa_s <- homa_s
ln_homa_s$ln_Baseline <- log(ln_homa_s$Baseline)
ln_homa_s$ln_Endpoint <- log(ln_homa_s$Endpoint)
ln_homa_s$ln_Change <- log(ln_homa_s$Endpoint)-log(ln_homa_s$Baseline)
ln_homa_s$Parameter <- "log(Homa_S)"

### Fasting Insulin
fasting_insulin_full <- efficacy %>% select (PATIENT, VISIT, TRT, PARACODE, VAL, VAL_B, VAL_C)
fasting_insulin_full <- fasting_insulin_full %>% filter(PARACODE == "TINSUL") 
fasting_insulin <- fasting_insulin_full[complete.cases(fasting_insulin_full),]
fasting_insulin <- fasting_insulin %>% group_by(PATIENT) %>%
  top_n(n=1, wt= VISIT)
names(fasting_insulin) <- c("Subject ID","Visit", "Treatment", "Parameter", "Endpoint", "Baseline", "Change")
fasting_insulin$ln_Baseline <- log(fasting_insulin$Baseline)
fasting_insulin$ln_Endpoint <- log(fasting_insulin$Endpoint)
fasting_insulin$ln_Change <- log(fasting_insulin$Endpoint)-log(fasting_insulin$Baseline)
fasting_insulin$Parameter <- "log(TINSUL)"

### Fasting Glucose
fast_glucose_full <-  efficacy %>% select (PATIENT, VISIT, TRT, PARACODE, VAL, VAL_B, VAL_C)
fast_glucose_full <- fast_glucose_full %>% filter(PARACODE == "FGLU")
fast_glucose <- fast_glucose_full[complete.cases(fast_glucose_full),]
fast_glucose <- fast_glucose %>% group_by(PATIENT) %>%
  top_n(n=1, wt= VISIT)
names(fast_glucose) <- c("Subject ID","Visit", "Treatment", "Parameter", "Endpoint", "Baseline", "Change")
fast_glucose$ln_Baseline <- log(fast_glucose$Baseline)
fast_glucose$ln_Endpoint <- log(fast_glucose$Endpoint)
fast_glucose$ln_Change <- log(fast_glucose$Endpoint)-log(fast_glucose$Baseline)
fast_glucose$Parameter <- "log(FGLU)"

### Insulin AUC
insulin_auc_full <-  testmeal %>% select (PATIENT, VISIT, TRT, ANALYTE, PARACODE, VAL, VAL_B, VAL_E, VAL_C)
insulin_auc_full <-  insulin_auc_full %>% filter(PARACODE == "AUC") 
insulin_auc_full <-  insulin_auc_full %>% filter(ANALYTE =="TINSUL")
insulin_auc_full$PARACODE <- "log(TINSUL-AUC)"
insulin_auc_full$ANALYTE <- NULL
insulin_auc_full$VAL <- NULL
insulin_auc <- insulin_auc_full[complete.cases(insulin_auc_full),]
insulin_auc <- insulin_auc %>% group_by(PATIENT) %>%
  top_n(n=1, wt= VISIT)
names(insulin_auc) <- c("Subject ID","Visit", "Treatment", "Parameter", "Baseline", "Endpoint", "Change")
insulin_auc$ln_Baseline <- log(insulin_auc$Baseline)
insulin_auc$ln_Endpoint <- log(insulin_auc$Endpoint)
insulin_auc$ln_Change <- log(insulin_auc$Endpoint)-log(insulin_auc$Baseline)

### Glucose AUC
glucose_auc_full <-  testmeal %>% select (PATIENT, VISIT, TRT, ANALYTE, PARACODE, VAL, VAL_B, VAL_E, VAL_C)
glucose_auc_full <-  glucose_auc_full %>% filter(PARACODE == "AUC") 
glucose_auc_full <-  glucose_auc_full %>% filter(ANALYTE =="FGLU")
glucose_auc_full$PARACODE <- "log(FGLU-AUC)"
glucose_auc_full$ANALYTE <- NULL
glucose_auc_full$VAL <- NULL
glucose_auc <- glucose_auc_full[complete.cases(glucose_auc_full),]
glucose_auc <- glucose_auc %>% group_by(PATIENT) %>%
  top_n(n=1, wt= VISIT)
names(glucose_auc) <- c("Subject ID","Visit", "Treatment", "Parameter", "Baseline", "Endpoint", "Change")
glucose_auc$ln_Baseline <- log(glucose_auc$Baseline)
glucose_auc$ln_Endpoint <- log(glucose_auc$Endpoint)
glucose_auc$ln_Change <- log(glucose_auc$Endpoint)-log(glucose_auc$Baseline)

#### Merge full sets into one
lnfull_set_testmeal <- merge(insulin_auc, glucose_auc, all=TRUE)
lnfull_set_testmeal <- merge(lnfull_set_testmeal, m_value, all=TRUE)
lnfull_set_efficacy <- merge( fast_glucose, fasting_insulin, all=TRUE)
lnfull_set_efficacy <- merge(lnfull_set_efficacy, ln_homa_s, all=TRUE)
lnfull_set <- merge(lnfull_set_testmeal, lnfull_set_efficacy, all=TRUE)
lnfull_set <- merge(lnfull_set, si, all=TRUE)
lnfull_set <- merge(lnfull_set, ogis, all=TRUE)
lnfull_set <- merge(lnfull_set, logpredim, all=TRUE)
lnfull_set <- merge(lnfull_set, quicki, all=TRUE)
all_trt_ln <- lnfull_set
lnfull_set_full <- lnfull_set %>% filter (Treatment=="Rosiglitazone")
lnfull_set <- subset(lnfull_set_full, select=-c(Baseline,Endpoint,Change,Visit))


## Baseline summary for all TRT
summary_baseline_all <-  all_trt %>% dplyr::group_by(Parameter, Treatment) %>%dplyr::summarize(n=n(), Baseline_Means= mean(ln_Baseline), SD_base =sd(ln_Baseline), effect_size_base = mean(ln_Baseline)/sd(ln_Baseline))
write.csv(summary_baseline_all, "Log_GWAY_insulin_sensitivity_baseline_alltrt.csv")
summary_baseline_all_comb <-  all_trt %>% dplyr::group_by(Parameter) %>%dplyr::summarize(n=n(), Baseline_Means= mean(ln_Baseline), SD_base =sd(ln_Baseline), effect_size_base = mean(ln_Baseline)/sd(ln_Baseline))
write.csv(summary_baseline_all_comb, "Log_GWAY_insulin_sensitivity_baseline_alltrtcomb.csv")

baseline_compins <- all_trt[,c(1,3:4,7)]
baseline_compins <- reshape(baseline_compins, idvar="Subject ID", timevar="Parameter", direction="wide")
baseline_compins$`Subject ID` <- NULL
baseline_compins$Composite_Insulin <- log(10000/(sqrt(baseline_compins[,4]*baseline_compins[,20]*baseline_compins[,2]/165*baseline_compins[,18]/165)))
baseline_compins <- baseline_compins[complete.cases(baseline_compins[,21]),]

summary_baseline_comp <-  baseline_compins %>%dplyr::summarize(n=101, Baseline_Means= mean(Composite_Insulin), SD_base =sd(Composite_Insulin), effect_size_base = mean(Composite_Insulin)/sd(Composite_Insulin))

## Baseline summary for Ros only
summary_baseline <-  lnfull_set %>% dplyr::group_by(Parameter) %>%dplyr::summarize(n=n(), Baseline_Means= mean(ln_Baseline), SD_base =sd(ln_Baseline), effect_size_base = mean(ln_Baseline)/sd(ln_Baseline))
write.csv(summary_baseline, "Log_GWAY_insulin_sensitivity_summary.csv")

summary_change <-  lnfull_set %>% dplyr::group_by(Parameter) %>%dplyr::summarize(n=n(), Change_Means= mean(ln_Change), SD_change =sd(ln_Change), effect_size_change = mean(ln_Change)/sd(ln_Change))
write.csv(summary_change, "Log_GWAY_insulin_sensitivity_summary_change.csv")

baseline_compins <- lnfull_set_full[,c(1,3:4,7)]
baseline_compins <- reshape(baseline_compins, idvar="Subject ID", timevar="Parameter", direction="wide")
baseline_compins$`Subject ID` <- NULL
baseline_compins$Composite_Insulin <- log(10000/(sqrt(baseline_compins[,4]*baseline_compins[,20]*baseline_compins[,2]/165*baseline_compins[,18]/165)))
baseline_compins <- baseline_compins[complete.cases(baseline_compins[,21]),]

endpoint_compins <- lnfull_set_full[,c(1,3,5,8)]
endpoint_compins <- reshape(endpoint_compins, idvar="Subject ID", timevar="Parameter", direction="wide")
endpoint_compins$`Subject ID` <- NULL
endpoint_compins$Composite_Insulin <- log(10000/(sqrt(endpoint_compins[,4]*endpoint_compins[,20]*endpoint_compins[,2]/165*endpoint_compins[,18]/165)))
endpoint_compins <- endpoint_compins[complete.cases(endpoint_compins[,21]),]

change_compins <- lnfull_set[,c(1,3,6)]
change_compins <- reshape(change_compins, idvar="Subject ID", timevar="Parameter", direction="wide")
change_compins$`Subject ID` <- NULL
change_compins <- change_compins[complete.cases(change_compins[,9]),]
change_compins$Composite_Insulin <- endpoint_compins$Composite_Insulin-baseline_compins$Composite_Insulin

summary_change_comp <-  change_compins %>%dplyr::summarize(n=33, Change_Means= mean(Composite_Insulin), SD_change =sd(Composite_Insulin), effect_size_change = mean(Composite_Insulin)/sd(Composite_Insulin))




library(plyr)
## Correlation Matrices /Add Composite Insulin
# Create data to do correlation on
baseline_corr <- lnfull_set_full[,c(1,3:4,7)]
baseline_corr <- reshape(baseline_corr, idvar="Subject ID", timevar="Parameter", direction="wide")
baseline_corr$`Subject ID` <- NULL
baseline_corr <- baseline_corr[complete.cases(baseline_corr[,c(1:11,17:20)]),]
baseline_corr$Composite_Insulin <- log(10000/(sqrt(baseline_corr[,4]*baseline_corr[,20]*baseline_corr[,2]/165*baseline_corr[,18]/165)))
baseline_corr <- baseline_corr[,c(1,3,5,7,9,11,13,15,17,19,21)]

endpoint_corr <- lnfull_set_full[,c(1,3,5,8)]
endpoint_corr <- reshape(endpoint_corr, idvar="Subject ID", timevar="Parameter", direction="wide")
endpoint_corr$`Subject ID` <- NULL
endpoint_corr <- endpoint_corr[complete.cases(endpoint_corr[,1:11]),]
#endpoint_corr$Composite_Insulin <- 10000/(sqrt(endpoint_corr$Baseline.log(FGLU)*endpoint_corr$Baseline.log(TINSUL)*endpoint_corr$Baseline.log(FGLU_AUC)/165*endpoint_corr$Baseline.log(TINSUL_AUC)/165))
endpoint_corr$Composite_Insulin <- log(10000/(sqrt(endpoint_corr[,4]*endpoint_corr[,20]*endpoint_corr[,2]/165*endpoint_corr[,18]/165)))
endpoint_corr <- endpoint_corr[,c(1,3,5,7,9,11,13,15,17,19,21)]

change_corr <- lnfull_set[,c(1,3,6)]
change_corr <- reshape(change_corr, idvar="Subject ID", timevar="Parameter", direction="wide")
change_corr$`Subject ID` <- NULL
change_corr <- change_corr[complete.cases(change_corr[]),]
change_corr$Composite_Insulin <- baseline_corr$Composite_Insulin-endpoint_corr$Composite_Insulin


names(baseline_corr) <- c("log_of_FGLU_AUC", "log_of_FGLU", "log_of_HOMA_S", "log_of_M_VALUE","log_of_OGIS","log_of_PREDIM","log_of_QUICKI","log_of_SI", "log_of_TINSUL_AUC","log_of_TINSUL", "log_of_Composite_Insulin")
names(endpoint_corr) <- c("log_of_FGLU_AUC", "log_of_FGLU", "log_of_HOMA_S", "log_of_M_VALUE","log_of_OGIS","log_of_PREDIM","log_of_QUICKI","log_of_SI", "log_of_TINSUL_AUC","log_of_TINSUL", "log_of_Composite_Insulin")
names(change_corr) <- c("log_of_FGLU_AUC", "log_of_FGLU", "log_of_HOMA_S", "log_of_M_VALUE","log_of_OGIS","log_of_PREDIM","log_of_QUICKI","log_of_SI", "log_of_TINSUL_AUC","log_of_TINSUL", "log_of_Composite_Insulin")


###################### Make Correlation Matrix
x <- data.frame(baseline_corr$log_of_M_VALUE)
y <- subset(baseline_corr, select=-c(log_of_M_VALUE))
temp_pvalue_base <- corr.test(x,y, adjust="none")
baseline_pvalue <- data.frame(temp_pvalue_base[[4]])
correlation_for_baseline <- data.frame(cor(x,y))

w <- data.frame(endpoint_corr$log_of_M_VALUE)
z <- subset(endpoint_corr, select=-c(log_of_M_VALUE))
temp_pvalue <- corr.test(w,z, adjust="none")
endpoint_pvalue <- data.frame(temp_pvalue[[4]])
correlation_for_endpoint <- data.frame(cor(w,z))

t <- data.frame(change_corr$log_of_M_VALUE)
v <- subset(change_corr, select=-c(log_of_M_VALUE))
temp_pvalue_change <- corr.test(t,v, adjust="none")
change_pvalue <- data.frame(temp_pvalue_change[[4]])
correlation_for_change <- data.frame(cor(t,v))

correlation_final <- merge(correlation_for_baseline, correlation_for_endpoint, all= TRUE)
correlation_final <- merge(correlation_final, correlation_for_change, all= TRUE)
row.names(correlation_final) <- c("ln_Endpoint","ln_Baseline","ln_Change")
library(ggplot2)


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
baseline_corr1 <- baseline_corr[,c(4,1,2,3,5,6,7,8,9,10,11)]
for(i in 2:length(baseline_corr1)){
  b1= ggplot(baseline_corr1, aes(x=log_of_M_VALUE, y=baseline_corr1[,i])) + geom_point()+
    geom_smooth(method='lm') + geom_text(data = labels, aes(x = max(baseline_corr1$log_of_M_VALUE)-2, y = max(baseline_corr1[,i]), label = lm_eqn_base(i)), parse = TRUE) +
    ggtitle(paste("Baseline Comparison for log_of_M-value v",colnames(baseline_corr1)[i])) + ylab(colnames(baseline_corr1)[i])
  b[[i]] <- b1
}
pdf("Log_Baseline_plots.pdf")
for(i in 2:length(baseline_corr1)){
  print(b[[i]])
}
dev.off()


#################### Endpoint graphs

e<- list()
endpoint_corr1 <- endpoint_corr[,c(4,1,2,3,5,6,7,8,9,10,11)]
for(i in 2:length(endpoint_corr1)){
  e1= ggplot(endpoint_corr1, aes(x=log_of_M_VALUE, y=endpoint_corr1[,i])) + geom_point()+
    geom_smooth(method='lm') + geom_text(data = labels, aes(x = max(endpoint_corr1$log_of_M_VALUE)-2, y = max(endpoint_corr1[,i]), label = lm_eqn_end(i)), parse = TRUE) +
    ggtitle(paste("Baseline Comparison for log_of_M-value v",colnames(endpoint_corr1)[i])) + ylab(colnames(endpoint_corr1)[i])
  e[[i]] <- e1
}
pdf("Log_Endpoint_plots.pdf")
for(i in 2:length(endpoint_corr1)){
  print(e[[i]])
}
dev.off()


#################### Change graphs

c<- list()
change_corr1 <- change_corr[,c(4,1,2,3,5,6,7,8,9,10,11)]
for(i in 2:length(change_corr1)){
  c1= ggplot(change_corr1, aes(x=log_of_M_VALUE, y=change_corr1[,i])) + geom_point()+
    geom_smooth(method='lm') + geom_text(data = labels, aes(x = max(change_corr1$log_of_M_VALUE)-2, y = max(change_corr1[,i]), label = lm_eqn_change(i)), parse = TRUE) +
    ggtitle(paste("Baseline Comparison for log_of_M-value v",colnames(change_corr1)[i])) + ylab(colnames(change_corr1)[i])
  c[[i]] <- c1
}
pdf("Log_Change_plots.pdf")
for(i in 2:length(change_corr1)){
  print(c[[i]])
}
dev.off()


### Summary of Variables



summary <- full_set %>% group_by(Parameter) %>% summarise( n=n(),mean_base= mean(Baseline), sd_base =sd(Baseline), effect_size_base = mean(Baseline)/sd(Baseline),
                                                           mean_end= mean(Endpoint), sd_end =sd(Endpoint),effect_size_end = mean(Endpoint)/sd(Endpoint),
                                                           mean_change= mean(Change), sd_change =sd(Change),effect_size_change = mean(Change)/sd(Change))
summary_by_trt_ln <- lnfull_set %>% group_by(Parameter,Treatment) %>% summarise( n=n(),mean_base= mean(Baseline), sd_base =sd(Baseline), effect_size_base = mean(Baseline)/sd(Baseline),
                                                                                 mean_end= mean(Endpoint), sd_end =sd(Endpoint),effect_size_end = mean(Endpoint)/sd(Endpoint),
                                                                                 mean_change= mean(Change), sd_change =sd(Change),effect_size_change = mean(Change)/sd(Change))
colnames(summary_by_trt_ln)[1]<-"ln(Parameter)"
write.csv(summary, "GWAY_insulin_sensitivity_summary.csv")
write.csv(summary_by_trt_ln, "Logof_GWAY_insulin_sensitivity_summary_by_treatment.csv")





