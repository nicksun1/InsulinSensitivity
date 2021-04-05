# Insulin Sensitivity
#
# Non-logged Data for FGLU, FGLU_AUC, HOMA_S, M_VALUE, OGIS, PREDIM, QUICKI, SI, TINSUL, TINSUL_AUC, Composite_Insulin
#
# Calculates Mastuda Index (composite insulin)
# Takes in larger study dataset, manipulates and prepares data for analysis.
#
# Data has been creating using rnorm(nrow(), mean=1, sd=.3) on real data and filtered down to variables of interest for this particular project.
#


suppressMessages(library(readr))
#library(plyr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(psych)
library(gridExtra)



#####################

new_params <- read.csv("IS_parameters.csv")
clamp <- read.csv("clamp_new.csv")
efficacy <- read.csv("efficacy_new.csv")
testmeal <- read.csv("testmeal_new.csv")


# Filter out non complete data from csv (SI, Quicki, OGIS, logPREDIM)
non_dup <-  data.frame(new_params[!duplicated(new_params$PATIENT),])
sort_dup <- new_params[order(new_params$PATIENT, -new_params$VISIT),]
dup <- data.frame(sort_dup[duplicated(new_params$PATIENT),])
diff <- setdiff(non_dup, dup)
final_new_params <- setdiff(new_params,diff)
final_new_params$THERAPY <- as.character(final_new_params$THERAPY)
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
full_set <- full_set %>% filter (Treatment=="TRT_A")
full_set_complete <- full_set[complete.cases(full_set),]


# Calculate Composite Insulin 
baseline_compins <- full_set[,c(1,3:4)]
baseline_compins <- reshape(baseline_compins, idvar="Subject ID", timevar="Parameter", direction="wide")
baseline_compins$`Subject ID` <- NULL
baseline_compins$Composite_Insulin <- 10000/(sqrt(baseline_compins[,1]*baseline_compins[,9]*baseline_compins[,2]/165*baseline_compins[,10]/165))
baseline_compins_trt_a <- baseline_compins[complete.cases(baseline_compins[,11]),]

endpoint_compins <- full_set[,c(1,3,5)]
endpoint_compins <- reshape(endpoint_compins, idvar="Subject ID", timevar="Parameter", direction="wide")
endpoint_compins$`Subject ID` <- NULL
endpoint_compins$Composite_Insulin <- 10000/(sqrt(endpoint_compins[,1]*endpoint_compins[,9]*endpoint_compins[,2]/165*endpoint_compins[,10]/165))
endpoint_compins_trt_a <- endpoint_compins[complete.cases(endpoint_compins[,11]),]

change_compins <- full_set[,c(1,3,6)]
change_compins <- reshape(change_compins, idvar="Subject ID", timevar="Parameter", direction="wide")
change_compins$`Subject ID` <- NULL
change_compins <- change_compins[complete.cases(change_compins[,c(1,2,9,10)]),]
change_compins$Composite_Insulin <- endpoint_compins_trt_a$Composite_Insulin-baseline_compins_trt_a$Composite_Insulin
change_compins_trt_a <- change_compins

# All Treatment Composite Insulin
baseline_compins <- all_trt[,c(1,3:4)]
baseline_compins <- reshape(baseline_compins, idvar="Subject ID", timevar="Parameter", direction="wide")
baseline_compins$`Subject ID` <- NULL
baseline_compins$Composite_Insulin <- 10000/(sqrt(baseline_compins[,1]*baseline_compins[,9]*baseline_compins[,2]/165*baseline_compins[,10]/165))
baseline_compins_all <- baseline_compins[complete.cases(baseline_compins[,11]),]

endpoint_compins <- all_trt[,c(1,3,5)]
endpoint_compins <- reshape(endpoint_compins, idvar="Subject ID", timevar="Parameter", direction="wide")
endpoint_compins$`Subject ID` <- NULL
endpoint_compins$Composite_Insulin <- 10000/(sqrt(endpoint_compins[,1]*endpoint_compins[,9]*endpoint_compins[,2]/165*endpoint_compins[,10]/165))
endpoint_compins_all <- endpoint_compins[complete.cases(endpoint_compins[,11]),]

change_compins <- all_trt[,c(1,3,6)]
change_compins <- reshape(change_compins, idvar="Subject ID", timevar="Parameter", direction="wide")
change_compins$`Subject ID` <- NULL
change_compins_all <- change_compins[complete.cases(change_compins[,c(1,2,9,10)]),]
change_compins_all$Composite_Insulin <- endpoint_compins_all$Composite_Insulin-baseline_compins_all$Composite_Insulin












