# Insulin Sensitivity
#
# Log Data (ln()) for FGLU, FGLU_AUC, HOMA_S, M_VALUE, OGIS, PREDIM, QUICKI, SI, TINSUL, TINSUL_AUC, Composite_Insulin
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
lnfull_set <- subset(lnfull_set, select=-c(Visit))
all_trt_ln <- lnfull_set
lnfull_set <- lnfull_set %>% filter (Treatment=="Rosiglitazone")




# Calculate Composite Insulin 
baseline_compins <- lnfull_set[,c(1,3:4,7)]
baseline_compins <- reshape(baseline_compins, idvar="Subject ID", timevar="Parameter", direction="wide")
baseline_compins$`Subject ID` <- NULL
baseline_compins$Composite_Insulin <- log(10000/(sqrt(baseline_compins[,4]*baseline_compins[,20]*baseline_compins[,2]/165*baseline_compins[,18]/165)))
baseline_compins_ros <- baseline_compins[complete.cases(baseline_compins[,21]),]

endpoint_compins <- lnfull_set[,c(1,3,5,8)]
endpoint_compins <- reshape(endpoint_compins, idvar="Subject ID", timevar="Parameter", direction="wide")
endpoint_compins$`Subject ID` <- NULL
endpoint_compins$Composite_Insulin <- log(10000/(sqrt(endpoint_compins[,4]*endpoint_compins[,20]*endpoint_compins[,2]/165*endpoint_compins[,18]/165)))
endpoint_compins_ros <- endpoint_compins[complete.cases(endpoint_compins[,21]),]

change_compins <- lnfull_set[,c(1,3,6)]
change_compins <- reshape(change_compins, idvar="Subject ID", timevar="Parameter", direction="wide")
change_compins$`Subject ID` <- NULL
change_compins <- change_compins[complete.cases(change_compins[,9]),]
change_compins$Composite_Insulin <- endpoint_compins_ros$Composite_Insulin-baseline_compins_ros$Composite_Insulin
change_compins_ros <- change_compins


# All Treatment Composite Insulin 
baseline_compins <- all_trt_ln[,c(1,3:4,7)]
baseline_compins <- reshape(baseline_compins, idvar="Subject ID", timevar="Parameter", direction="wide")
baseline_compins$`Subject ID` <- NULL
baseline_compins$Composite_Insulin <- log(10000/(sqrt(baseline_compins[,4]*baseline_compins[,20]*baseline_compins[,2]/165*baseline_compins[,18]/165)))
baseline_compins_all <- baseline_compins[complete.cases(baseline_compins[,21]),]

endpoint_compins <- all_trt_ln[,c(1,3,5,8)]
endpoint_compins <- reshape(endpoint_compins, idvar="Subject ID", timevar="Parameter", direction="wide")
endpoint_compins$`Subject ID` <- NULL
endpoint_compins$Composite_Insulin <- log(10000/(sqrt(endpoint_compins[,4]*endpoint_compins[,20]*endpoint_compins[,2]/165*endpoint_compins[,18]/165)))
endpoint_compins_all <- endpoint_compins[complete.cases(endpoint_compins[,21]),]

change_compins <- all_trt_ln[,c(1,3,6)]
change_compins <- reshape(change_compins, idvar="Subject ID", timevar="Parameter", direction="wide")
change_compins$`Subject ID` <- NULL
change_compins <- change_compins[complete.cases(change_compins[,c(2,9)]),]
change_compins$Composite_Insulin <- endpoint_compins_all$Composite_Insulin-baseline_compins_all$Composite_Insulin
change_compins_all <- change_compins



