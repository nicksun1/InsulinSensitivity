# Insulin Sensitivity
#
# For All Treatments
#
# Calculates and outputs correlation and corresponding p-value between M-value and various indices.
# Creates and outputs plots to show each relationship over time and their variability. 
# Above done for baseline, endpoint, and change values.
# 
# Data was created using rnorm(nrow(), mean=1, sd=.3) on real data and filtered down to variables of interest for this particular project.
#


source("H:/R/InsulinSensitivity/R/Indices_Data.R")


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

write.csv(pvalue_final_all, "ALL_pvalues.csv")
write.csv(correlation_final_all, "ALL_correlation.csv")



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


