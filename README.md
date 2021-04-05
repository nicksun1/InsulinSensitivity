# InsulinSensitivity

Objective: 
  -To evaluate the correlation between commmonly used insulin sensitivity indices with the M-value from euglycemic clamp test
  -Understand the magnitude of change and the variability in insulin sensitivity indices for patients with Type-2 Diabetes
  -Evaluate the time trend of insulin sensitivity change over the course of treatment

Refer to Insulin Sensitivity Reference.png for summary of measures and abbreviations.
- OGTT -> oral glucose tolerance test; MMTT -> mixed meal tolerance test.
- M-value, used to quantify whole-body insulin sensitivity and considered the gold standard, is the measurement of insulin stimulated glucose disposal during the euglycemic clamp. 

Data recreated through rnorm(nrow(),mean=1,sd=0.2) results include:
  -Calculation of Mastuda Index (composite insulin)
  -Correlation and corresponding p-value between M-value and various insulin sensitivity indices
  -Regression graphs of M-value against various insulin sensitivity indices with significance and variability factors (p-value and R^2)
  -Summary of analysis data with N, Mean, SD, and effect size for by treatment and combined.
  
Indices_Data.R and Indices_Ln_Data.R are the starting points which manipulate data and arrange data into workable format.