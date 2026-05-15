#-------cox regression----------------------------------------------------
independent_vars<- c("sar_base_diagnose","ALM_SMI", "UWP_i0","HGS")
covariates <- c("sex", "age_i0", "bmi_i0","Mixed","Asian","Black","Chinese","others")
data$UWP_i0 <- as.factor(data$UWP_i0)
 
for (ind_var in independent_vars) {
  formula <- as.formula(paste("Surv(time_diff, status) ~", ind_var, "+", paste(covariates, collapse = "+")))
  cox_model <- coxph(formula, data = data)  
  print(summary(cox_model))
}

results_matrix <- matrix(nrow = length(independent_vars), ncol = 9)

for (i in 1:length(independent_vars)) {
  formula <- as.formula(paste("Surv(time_diff, status) ~", independent_vars[i], "+", paste(covariates, collapse = "+")))
  print(formula) 
  cox_model <- coxph(formula, data = data) 
  model_summary <- summary(cox_model)
  col <- c(colnames(model_summary$coefficients), colnames(model_summary$conf.int))
  colnames(results_matrix) <- col
  results_matrix[i,1:5] <- model_summary$coefficients[1,]
  results_matrix[i,6:9] <- model_summary$conf.int[1,]
}
rownames(results_matrix) <- independent_vars
write.csv(results_matrix,"res_cox_ost_i0.csv")

#-------logistic regression----------------------------------------------------

result_df <- data.frame(
  Variable = character(),      
  Estimate = numeric(),        
  StdError = numeric(),       
  ZValue = numeric(),         
  PValue = numeric(),          
  OR = numeric(),             
  OR_lower = numeric(),        
  OR_upper = numeric(),        
  stringsAsFactors = FALSE     
)

covariates <- c("time_diff", "sex", "age_i0", "bmi_i0", "Mixed", 
                "Asian", "Black", "Chinese", "others")

for (var in predictors) {
  formula <- as.formula(paste("ProbSAR_i2~", 
                              paste(c(var, covariates), collapse = " + ")))
  

  model <- glm(formula, data = data, family = binomial)
  model_summary <- summary(model)
  coef_summary <- model_summary$coefficients[var, ]
  OR <- exp(coef_summary["Estimate"]) 
  conf_int <- exp(confint(model))[var, ]  

  result_df <- rbind(result_df, data.frame(
    Variable = var,
    Estimate = coef_summary["Estimate"],
    StdError = coef_summary["Std. Error"],
    ZValue = coef_summary["z value"],
    PValue = coef_summary["Pr(>|z|)"],
    OR = OR,
    OR_lower = conf_int[1],
    OR_upper = conf_int[2]
  ))
}

print(result_df)
write.csv(result_df,"res_logistic_ProbSAR.csv")


#---- interaction--------------------------------------

library(survival)
library(dplyr)


data$ALM_z <- scale(data$ALM0_height)


cox_model <- coxph(Surv(time_diff, status) ~ 
                     sex + age_i0 + bmi_i0 + Mixed + Asian + Black + Chinese + others + 
                     ALM_z +
                     sex:bmi_i0+
                     sex:age_i0+
                     sex:Mixed +
                     sex:Asian+
                     sex:Black+
                     sex:Chinese+
                     sex:others+
                     sex:ALM_z,
                   data = data)

cox_summary <- summary(cox_model)

result <- data.frame(
  coef = cox_summary$coefficients[, "coef"],
  exp_coef = cox_summary$coefficients[, "exp(coef)"],
  se_coef = cox_summary$coefficients[, "se(coef)"],
  z = cox_summary$coefficients[, "z"],
  p_value = cox_summary$coefficients[, "Pr(>|z|)"],
  lower_95 = cox_summary$conf.int[, "lower .95"],
  upper_95 = cox_summary$conf.int[, "upper .95"]
)

#-----additive interaction------------------------------------

library(interactionR)
cleaned_data <- data %>%
  filter(HGS_classify3 != 3, UWP_i0 != 3)
updated_data <- data %>%
  mutate(
    HGS_classify3 = as.numeric(HGS_classify3),   
    UWP_i0 = as.numeric(UWP_i0),                 
    HGS_classify3 = if_else(HGS_classify3 == 3, 1, HGS_classify3),
    UWP_i0 = if_else(UWP_i0 == 3, 1, UWP_i0)
  )
updated_data$ALM_classify3 <- as.factor(updated_data$ALM_classify3)
updated_data$HGS_classify3<- as.factor(updated_data$HGS_classify3)
updated_data$UWP_i0 <- as.factor(updated_data$UWP_i0)
fit1 <- coxph(Surv(time_diff, status) ~ HGS_classify3*UWP_i0 + 
                sex + age_i0 + bmi_i0 + Mixed + Asian + Black + Chinese + others, 
              data = updated_data)
fit1
summary(fit1)
library(car)
 
out<- interactionR(fit1, 
                   exposure_names =c("HGS_classify32", "UWP_i02"), 
                   ci.type ="delta", ci.level = 0.95,em = F, recode = F)

out$dframe
interactionR_table(out) 


#----------------PAF---------------------------------

data <- data %>%
  mutate(HGS_inverted3 = case_when(
    HGS_classify3 == 1 ~ 3,
    HGS_classify3 == 2 ~ 2,
    HGS_classify3 == 3 ~ 1,
    TRUE ~ NA_real_   
  ))
table(data$HGS_classify3)
table(data$HGS_inverted3)
 
data$HGS_inverted3 <- as.factor(data$HGS_inverted3) 
cox_model <- coxph(Surv(time_diff, status) ~ HGS_inverted3 + 
                     sex + age_i0 + bmi_i0 + Mixed + Asian + Black + Chinese + others, 
                   data = data)
 
summary(cox_model)
cox_summary <- summary(cox_model)
 
coef_HGS3 <-  cox_model$coef[which(names(cox_model$coef) == "HGS_inverted33")]
se_HGS3 <- summary(cox_model)$coefficients[which(rownames(summary(cox_model)$coefficients) == "HGS_inverted33"), "se(coef)"]
HR_HGS3 <- exp(cox_model$coef[which(names(cox_model$coef) == "HGS_inverted33")])  # 获取 HGS_inverted3 = 2 的风险比
 
z_score <- qnorm(0.975)  
lower_HR_HGS3 <- exp(coef_HGS3 - z_score * se_HGS3)
upper_HR_HGS3 <- exp(coef_HGS3 + z_score * se_HGS3)
 
P_e <- mean(data$HGS_inverted3 == 3, na.rm = TRUE)  
PAF <- (P_e * (HR_HGS3 - 1)) / (1 + P_e * (HR_HGS3 - 1))
 
PAF_lower <- (P_e * (lower_HR_HGS3 - 1)) / (1 + P_e * (lower_HR_HGS3 - 1))
PAF_upper <- (P_e * (upper_HR_HGS3 - 1)) / (1 + P_e * (upper_HR_HGS3 - 1))

cat("PAF: ", PAF, "\n")
cat("95% CI for PAF: ", PAF_lower, " to ", PAF_upper, "\n")


#---------------- SAR-OP mediation---------------------------------
library(data.table)
library(survival)
library(tidyr)
library(foreach)
library(doParallel)
registerDoParallel(cores = 50)


boot_function = function(boot_sample){
  
  data = data0[boot_sample,]
  
  lm_formula <- as.formula(paste(mediator, "~", exposure, " + sex + age_i0 + bmi_i0 + Mixed + Asian + Black + Chinese + others"))
  alpha.temp <- coefficients(lm(lm_formula, data))[2]  
  
  cox_formula1 <- as.formula(paste("Surv(time_diff, status == 1)  ~", mediator, "+", exposure, " + sex + age_i0 + bmi_i0 + Mixed + Asian + Black + Chinese + others"))
  coeff = coefficients(coxph(cox_formula1, data))
  beta.temp <- coeff[1]
  c_prime.temp <- coeff[2]
  
  cox_formula3 <- as.formula(paste("Surv(time_diff, status == 1)  ~", exposure, " + sex + age_i0 + bmi_i0 + Mixed + Asian + Black + Chinese + others"))
  c.temp <- coefficients(coxph(cox_formula3, data))[1]
  
  IE1 <- alpha.temp * beta.temp
  DE <- c_prime.temp
  TOT <- c.temp
  
  results <- c(IE1,  DE, TOT)
  return(results)
  
}

exposure = 'ALM0_SMI'

result = numeric()
for(mediator in proteins){
  est_obs = boot_function(1:n)
  
  est_boot = foreach(j= 1:1000, .packages='foreach')%dopar%{
    set.seed(j)
    boot_sample = sample(1:n, replace = T)
    res = boot_function(boot_sample)[1]
    return(res)
  }
  
  
  est_boot =  sort(est_boot, decreasing = F)
  ci = est_boot[c(25, 975)]
  se = sd(est_boot)
  p.approx = pnorm(abs(est_obs[1]/se), lower.tail = F) * 2 
  
  result = rbind(result, c(mediator,est_obs,se,ci,p.approx))  
}

#---------------- lifestyle-disease mediation---------------------------------

boot_function = function(boot_sample){
  
  data = data0[boot_sample,]
  
  formula <- as.formula(paste(mediator, "~", exposure, " + sex + age_i0 + bmi_i0 + Mixed + Asian + Black + Chinese + others"))
  alpha.temp <- coefficients(lm(formula, data))[2]  
  
  formula <- as.formula(paste(outcome, " ~", mediator, "+", exposure, " + sex + age_i0 + bmi_i0 + Mixed + Asian + Black + Chinese + others"))
  coeff = coefficients(glm(formula, data, family=binomial(link = "logit")))
  beta.temp <- coeff[2]
  c_prime.temp <- coeff[3]
  
  formula <- as.formula(paste(outcome, " ~", exposure, " + sex + age_i0 + bmi_i0 + Mixed + Asian + Black + Chinese + others"))
  c.temp <- coefficients(glm(formula, data, family=binomial(link = "logit")))[2]
  
  IE1 <- alpha.temp * beta.temp
  DE <- c_prime.temp
  TOT <- c.temp
  
  results <- c(IE1,  DE, TOT)
  return(results)
  
}

exposure = 'sleep'
outcome = 'OP'

result = numeric()
for(mediator in proteins){
  est_obs = boot_function(1:n)
  
  est_boot = foreach(j= 1:1000, .packages='foreach')%dopar%{
    set.seed(j)
    boot_sample = sample(1:n, replace = T)
    res = boot_function(boot_sample)[1]
    return(res)
  }
  
  
  est_boot =  sort(est_boot, decreasing = F)
  ci = est_boot[c(25, 975)]
  se = sd(est_boot)
  p.approx = pnorm(abs(est_obs[1]/se), lower.tail = F) * 2 
  
  result = rbind(result, c(mediator,est_obs,se,ci,p.approx))  
}

















