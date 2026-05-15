library(data.table)
library(dplyr)


#------------- LD clumping ---------------------

# phen = 'SAR'
phen = 'OP'
dir = paste0("/fastGWA_result/",phen)
setwd(dir)
gwas = data.table:: fread(paste0(phen,".fastGWA"))


library(ieugwasr)
bfile = "/data2/xsq/1kg.v3/EUR"  # path of reference file for LD computation

gwas = gwas[gwas$P<0.05,]
gwas = dplyr::rename(gwas, rsid=SNP, pval=P)
gwas$SNP = paste0(gwas$rsid, '_', gwas$A1,gwas$A2)

gwas_clump = ld_clump(dat = gwas, 
                    clump_kb = 5000,
                    clump_r2 = 0.01,
                    clump_p = 1,
                    plink_bin = genetics.binaRies::get_plink_binary(),
                    bfile = bfile,
                    pop = "EUR"
)


write.table(gwas_clump$SNP, file="clump_snp.txt",row.names = F,col.names = F, quote = F)
write.table(gwas_clump, file="gwas_clump.txt",row.names = F, quote = F)



#------------- PRS ---------------------

prs_sum=numeric()
for(chr in 1:22){
  
  filename = paste0("prs/",phen,"/chr",chr,".sscore")
  
  if(!file.exists(filename)){next}
  prs <- read.delim(filename)
  
  prs_sum = cbind(prs_sum, prs$SCORE1_SUM)
}

prs_sum = cbind(FID=prs$X.FID, IID=prs$IID, PRS=rowSums(prs_sum))


#-------------  mediation ---------------------
library(data.table)
library(survival)
library(tidyr)
library(foreach)
library(doParallel)
registerDoParallel(cores = 50)


exposure = 'PRS_OP'
outcome = 'OP'

boot_function = function(boot_sample,mediator,exposure,outcome){
  
  data = data0[boot_sample,]
  
  formula <- as.formula(paste(mediator, "~", exposure, " + sex + age_i0 + bmi_i0 + Mixed + Asian + Black + Chinese + others"))
  alpha.temp <- coefficients(lm(formula, data))[2]  
  
  formula <- as.formula(paste(outcome," ~", mediator, "+", exposure, " + sex + age_i0 + bmi_i0 + Mixed + Asian + Black + Chinese + others"))
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

result = numeric()
for(mediator in proteins){
  est_obs = boot_function(1:n,mediator,exposure,outcome)
  
  est_boot = foreach(j= 1:1000, .packages='foreach')%dopar%{
    set.seed(j)
    boot_sample = sample(1:n, replace = T)
    res = boot_function(boot_sample,mediator,exposure,outcome)[1]
    return(res)
  }
  
  
  est_boot =  sort(est_boot, decreasing = F)
  ci = est_boot[c(25, 975)]
  se = sd(est_boot)
  p.approx = pnorm(abs(est_obs[1]/se), lower.tail = F) * 2 
  
  result = rbind(result, c(mediator,est_obs,se,ci,p.approx))  
}






















