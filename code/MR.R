library(TwoSampleMR)
library(data.table)
library(dplyr)
library(tidyr)


direction = 'sar_op' 
# direction = 'op_sar'

#------------------------------- read data ----------------------------------------------

load(file='data/mr_processed/op_sar/expo_sar.RData')
expo_sar <- expo_sar %>% filter(mr_keep.exposure=="TRUE")

load(file='data/mr_processed/op_sar/expo_op.RData')
expo_op <- expo_op %>% filter(mr_keep.exposure=="TRUE")

load(file='data/mr_processed/op_sar/outcome_op.RData')
outcome_op <- outcome_op %>% filter(mr_keep.outcome=="TRUE")

load(file='data/mr_processed/op_sar/outcome_sar.RData')
outcome_sar <- outcome_sar %>% filter(mr_keep.outcome=="TRUE")

if(direction == 'sar_op'){
  dat <- harmonise_data( 
    exposure_dat = expo_sar, 
    outcome_dat = outcome_op, 
    action = 1)
}else{
  dat <- harmonise_data( 
    exposure_dat = expo_op, 
    outcome_dat = outcome_sar, 
    action = 1)
}



#----------------------------- MR analysis -----------------------------------------------

mr_res <- mr(dat, method_list = c("mr_ivw"))
mr_res$sig = mr_res$pval < 0.05/9
median_res <- mr(dat, method_list=c("mr_weighted_median"))   


#----------------------------- Steiger filtering  -------------------------------------

N_UWP = 459915
N_WFFM = 454850
N_RLFFM = 454835
N_LLFFM = 454805
N_RAFFM = 454753
N_LAFFM = 454672
N_TFFM = 454508
N_RHGS = 461089
N_LHGS = 461026
N_LS = 33297 

mr_res <- mr_res %>% 
  unite("test",c("outcome","exposure"),sep="_",remove=F)
sig_test <- mr_res %>% filter(sig=TRUE)


dat <- dat %>% 
  unite("test",c("outcome","exposure"),sep="_",remove=F)


dat_steiger <- dat %>% 
  filter(test %in% sig_test$test) %>% 
  mutate(n_exp = case_when(exposure == "UWP" ~ N_UWP,
                           exposure == "WFFM" ~ N_WFFM,
                           exposure == "RLFFM" ~ N_RLFFM,
                           exposure == "LLFFM" ~ N_LLFFM,
                           exposure == "RAFFM" ~ N_RAFFM,
                           exposure == "LAFFM" ~ N_LAFFM,
                           exposure == "TFFM" ~ N_TFFM, 
                           exposure == "RHGS" ~ N_RHGS,
                           exposure == "LHGS" ~ N_LHGS,
                           exposure == "LS" ~ N_LS),
         n_out = case_when(outcome == "UWP" ~ N_UWP,
                           outcome == "WFFM" ~ N_WFFM,
                           outcome == "RLFFM" ~ N_RLFFM,
                           outcome == "LLFFM" ~ N_LLFFM,
                           outcome == "RAFFM" ~ N_RAFFM,
                           outcome == "LAFFM" ~ N_LAFFM,
                           outcome == "TFFM" ~ N_TFFM, 
                           outcome == "RHGS" ~ N_RHGS,
                           outcome == "LHGS" ~ N_LHGS,
                           outcome == "LS" ~ N_LS))
dat_bytest = sapply(
  unique(dat_steiger$test),
  FUN = function(x) {
    dat_steiger %>% filter(test == x)
  },
  simplify = FALSE,
  USE.NAMES = TRUE
)

r_exp <- list(NA)
r_out <- list(NA)
for(i in 1:length(dat_bytest)){
  r_exp[[i]] <- get_r_from_bsen(dat_bytest[[i]]$beta.exposure,
                                    dat_bytest[[i]]$se.exposure, 
                                    dat_bytest[[i]]$n_exp)
  r_out[[i]] <- get_r_from_bsen(dat_bytest[[i]]$beta.outcome,
                                    dat_bytest[[i]]$se.outcome, 
                                    dat_bytest[[i]]$n_out)
}

steiger_res <- list(NA)
for(i in 1:length(dat_bytest)){
  steiger_res[[i]] <- mr_steiger(p_exp = dat_bytest[[i]]$pval.exposure, 
                                 p_out = dat_bytest[[i]]$pval.outcome, 
                                 n_exp = dat_bytest[[i]]$n_exp, 
                                 n_out = dat_bytest[[i]]$n_out, 
                                 r_exp = r_exp[[i]], 
                                 r_out = r_out[[i]], 
                                 r_xxo = 1, r_yyo = 1)
}
names(steiger_res) <- names(dat_bytest)
steiger_res2 <- as.data.frame(array(NA,c(length(steiger_res),3)))
for(i in 1:length(steiger_res)){
  steiger_res2[,1] <- names(steiger_res)
  steiger_res2[i,2] <- steiger_res[[i]]$correct_causal_direction_adj
  steiger_res2[i,3] <- steiger_res[[i]]$steiger_test_adj
}
colnames(steiger_res2) <- c("test","direction","steiger_p")
steiger_res2 = steiger_res2 %>% separate(test, into = c('outcome', 'exposure'))







