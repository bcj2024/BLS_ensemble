source("specialist_source_pcaonly.R")
library(mice)
library(dplyr)
library(caret)
library(lme4)
library(pROC)
library(e1071)

emo_ema_lst = readRDS("./Robj/all_data.rds")

# define events anxious
pre_lst = lapply(emo_ema_lst, function(x){
  mean_score = mean(x$Anxious, na.rm = T)
  event = (x$Anxious - mean_score)>0.5
  x$event = event
  x[!is.na(event),]
})

final_lst = pre_lst[sapply(pre_lst, function(x) mean(x$event))>0.1 & sapply(pre_lst, nrow)>10]
final_lst_merge = do.call(rbind, final_lst)

final_pred_cov = final_lst_merge %>% select(starts_with("GPS"), starts_with("Phone"), pow_daily, acc_daily,
                                            weekday)
final_pred_imp = mice(final_pred_cov)
final_pred = complete(final_pred_imp)
final_pred$ID = final_lst_merge$Subject
final_pred$event = (final_lst_merge$event)
final_pred_lst = split(final_pred, final_pred$ID)

formula.pca = as.formula(sprintf("~%s", paste(names(final_pred_cov), collapse = "+")))

set.seed(1)
anx.all.res = analysis.wrapper("event", formula.pca, final_pred_lst, cv = "regular", n.cpus = 8)

# define events irritable
pre_lst_ir = lapply(emo_ema_lst, function(x){
  mean_score = mean(x$Irritable, na.rm = T)
  event = (x$Irritable - mean_score)>0.5
  x$event = event
  x[!is.na(event),]
})
final_lst_ir = pre_lst_ir[sapply(pre_lst_ir, function(x) mean(x$event))>0.1 & sapply(pre_lst_ir, nrow)>10]
final_lst_merge_ir = do.call(rbind, final_lst_ir)

final_pred_cov_ir = final_lst_merge_ir %>% select(starts_with("GPS"), starts_with("Phone"),
                                                  pow_daily, acc_daily, weekday)
final_pred_imp_ir = mice(final_pred_cov_ir)
final_pred_ir = complete(final_pred_imp_ir)
final_pred_ir$ID = final_lst_merge_ir$Subject
final_pred_ir$event = final_lst_merge_ir$event
final_pred_ir_lst = split(final_pred_ir, final_pred_ir$ID)

ir.all.res = analysis.wrapper("event", formula.pca, final_pred_ir_lst, cv = "regular", n.cpus = 8)

#### define upset ####
pre_lst_up = lapply(emo_ema_lst, function(x){
  mean_score = mean(x$Upset, na.rm = T)
  event = (x$Upset - mean_score)>0.5
  x$event = event
  x[!is.na(event),]
})
final_lst_up = pre_lst_up[sapply(pre_lst_up, function(x) mean(x$event))>0.1 & sapply(pre_lst_up, nrow)>10]
final_lst_merge_up = do.call(rbind, final_lst_up)
final_pred_cov_up = final_lst_merge_up %>% select(starts_with("GPS"), starts_with("Phone"),
                                                  pow_daily, acc_daily, weekday)

final_pred_imp_up = mice(final_pred_cov_up)
final_pred_up = complete(final_pred_imp_up)
final_pred_up$ID = final_lst_merge_up$Subject
final_pred_up$event = (final_lst_merge_up$event)
final_pred_up_lst = split(final_pred_up, final_pred_up$ID)

up.all.res = analysis.wrapper("event", formula.pca, final_pred_up_lst, cv = "regular", n.cpus = 8)

#### define Lonely ####
pre_lst_lo = lapply(emo_ema_lst, function(x){
  mean_score = mean(x$Lonely, na.rm = T)
  event = (x$Lonely - mean_score)>0.5
  x$event = event
  x[!is.na(event),]
})
final_lst_lo = pre_lst_lo[sapply(pre_lst_lo, function(x) mean(x$event))>0.1 & sapply(pre_lst_lo, nrow)>10]
final_lst_merge_lo = do.call(rbind, final_lst_lo)
final_pred_cov_lo = final_lst_merge_lo %>% select(starts_with("GPS"), starts_with("Phone"),
                                                  pow_daily, acc_daily, weekday)

library(mice)
final_pred_imp_lo = mice(final_pred_cov_lo)
final_pred_lo = complete(final_pred_imp_lo)
final_pred_lo$ID = final_lst_merge_lo$Subject
final_pred_lo$event = (final_lst_merge_lo$event)
final_pred_lo_lst = split(final_pred_lo, final_pred_lo$ID)

lo.all.res = analysis.wrapper("event", formula.pca, final_pred_lo_lst, cv = "regular", n.cpus = 8)

#save results
saveRDS(anx.all.res, "./Robj/anx_pcaonly.rds")
saveRDS(ir.all.res, "./Robj/ir_pcaonly.rds")
saveRDS(up.all.res, "./Robj/up_pcaonly.rds")
saveRDS(lo.all.res, "./Robj/lo_pcaonly.rds")