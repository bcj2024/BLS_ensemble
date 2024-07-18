source("specialist_source_pcaonly.R")
source("specialist_source_featureEnsemble.R")

library(mice)
library(dplyr)
library(caret)
library(lme4)

pre_proc = function(outcome){
  browser()
  pre_lst = lapply(emo_ema_lst, function(x){
    mean_score = mean(x[[outcome]], na.rm = T)
    event = (x[[outcome]] - mean_score)>0.5
    x$event = event
    x[!is.na(event),]
  })
  
  final_lst = pre_lst[sapply(pre_lst, function(x) mean(x$event))>0.1 & sapply(pre_lst, nrow)>10 & sapply(pre_lst, function(x) sum(x$event)) > 4]
  final_lst_merge = do.call(rbind, final_lst)
  
  final_pred_cov = final_lst_merge %>% select(starts_with("GPS"), starts_with("Phone"), pow_daily, acc_daily,
                                              weekday, starts_with("Watch")) %>% select(-Watch_acceptedDays)
  formula.pca = as.formula(sprintf("~%s", paste(names(final_pred_cov), collapse = "+")))
  
  final_pred_imp = mice(final_pred_cov)
  final_pred = complete(final_pred_imp)
  final_pred$ID = final_lst_merge$Subject
  final_pred$event = (final_lst_merge$event)
  final_pred = final_pred[final_pred$Watch_ActivityDifferenceWakeSleep != Inf,]
  final_pred_lst = split(final_pred, final_pred$ID)
  final_pred_lst = final_pred_lst[sapply(final_pred_lst, function(x) sum(x$event) > 4)]
  
  final_pred_lst
}

get_comparison = function(res_w, file_raw){
  anx.raw.names = colnames(readRDS(file_raw)[[1]]$w)
  anx.raw.auc = summarize.analysis(readRDS(file_raw))
  anx.raw = anx.raw.auc$ind[[4]]
  names(anx.raw) = anx.raw.names
  
  cbind(res_w, anx.raw[names(res_w)])
}

emo_ema_lst = readRDS("~/../Dropbox/McLean/Christian_K23/Justin_data_res/Robj/all_data.rds")

use_idx = sapply(emo_ema_lst, function(x){
  mean(is.na(x %>% select(starts_with("Watch"))))
}) < 0.8

emo_ema_lst = emo_ema_lst[use_idx]

# anxious
anx_data = pre_proc("Anxious")
anx.w.res = analysis.wrapper.both("event", formula.pca, anx_data, cv = "regular", n.cpus = 8)
anx.w.auc = summarize.analysis(anx.w.res)$ind[[4]]
names(anx.w.auc) = names(anx_data)

plot( get_comparison(anx.w.auc, "Robj/anx.rds"), xlab = "PDEM with watch data",
      ylab = "PDEM without watch data", main = "Anxious")
abline(a = 0, b = 1, lty = 2)


# irritable
ir_data = pre_proc("Irritable")
ir.w.res = analysis.wrapper.both("event", formula.pca, ir_data, cv = "regular", n.cpus = 8)
ir.w.auc = summarize.analysis(ir.all.res)$ind[[4]]
names(ir.w.auc) = names(ir_data)

plot( get_comparison(ir.w.auc, "Robj/ir.rds"))
abline(a = 0, b = 1, lty = 2)
 
# upset
upset_data = pre_proc("Upset")
up.w.res = analysis.wrapper.both("event", formula.pca, upset_data, cv = "regular", n.cpus = 8)
up.w.auc = summarize.analysis(up.w.res)$ind[[4]]
names(up.w.auc) = names(upset_data)

plot( get_comparison(up.w.auc, "Robj/up.rds"))
abline(a = 0, b = 1, lty = 2)

# loneliness
lonely_data = pre_proc("Lonely")
lo.w.res = analysis.wrapper.both("event", formula.pca, lonely_data, cv = "regular", n.cpus = 8)
lo.w.auc = summarize.analysis(lo.w.res)$ind[[4]]
names(lo.w.auc) = names(lonely_data)

pdf("FeatureEnsemble/watch_data_res.pdf", width = 10, height = 8)
par(mfrow = c(2,2))
plot( get_comparison(anx.w.auc, "Robj/anx.rds"), xlab = "PDEM with watch data",
      ylab = "PDEM without watch data", main = "Subject-specific AUC: Anxious")
abline(a = 0, b = 1, lty = 2)
plot( get_comparison(ir.w.auc, "Robj/ir.rds"), xlab = "PDEM with watch data",
      ylab = "PDEM without watch data", main = "Subject-specific AUC: Irritable")
abline(a = 0, b = 1, lty = 2)
plot( get_comparison(up.w.auc, "Robj/up.rds"), xlab = "PDEM with watch data",
      ylab = "PDEM without watch data", main = "Subject-specific AUC: Upset")
abline(a = 0, b = 1, lty = 2)
plot( get_comparison(lo.w.auc, "Robj/lo.rds"), xlab = "PDEM with watch data",
      ylab = "PDEM without watch data", main = "Subject-specific AUC: Lonely")
abline(a = 0, b = 1, lty = 2)
dev.off()
