source("specialist_source_pcaonly.R")
source("specialist_source_featureEnsemble.R")
library(dplyr)
library(mice)
library(dplyr)
library(caret)
library(lme4)
library(pROC)
library(e1071)

emo_ema_lst = readRDS("./Robj/all_data.rds")

# deidentify
names_map = as.list(paste("SUB", seq_along(emo_ema_lst), sep = ""))
names(names_map) = names(emo_ema_lst)

# define events anxious
pre_lst = lapply(emo_ema_lst, function(x){
  mean_score = mean(x$Anxious, na.rm = T)
  event = (x$Anxious - mean_score)>0.5
  x$event = event
  x[!is.na(event),]
})

# filtering based on HNA rate
final_lst = pre_lst[sapply(pre_lst, function(x) mean(x$event))>0.1 & sapply(pre_lst, nrow)>10 & sapply(pre_lst, function(x) sum(x$event)) > 4]
# summary
summary_raw(pre_lst, "Anxious")

# imputation
final_lst_merge = do.call(rbind, final_lst)
final_pred_cov = final_lst_merge %>% select(starts_with("GPS"), starts_with("Phone"), pow_daily, acc_daily,
                                            weekday)
formula.pca = as.formula(sprintf("~%s", paste(names(final_pred_cov), collapse = "+")))
final_pred_imp = mice(final_pred_cov)
final_pred = complete(final_pred_imp)
final_pred$ID = sapply(final_lst_merge$Subject, function(x) names_map[[x]])
final_pred$event = (final_lst_merge$event)
final_pred$Anxious = final_lst_merge$Anxious
final_pred_lst = split(final_pred, final_pred$ID)

# raw cor plot
p1 = raw_cor_visualization(final_pred_lst, "Anxious")

# speciliast prediction
set.seed(1)
anx.all.res = analysis.wrapper.both("event", formula.pca, final_pred_lst, cv = "regular", n.cpus = 8)
anx.feature = get_pca_res(final_pred_lst, anx.all.res, anx.only.res, "Anxious", nclust = 2:5)

pdf("./FeatureEnsemble/anx_feature_importance.pdf", width = 9, height = 6)
anx.feature$plot + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

pdf("./FeatureEnsemble/anx_ens_weights_1.pdf", width = 14, height = 7)
anx.feature$plot_weight
dev.off()

# define events irritable
pre_lst_ir = lapply(emo_ema_lst, function(x){
  mean_score = mean(x$Irritable, na.rm = T)
  event = (x$Irritable - mean_score)>0.5
  x$event = event
  x[!is.na(event),]
})

# summary
summary_raw(pre_lst_ir, "Irritable")

final_lst_ir = pre_lst_ir[sapply(pre_lst_ir, function(x) mean(x$event))>0.1 & sapply(pre_lst_ir, nrow)>10 & sapply(pre_lst_ir, function(x) sum(x$event)) > 4]
final_lst_merge_ir = do.call(rbind, final_lst_ir)

final_pred_cov_ir = final_lst_merge_ir %>% select(starts_with("GPS"), starts_with("Phone"),
                                                  pow_daily, acc_daily, weekday)
final_pred_imp_ir = mice(final_pred_cov_ir)
final_pred_ir = complete(final_pred_imp_ir)
final_pred_ir$ID = sapply(final_lst_merge_ir$Subject, function(x) names_map[[x]])
final_pred_ir$event = final_lst_merge_ir$event
final_pred_ir$Irritable = final_lst_merge_ir$Irritable
final_pred_ir_lst = split(final_pred_ir, final_pred_ir$ID)

# raw cor plot
p2 = raw_cor_visualization(final_pred_ir_lst, "Irritable")

# specialist prediction
ir.all.res = analysis.wrapper.both("event", formula.pca, final_pred_ir_lst, cv = "regular", n.cpus = 8)
ir.feature = get_pca_res(final_pred_ir_lst, ir.all.res, ir.only.res, "Irritable", 3)
pdf("./FeatureEnsemble/ir_feature_importance.pdf", width = 9, height = 6)
ir.feature$plot + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

pdf("./FeatureEnsemble/ir_ens_weights.pdf", width = 14, height = 7)
ir.feature$plot_weight
dev.off()

#### define upset ####
pre_lst_up = lapply(emo_ema_lst, function(x){
  mean_score = mean(x$Upset, na.rm = T)
  event = (x$Upset - mean_score)>0.5
  x$event = event
  x[!is.na(event),]
})

# summary
summary_raw(pre_lst_up, "Upset")

final_lst_up = pre_lst_up[sapply(pre_lst_up, function(x) mean(x$event))>0.1 & sapply(pre_lst_up, nrow)>10 & sapply(pre_lst_up, function(x) sum(x$event)) > 4]
final_lst_merge_up = do.call(rbind, final_lst_up)
final_pred_cov_up = final_lst_merge_up %>% select(starts_with("GPS"), starts_with("Phone"),
                                                  pow_daily, acc_daily, weekday)

final_pred_imp_up = mice(final_pred_cov_up)
final_pred_up = complete(final_pred_imp_up)
final_pred_up$ID = sapply(final_lst_merge_up$Subject, function(x) names_map[[x]])
final_pred_up$event = (final_lst_merge_up$event)
final_pred_up$Upset = final_lst_merge_up$Upset
final_pred_up_lst = split(final_pred_up, final_pred_up$ID)

# raw cor plot
p3 = raw_cor_visualization(final_pred_up_lst, "Upset")

# specialist prediction
up.all.res = analysis.wrapper.both("event", formula.pca, final_pred_up_lst, cv = "regular", n.cpus = 8)
up.feature = get_pca_res(final_pred_up_lst, up.all.res, up.only.res, "Upset", 2)
pdf("./FeatureEnsemble/up_feature_importance.pdf", width = 9, height = 6)
up.feature$plot + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

pdf("./FeatureEnsemble/up_ens_weights.pdf", width = 14, height = 7)
up.feature$plot_weight
dev.off()

#### define Lonely ####
pre_lst_lo = lapply(emo_ema_lst, function(x){
  mean_score = mean(x$Lonely, na.rm = T)
  event = (x$Lonely - mean_score)>0.5
  x$event = event
  x[!is.na(event),]
})

# summary
summary_raw(pre_lst_lo, "Lonely")

final_lst_lo = pre_lst_lo[sapply(pre_lst_lo, function(x) mean(x$event))>0.1 & sapply(pre_lst_lo, nrow)>10 & sapply(pre_lst_lo, function(x) sum(x$event)) > 4]
final_lst_merge_lo = do.call(rbind, final_lst_lo)
final_pred_cov_lo = final_lst_merge_lo %>% select(starts_with("GPS"), starts_with("Phone"),
                                                  pow_daily, acc_daily, weekday)

final_pred_imp_lo = mice(final_pred_cov_lo)
final_pred_lo = complete(final_pred_imp_lo)
final_pred_lo$ID = sapply(final_lst_merge_lo$Subject, function(x) names_map[[x]])
final_pred_lo$event = (final_lst_merge_lo$event)
final_pred_lo$Lonely = final_lst_merge_lo$Lonely
final_pred_lo_lst = split(final_pred_lo, final_pred_lo$ID)

# raw cor plot
p4 = raw_cor_visualization(final_pred_lo_lst, "Lonely")

# Specialist predictions
lo.all.res = analysis.wrapper.both("event", formula.pca, final_pred_lo_lst, cv = "regular", n.cpus = 8)
lo.feature = get_pca_res(final_pred_lo_lst, lo.all.res, lo.only.res, "Lonely")
pdf("./FeatureEnsemble/lo_feature_importance.pdf", width = 9, height = 6)
lo.feature$plot + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

pdf("./FeatureEnsemble/lo_ens_weights.pdf", width = 14, height = 7)
lo.feature$plot_weight
dev.off()

#save results
saveRDS(anx.all.res, "./Robj/anx.rds")
saveRDS(ir.all.res, "./Robj/ir.rds")
saveRDS(up.all.res, "./Robj/up.rds")
saveRDS(lo.all.res, "./Robj/lo.rds")

#### Additional analysis ####
anx.all.auc = summarize.analysis(anx.all.res)
ir.all.auc = summarize.analysis(ir.all.res)
up.all.auc = summarize.analysis(up.all.res)
lo.all.auc = summarize.analysis(lo.all.res)

# optimal cut-off sen and spe
pdf("./FeatureEnsemble/final_anx_roc_featureEns.pdf", width = 6, height = 6)
plot.auc(summarize.analysis(anx.all.res)) + ggtitle("Anxious")
dev.off()

pdf("./FeatureEnsemble/final_ir_roc_featureEns.pdf", width = 6, height = 6)
plot.auc(summarize.analysis(ir.all.res)) + ggtitle("Irritable")
dev.off()

pdf("./FeatureEnsemble/final_up_roc_featureEns.pdf", width = 6, height = 6)
plot.auc(summarize.analysis(up.all.res)) + ggtitle("Upset")
dev.off()

pdf("./FeatureEnsemble/final_lo_roc_featureEns.pdf", width = 6, height = 6)
plot.auc(summarize.analysis(lo.all.res)) + ggtitle("Lonely")
dev.off()

# AUC/Acc table
rbind( table.auc(summarize.analysis(anx.all.res)),
       table.auc(summarize.analysis(ir.all.res)),
       table.auc(summarize.analysis(up.all.res)),
       table.auc(summarize.analysis(lo.all.res)) )

rbind( table.acc(anx.all.res),
       table.acc(ir.all.res),
       table.acc(up.all.res),
       table.acc(lo.all.res))

# subject specific auc
AUC_compare = rbind( get.ind.metric(anx.all.res[1:3], "Anxious", metric = "AUC"),
                     get.ind.metric(ir.all.res[1:3], "Irritable", metric = "AUC"),
                     get.ind.metric(up.all.res[1:3], "Upset", metric = "AUC"),
                     get.ind.metric(lo.all.res[1:3], "Lonely", metric = "AUC"))
AUC_compare$Emotion = factor(AUC_compare$Emotion, levels = c("Anxious", "Irritable", "Upset", "Lonely"))
AUC_compare$Method = factor(AUC_compare$Method, levels = c("ENet", "SVM", "RF"))

pdf("./FeatureEnsemble/PEM_ID_compare_AUC_featureEns.pdf", width = 6, height = 3)
ggplot(AUC_compare, aes(x = Emotion, y = Diff, fill = Method)) + geom_boxplot(outlier.shape = NA) +
  geom_hline(yintercept = 0, linetype = 2) + theme_bw() + scale_fill_discrete() +
  ylab("AUC(PEM) - AUC(ID)") + ggtitle("Difference in AUC for PEM vs. Idiosyncratic Models")
dev.off()

MSE_compare = rbind( get.ind.metric(anx.all.res[1:3], "Anxious", metric = "MSE"),
                     get.ind.metric(ir.all.res[1:3], "Irritable", metric = "MSE"),
                     get.ind.metric(up.all.res[1:3], "Upset", metric = "MSE"),
                     get.ind.metric(lo.all.res[1:3], "Lonely", metric = "MSE"))
MSE_compare$Emotion = factor(MSE_compare$Emotion, levels = c("Anxious", "Irritable", "Upset", "Lonely"), ordered = T)
MSE_compare$Method = factor(MSE_compare$Method, levels = c("ENet", "SVM", "RF"), ordered = T)

pdf("./FeatureEnsemble/PEM_ID_compare_Brier.pdf", width = 6, height = 3)
ggplot(MSE_compare, aes(x = Emotion, y = -Diff, fill = Method)) + geom_boxplot(outlier.shape = NA) +
  geom_hline(yintercept = 0, linetype = 2) + theme_bw() + scale_fill_discrete() + ylim(c(-0.05,0.1)) +
  ylab("Brier(ID) - Brier(PEM)") + ggtitle("Difference in Brier Score for PEM vs. Idiosyncratic Models")
dev.off()

pdf("./FeatureEnsemble/raw_correlation.pdf", width = 18, height = 18)
gridExtra::grid.arrange(p1, p2, p3, p4, ncol = 1)
dev.off()
