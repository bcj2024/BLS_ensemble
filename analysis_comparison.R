source("specialist_source_pcaonly.R")

# readin both datasets
anx.pca = readRDS("./Robj/anx_pcaonly.rds")
ir.pca = readRDS("./Robj/ir_pcaonly.rds")
up.pca = readRDS("./Robj/up_pcaonly.rds")
lo.pca = readRDS("./Robj/lo_pcaonly.rds")

anx.fe = readRDS("./Robj/anx.rds")
ir.fe = readRDS("./Robj/ir.rds")
up.fe = readRDS("./Robj/up.rds")
lo.fe = readRDS("./Robj/lo.rds")

table.auc(summarize.analysis(anx.pca))
table.auc(summarize.analysis(anx.fe))

# get individual AUCs, only focus on RF
two_methods_compare_wrapper = function(pca_res, fe_res, emotion){
  auc_anx = get.ind(pca_res[1:3], "Anxious", metric = "AUC")
  auc_fe_anx = get.ind(fe_res[1:3], "Anxious", metric = "AUC")
  
  auc_both = left_join( auc_anx %>% filter(Method == "RF") %>% select(ID, ens), 
                        auc_fe_anx %>% filter(Method == "RF") %>% select(ID, ens), by = "ID" )
  auc_both$diff = auc_both$ens.y - auc_both$ens.x
  auc_both$ID = sapply(auc_both$ID, function(x) names_map[[x]])
  
  p1 = ggplot(auc_both, aes(x = ens.x, y = ens.y)) + geom_point() +
    geom_abline(slope = 1, intercept = 0, linetype = 2) + theme_bw() +
    xlab("PCA only") + ylab("Feature Ensemble") + geom_text(aes(x = ens.x, y = ens.y*1.02, label = ID)) +
    ggtitle(sprintf("%s", emotion))
  
  # pca only weights
  pca_only_w = do.call(rbind, lapply( seq_along(pca_res[[3]]$ind), function(x){
    w_raw = pca_res[[3]]$ind[[x]]$w
    w = rep(0, length(w_raw))
    w[x] = w_raw[1]
    w[-x] = w_raw[-1]
    w
  } ))
  id = sapply(pca_res[[3]]$ind, function(x) x$data$ID[1])
  row.names(pca_only_w) = id
  colnames(pca_only_w) = id
  plt_pcaonly = reshape2::melt(pca_only_w)
  plt_pcaonly$type = "PCA only"
  colnames(plt_pcaonly)[1:3] = c("Ensembled", "Targeted", "Weight")
  
  # feature ensemble
  w_converted = lapply(seq_along(fe_res[[3]]$ind), function(idx){
    w_all = fe_res[[3]]$ind[[idx]]$w
    tmp1 = w_all[seq(1,length(w_all),2)]
    tmp2 = w_all[-seq(1,length(w_all),2)]
    
    pca_w = tmp1
    pca_w[-idx] = tmp1[-1]
    pca_w[idx] = tmp1[1]
    
    ind_w = tmp2
    ind_w[-idx] = tmp2[-1]
    ind_w[idx] = tmp2[1]
    
    list(pca = pca_w, ind = ind_w, ID = fe_res[[3]]$ind[[idx]]$data$ID[1])
  })
  pca_w_all = sapply(w_converted, '[[', "pca")
  ind_w_all = sapply(w_converted, '[[', "ind")
  names_all = sapply(w_converted, '[[', 'ID')
  
  colnames(pca_w_all) = names_all
  row.names(pca_w_all) = names_all
  colnames(ind_w_all) = names_all
  row.names(ind_w_all) = names_all
  plt_pca = melt(pca_w_all)
  names(plt_pca) = c("Ensembled", "Targeted", "Weight")
  plt_pca$type = "PCA"
  plt_ind = melt(ind_w_all)
  names(plt_ind) = c("Ensembled", "Targeted", "Weight")
  plt_ind$type = "All Features"
  plt_all = rbind(plt_pcaonly, plt_pca, plt_ind)
  plt_all$type = factor(plt_all$type, levels = c("PCA only", "PCA", "All Features"))
  
  p2 = ggplot(plt_all, aes(x = Ensembled, y = Targeted, fill = Weight)) + geom_tile(color = "black") + facet_grid(.~type) +
    coord_fixed() + theme_bw() + scale_fill_gradient2(low = "blue", high = "red") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlab("Ensembled Subject") + ylab("Targeted Subject") +
    ggtitle(sprintf("Ensemble Weights: %s", emotion)) + geom_abline(slope = 1, intercept = 0, linetype = 2) +
    scale_x_discrete(expand = c(0,0)) + scale_y_discrete(expand = c(0,0))
  list( plot = list(p1, p2),
        weights = plt_all )
}

res_anx = two_methods_compare_wrapper(anx.pca, anx.fe, "Anxious")
res_ir = two_methods_compare_wrapper(ir.pca, ir.fe, "Irritable")
res_up = two_methods_compare_wrapper(up.pca, up.fe, "Upset")
res_lo = two_methods_compare_wrapper(lo.pca, lo.fe, "Lonely")

pdf("FeatureEnsemble/all_twomethods_compare.pdf", width = 12, height = 10)
gridExtra::grid.arrange(res_anx$plot[[1]], res_ir$plot[[1]], res_up$plot[[1]],
                        res_lo$plot[[1]], nrow = 2)
dev.off()

ggsave("FeatureEnsemble/anx_twomethods_compare.pdf", res_anx$plot[[1]], width = 6, height = 5)
ggsave("FeatureEnsemble/anx_ens_weights.pdf", res_anx$plot[[2]], width = 16, height = 6.5)
ggsave("FeatureEnsemble/ir_twomethods_compare.pdf", res_ir$plot[[1]], width = 6, height = 5)
ggsave("FeatureEnsemble/ir_ens_weights.pdf", res_ir$plot[[2]], width = 16, height = 6.5)
ggsave("FeatureEnsemble/up_twomethods_compare.pdf", res_up$plot[[1]], width = 6, height = 5)
ggsave("FeatureEnsemble/up_ens_weights.pdf", res_up$plot[[2]], width = 16, height = 6.5)
ggsave("FeatureEnsemble/lo_twomethods_compare.pdf", res_lo$plot[[1]], width = 6, height = 5)
ggsave("FeatureEnsemble/lo_ens_weights.pdf", res_lo$plot[[2]], width = 16, height = 6.5)
