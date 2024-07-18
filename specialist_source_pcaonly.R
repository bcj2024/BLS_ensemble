require(glmnet)
library(e1071)
library(glmnetUtils)
library(alabama)
library(reshape2)
library(lme4)
library(tidyverse)
library(dplyr)
library(mice)
library(caret)
library(PRROC)
library(pROC)

summary_raw = function(data, outcome){
  # event proportion
  sapply(data, function(x) mean(x$event))
  
  # deviation of raw outcome between HNA and normal
  rowMeans(sapply(data, function(x){
    mean_score = mean(x[[outcome]], na.rm = T)
    c(sum(!is.na(x[[outcome]])),
      mean(x$event, na.rm = T),
      mean((x[[outcome]] - mean_score)[x$event==1], na.rm = T),
      mean((x[[outcome]] - mean_score)[x$event==1], na.rm = T)/sd(x[[outcome]], na.rm = T),
      sd((x[[outcome]] - mean_score)[x$event==1], na.rm = T))
  }), na.rm = T)
}

raw_cor_visualization = function(data, outcome){
  cor_all = sapply(data, function(x){
    cor(x[[outcome]], x[,1:15])
  })
  row.names(cor_all) = colnames(final_pred_lst[[1]])[1:15]
  import_cor = t(t(cor_all)/colSums(abs(cor_all), na.rm = T))
  
  dist_cor = dist(t(import_cor))
  cluster_res = sapply(2:5, function(k) cluster::pam(dist_cor, k = k)$silinfo$avg.width)
  cluster_id = cluster::pam(dist_cor, k = (2:5)[order(cluster_res, decreasing = T)[1]])$cluster
  
  import_cor = import_cor[order(rowMeans(abs(import_cor), na.rm = T)), order(cluster_id)]
  
  ggplot(reshape2::melt(import_cor), aes(x = Var2, y = Var1, fill = value) ) + geom_tile(color = "black") +
    scale_fill_gradient2(low = "darkblue", high = "darkred") + coord_fixed() + 
    scale_x_discrete(expand = c(0,0)) + scale_y_discrete(expand = c(0,0)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlab("") + ylab("") + 
    ggtitle(sprintf("%s: original correlation"), outcome) + 
    geom_vline(xintercept = cumsum(table(cluster_id) + 0.5), linewidth = 2)
}

#pca res
get_pca_res = function(data, res, res_only, title, nclust = 2:5){
  # browser()
  library(snowfall)
  sfInit(parallel = T, cpus = 8)
  pca_anx = pca.preproc.import(data, "event", formula.pca)
  
  # In model feature importance
  pca_import_in = sfLapply(pca_anx$data, function(ind_x){
    rf_fit = randomForest::randomForest(as.factor(event) ~ .-ID, data = ind_x)
    import_sign = sign(sapply(1:5, function(idx) cor(pdp::partial(rf_fit, idx, which.class = "TRUE"))[1,2]))
    import_sign[is.na(import_sign)] = 0
    import_norm = c(rf_fit$importance/sum(rf_fit$importance))*import_sign
    final_import = rowSums(t(t(pca_anx$loadings)/colSums(abs(pca_anx$loadings))*import_norm))
    final_import/sum(abs(final_import))
  })
  
  ind_import_in = sfLapply(data, function(x){
    rf_fit = randomForest::randomForest(as.factor(event) ~ .-ID, data = x)
    import_sign = sign(sapply(1:15, function(idx) cor(pdp::partial(rf_fit, idx, which.class = "TRUE"))[1,2]))
    import_sign[is.na(import_sign)] = 0
    import_norm = c(rf_fit$importance/sum(rf_fit$importance))*import_sign
    import_norm
  })
  
  subject_lst = sapply(res[[1]]$ind, function(x){
    names_map[[x$data$ID[1]]]
  })
  
  pca_im = do.call(cbind, pca_import_in[subject_lst])
  ind_im = do.call(cbind, ind_import_in[subject_lst])
  
  in_im_combined = matrix(nrow = nrow(pca_im), ncol = 2*ncol(pca_im))
  in_im_combined[,seq(1,ncol(in_im_combined),2)] = pca_im
  in_im_combined[,-seq(1,ncol(in_im_combined),2)] = ind_im
  
  final_imp = sapply(seq_along(res[[3]]$ind), function(idx){
    in_idx = (idx-1)*2 + 1:2
    tmp_im = cbind(in_im_combined[,in_idx], in_im_combined[,-in_idx])
    tmp_im %*% res[[3]]$ind[[idx]]$w
  })
  colnames(final_imp) = subject_lst
  row.names(final_imp) = colnames(data[[1]])[1:15]
  
  import_dist = dist(t(final_imp))
  cluster_select = sapply( nclust, function(k) cluster::pam(import_dist, k = k)$silinfo$avg.width )
  cluster_res = cluster::pam(import_dist, k = nclust[order(cluster_select, decreasing = T)][1])
  
  df_cluster = data.frame(study_id = names(cluster_res$clustering), cluster = cluster_res$clustering)
  
  order_col = order(cluster_res$clustering)
  order_row = order(rowMeans(abs(final_imp)))
  
  # also collect weights
  w_converted = lapply(seq_along(res[[3]]$ind), function(idx){
    w_all = res[[3]]$ind[[idx]]$w
    tmp1 = w_all[seq(1,length(w_all),2)]
    tmp2 = w_all[-seq(1,length(w_all),2)]
    
    pca_w = tmp1
    pca_w[-idx] = tmp1[-1]
    pca_w[idx] = tmp1[1]
    
    ind_w = tmp2
    ind_w[-idx] = tmp2[-1]
    ind_w[idx] = tmp2[1]
    
    list(pca = pca_w, ind = ind_w, ID = names_map[[res[[3]]$ind[[idx]]$data$ID[1]]])
  })
  
  w_only = lapply(seq_along(res_only[[3]]$ind), function(idx){
    tmp1 = res_only[[3]]$ind[[idx]]$w
    pca_w = tmp1
    pca_w[-idx] = tmp1[-1]
    pca_w[idx] = tmp1[1]
    
    list(pca = pca_w, ID = names_map[[res_only[[3]]$ind[[idx]]$data$ID[1]]])
  })
  
  # cols are targeted ind, rows are ensembled ind
  pca_w_all = sapply(w_converted, '[[', "pca")
  ind_w_all = sapply(w_converted, '[[', "ind")
  names_all = sapply(w_converted, '[[', 'ID')
  
  pca_w_only = sapply(w_only, '[[', 'pca')
  
  colnames(pca_w_all) = names_all
  row.names(pca_w_all) = names_all
  colnames(ind_w_all) = names_all
  row.names(ind_w_all) = names_all
  colnames(pca_w_only) = names_all
  row.names(pca_w_only) = names_all
  
  plt_pca = melt(pca_w_all)
  names(plt_pca) = c("Ensembled", "Targeted", "Weight")
  plt_pca$type = "PCA-based Models"
  plt_ind = melt(ind_w_all)
  names(plt_ind) = c("Ensembled", "Targeted", "Weight")
  plt_ind$type = "Raw-feature-based Models"
  plt_only = melt(pca_w_only)
  names(plt_only) = c("Ensembled", "Targeted", "Weight")
  plt_only$type = "PCA only"
  plt_all = rbind(plt_pca, plt_ind, plt_only)
  plt_all$type = factor(plt_all$type, levels = c("PCA only", "PCA-based Models", "Raw-feature-based Models"))
  sfStop()
  
  return( list(plot = ggplot(melt(final_imp[order_row, order_col]), aes(y = Var1, x = Var2, fill = value)) + geom_tile(color = "black") +
                 scale_fill_gradient2(low = "darkblue", high = "darkred") + coord_fixed() +
                 theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_x_discrete(expand = c(0,0)) + scale_y_discrete(expand = c(0,0)) +
                 geom_vline(xintercept = cumsum(table(cluster_res$clustering))+0.5, linewidth = 1.5) +
                 xlab("") + ylab("") + ggtitle(sprintf("Feature Importance: %s", title)) + theme_bw(),
               plot_weight = ggplot(plt_all, aes(x = Ensembled, y = Targeted, fill = Weight)) + geom_tile(color = "black") + facet_grid(.~type) +
                 coord_fixed() + theme_bw() + scale_fill_gradient2(low = "blue", high = "darkred") +
                 theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlab("Ensembled Subject") + ylab("Targeted Subject") +
                 ggtitle(sprintf("Ensemble Weights: %s", title)) + geom_abline(slope = 1, intercept = 0, linetype = 2),
               df_cluster = df_cluster) )
}

logistic.loss = function(y, p.pred){
  n.pos = sum(y==1)
  n.neg = sum(y==0)
  p.pred[p.pred<=0] = 1e-3
  p.pred[p.pred>=1] = 1 - 1e-3
  n.pos = sum(y==1)
  -mean(y*log(p.pred) + (1-y)*log(1-p.pred))
}

w.simplex = function(p.stacking, y, L, lambda){
  # browser()
  K = ncol(p.stacking)/L
  target.fn = function(x){
    # try to maximize C-statistics
    p.pred = p.stacking%*%x
    # p.pos = p.pred[y==1]
    # p.neg = p.pred[y==0]
    # n.pos = sum(y==1)
    # n.neg = sum(y==0)
    # diff = outer(p.pos, p.neg, "-")
    # sum(-((diff>gamma)*(diff-gamma))^2)/length(p.pos)/length(p.neg)
    # sum(log(1+exp(-3*diff)))/length(p.pos)/length(p.neg)
    # sum((1-diff)^2)/length(p.pos)/length(p.neg)
    # sum((y - p.pred)^2)
    logistic.loss(y, p.pred) #+ lambda*sum((x - w.gen.new)^2)
  }
  res = auglag(c(rep(0.9/L,L), rep(0.1/(K-1)/L,(K-1)*L)),
               target.fn,
               hin = function(x) x,
               heq = function(x) sum(x)-1,
               control.outer = list(trace = F))
  res$par
}

enet.wrapper = function(X, y){
  if(sum(y==1)<=3){
    y.tmp = sample(c(1,1,rep(0, length(y)-2)))
    glm(y~., family = "binomial", data = data.frame(y = y.tmp, X))
  }
  glm(y~., family = "binomial", data = data.frame(y = y, X))
}

svm.wrapper = function(X, y){
  svm(x = X, y = y, type = "C-classification", probability = T)
}

rf.wrapper = function(X, y){
  randomForest::randomForest(x = X, y = y)
}

specialist.wrapper = function(X, y, learner.lst, other.models, lambda = 0, nfolds = 10){
  # browser()
  folds = createFolds(y, k = nfolds)

  L = length(learner.lst)

  # create
  p.stacking = do.call(rbind, lapply(folds, function(fold){
    in.models = lapply( learner.lst, function(learner){
      learner(X[-fold,,drop = F], y[-fold])
    })
    all.models = c(in.models, other.models)

    sapply(all.models, function(ind.model){
      if(class(ind.model)[1] == "glm"){
        predict(ind.model, data.frame(X[fold,,drop = F]), type = "response")
      }else if(class(ind.model)[1] == "svm"){
        attributes(predict(ind.model, X[fold,,drop = F], probability = T))$probabilities[,2]
      }else if(class(ind.model)[1] == "randomForest"){
        predict(ind.model, X[fold,,drop = F], type = "prob")[,2]
      }
    })
  }))

  y.obs = y[unlist(folds)]
  # cv.res = cv.glmnet(x = p.stacking, y = y.obs, family = "binomial", alpha = 1)
  # glm.res = glmnet(x = p.stacking, y = y.obs, family = "binomial", alpha = 1, lambda = cv.res$lambda.min)
  # c(glm.res$a0, as.numeric(glm.res$beta))
  # nnls(A = p.stacking, b = as.numeric(y.obs)-1)$x
  w.simplex( p.stacking, as.numeric(y.obs) - 1, L, lambda = lambda)
}

penalized.stack.wrapper = function(data.ls, formula, outcome, learner.lst, n.cpus = 4){
  library(snowfall)
  sfInit(parallel = T, cpus = n.cpus)
  # browser()
  # L is the number of algorithms in consideration
  L = length(learner.lst)
  all.models = do.call(c, lapply(data.ls, function(data.ind){
    # get rid of the intercept
    X = model.matrix(formula, data = data.ind)[,-1]
    y = as.factor(data.ind[[outcome]])
    model.res = lapply(learner.lst, function(learner){
      learner(X, y)
    })
    model.res
  }))

  # specialist training
  sfExportAll()
  sfLibrary(caret)
  sfLibrary(alabama)
  sfLibrary(e1071)
  sfLibrary(randomForest)
  specialist.auc.all = sfLapply(seq_along(data.ls), function(idx){
    print(idx)
    X = model.matrix(formula, data = data.ls[[idx]])[,-1]
    y = as.factor(data.ls[[idx]][[outcome]])

    folds = createFolds(y, k = 10)
    remove.idx = (idx-1)*L + (1:L)

    res.stack = lapply(folds, function(fold){
      print(fold)
      in.models = lapply(learner.lst, function(learner){
        learner(X[-fold,,drop=F], y[-fold])
      })
      w.spe = specialist.wrapper(X[-fold,,drop=F], y[-fold], learner.lst, all.models[-remove.idx], nfolds = 10)

      all.pred = sapply(c(in.models, all.models[-remove.idx]), function(ind.model){
        if(class(ind.model)[1] == "glm"){
          predict(ind.model, data.frame(X[fold,,drop = F]), type = "response")
        }else if(class(ind.model)[1] == "svm"){
          attributes(predict(ind.model, X[fold,,drop = F], probability = T))$probabilities[,2]
        }else if(class(ind.model)[1] == "randomForest"){
          predict(ind.model, X[fold,,drop = F], type = "prob")[,2]
        }
      })
      if(is.null(dim(all.pred))){
        p.in = matrix(all.pred[1:L], nrow = 1)
      }else{
        p.in = all.pred[,1:L,drop = F]
      }
      list( p.pred = all.pred%*%w.spe, w = w.spe, p.in = p.in)
    })

    p.stack = unlist(lapply(res.stack, function(x) x$p.pred))
    p.in = do.call(rbind, lapply(res.stack, function(x) x$p.in))
    w.all = sapply(res.stack, function(x) x$w)
    y.obs = y[unlist(folds)]

    list(data = data.frame(ID = names(data.ls)[idx],
                           y = y.obs, p.stack = p.stack, p.in = p.in), w = rowMeans(w.all))
  })
  sfStop()
  all.pred = do.call(rbind, lapply(specialist.auc.all, function(x) x$data))

  w.all = sapply(seq_along(specialist.auc.all), function(i){
    tmp = specialist.auc.all[[i]]$w
    weight = tmp
    weight[i] = tmp[1]
    weight[-i] = tmp[-1]
    weight
  })
  colnames(w.all) = names(data.ls)
  return(list(res = all.pred, w = w.all, ind = specialist.auc.all))
}

penalized.stack.tsCV.wrapper = function(data.ls, formula, outcome, learner.lst, test.prop = 0.75){
  data.ls = lapply(data.ls, function(data.ind){
    data.ind[order(data.ind$day),]
  })

  # L is the number of algorithms in consideration
  L = length(learner.lst)
  all.models = do.call(c, lapply(data.ls, function(data.ind){
    # get rid of the intercept
    X = model.matrix(formula, data = data.ind)[,-1]
    y = as.factor(data.ind[[outcome]])
    model.res = lapply(learner.lst, function(learner){
      learner(X, y)
    })
    model.res
  }))

  # specialist training
  specialist.auc.all = lapply(seq_along(data.ls), function(idx){
    print(idx)
    X = model.matrix(formula, data = data.ls[[idx]])[,-1]
    y = as.factor(data.ls[[idx]][[outcome]])

    n.all = nrow(X)
    test.idx = (round(test.prop*n.all)):n.all
    remove.idx = (idx-1)*L + (1:L)

    res.stack = lapply(test.idx, function(test.ind){
      train.idx = 1:(test.ind-1)
      in.models = lapply(learner.lst, function(learner){
        learner(X[train.idx,,drop=F], y[train.idx])
      })
      w.spe = specialist.wrapper(X[train.idx,,drop=F], y[train.idx], learner.lst, all.models[-remove.idx], nfolds = 10)

      all.pred = sapply(c(in.models, all.models[-remove.idx]), function(ind.model){
        if(class(ind.model)[1] == "glm"){
          predict(ind.model, data.frame(X[test.ind,,drop = F]), type = "response")
        }else if(class(ind.model)[1] == "svm"){
          attributes(predict(ind.model, X[test.ind,,drop = F], probability = T))$probabilities[,2]
        }else if(class(ind.model)[1] == "randomForest"){
          predict(ind.model, X[test.ind,,drop = F], type = "prob")[,2]
        }
      })
      if(is.null(dim(all.pred))){
        p.in = matrix(all.pred[1:L], nrow = 1)
      }else{
        p.in = all.pred[,1:L,drop = F]
      }
      list( p.pred = all.pred%*%w.spe, w = w.spe, p.in = p.in)
    })

    p.stack = unlist(lapply(res.stack, function(x) x$p.pred))
    p.in = do.call(rbind, lapply(res.stack, function(x) x$p.in))
    w.all = sapply(res.stack, function(x) x$w)
    y.obs = y[test.idx]

    list(data = data.frame(ID = names(data.ls)[idx],
                           y = y.obs, p.stack = p.stack, p.in = p.in), w = rowMeans(w.all))
  })
  all.pred = do.call(rbind, lapply(specialist.auc.all, function(x) x$data))
  all.pred
}

glmer.wrapper = function(data.use, outcome){
  all.folds = lapply(data.use, function(x) createFolds(1:nrow(x), 10))
  vars = colnames(data.use[[1]])
  vars = vars[-c(length(vars)-2, length(vars)-1, length(vars))]
  glmer.formula = as.formula( sprintf( "%s~(1|ID) + %s", outcome, paste(vars, collapse = "+") ) )
  glmer.res = do.call(rbind, lapply(1:10, function(fold){
    print(fold)
    data.train = do.call(rbind, lapply(seq_along(data.use), function(study){
      fold.idx = all.folds[[study]][[fold]]
      data.use[[study]][-fold.idx,]
    }))
    data.test = do.call(rbind, lapply(seq_along(data.use), function(study){
      fold.idx = all.folds[[study]][[fold]]
      data.use[[study]][fold.idx,]
    }))

    glmer.res = glmer(glmer.formula, data = data.train, family = "binomial")
    p.glmer = predict(glmer.res, data.test, type = "response")
    y = data.test[[outcome]]
    data.frame(y = y, p.glmer = p.glmer)
  }))
  glmer.res
}

glmer.tscv.wrapper = function(data.use, outcome, test.prop = 0.75){
  formula.pca = as.formula("~ActScore + GPS_freq + ActScore + missing + homDist+
    radiusMobility + percentHome + numPlaces + Power +
    SleepOnset_time + Wake_time + Sleep_Duration +
    mphuse + hphuse + TimeDay")
  df.ls = pca.preproc(data.use, outcome, formula.pca)
  data.ls = df.ls[sapply(df.ls, function(x) sum(x[[outcome]]))>4]

  data.ls = lapply(data.ls, function(data.ind){
    data.ind[order(data.ind$day),]
  })
  vars = colnames(data.ls[[1]])[1:5]
  glmer.formula = as.formula( sprintf( "%s~(1|ID) + %s", outcome, paste(vars, collapse = "+") ) )

  glmer.res = lapply(seq_along(data.ls), function(idx){
    print(idx)
    data.tmp = data.ls
    n.all = nrow(data.ls[[idx]])
    test.idx = (round(test.prop*n.all)):n.all

    res.cv = sapply(test.idx, function(test.ind){
      train.idx = 1:(test.ind-1)
      data.tmp[[idx]] = data.ls[[idx]][train.idx,]
      data.tmp.merge = do.call(rbind, data.tmp)

      glmer.res = glmer(glmer.formula, family = "binomial", data = data.tmp.merge)
      predict(glmer.res, newdata = data.ls[[idx]][test.ind,,drop = F], type = "response")
    })

    data.frame(ID = names(data.ls)[idx], day = data.ls[[idx]]$day[test.idx],
               y = data.ls[[idx]][[outcome]][test.idx], p.pred = res.cv)
  })
  do.call(rbind, glmer.res)
}

pca.preproc = function(data.use, outcome, formula.pca){
  data.all = do.call(rbind, data.use)
  data.all.mt = model.matrix(formula.pca, data.all)[,-1]
  data.pca = princomp(data.all.mt)
  # n.pcs = sum(cumsum(data.pca$sdev^2/sum(data.pca$sdev^2))<0.9) + 1
  n.pcs = 5
  X.data.pca = data.pca$scores[,1:n.pcs]
  df.tmp = data.frame(X.data.pca, ID = data.all$ID, data.all[[outcome]])
  colnames(df.tmp)[ncol(df.tmp)] = outcome
  return( split(df.tmp, df.tmp$ID) )
}

pca.preproc.import = function(data.use, outcome, formula.pca){
  data.all = do.call(rbind, data.use)
  data.all.mt = model.matrix(formula.pca, data.all)[,-1]
  data.pca = princomp(data.all.mt)
  # n.pcs = sum(cumsum(data.pca$sdev^2/sum(data.pca$sdev^2))<0.9) + 1
  n.pcs = 5
  X.data.pca = data.pca$scores[,1:n.pcs]
  df.tmp = data.frame(X.data.pca, ID = data.all$ID, data.all[[outcome]])
  colnames(df.tmp)[ncol(df.tmp)] = outcome
  return( list(data = split(df.tmp, df.tmp$ID),
               loadings = data.pca$loadings[,1:n.pcs]) )
}

analysis.wrapper = function(outcome, formula.pca, data.use, n.cpus = 4, cv = c("regular", "ts")){
  # browser()
  df.ls = pca.preproc(data.use, outcome, formula.pca)
  formula.data = as.formula(paste(outcome,paste(grep("Comp", colnames(df.ls[[1]]), value = T), collapse = "+"),
                                  sep = "~"))
  data.ls = df.ls[sapply(df.ls, function(x) sum(x[[outcome]]))>4]

  if(cv == "regular"){
    print("glmer")
    glmer.res = glmer.wrapper(data.ls, outcome)
    print("enet")
    enet.res = penalized.stack.wrapper(data.ls, formula.data, outcome, c(enet.wrapper), n.cpus = n.cpus)
    print("svm")
    svm.res = penalized.stack.wrapper(data.ls, formula.data, outcome, c(svm.wrapper), n.cpus = n.cpus)
    print("rf")
    rf.res = penalized.stack.wrapper(data.ls, formula.data, outcome, c(rf.wrapper), n.cpus = n.cpus)
    print("DE")
    de.res = penalized.stack.wrapper(data.ls, formula.data, outcome, c(enet.wrapper, svm.wrapper, rf.wrapper), n.cpus = n.cpus)
  }else if(cv == "ts"){
    glmer.res = NA
    enet.res = penalized.stack.tsCV.wrapper(data.ls, formula.data, outcome, c(enet.wrapper))
    svm.res = penalized.stack.tsCV.wrapper(data.ls, formula.data, outcome, c(svm.wrapper))
    rf.res = penalized.stack.tsCV.wrapper(data.ls, formula.data, outcome, c(rf.wrapper))
    de.res = penalized.stack.tsCV.wrapper(data.ls, formula.data, outcome, c(enet.wrapper, svm.wrapper, rf.wrapper))
  }

  list(enet.res, svm.res, rf.res, de.res, glmer.res)
}

summarize.analysis = function( all.res ){
  # get auc
  # browser()
  auc.merge = lapply(seq_along(all.res), function(i){
    if(i==5){
      roc(all.res[[i]]$y, all.res[[i]]$p.glmer, direction = "<")
    }else if(i==4){
      roc(all.res[[i]]$res$y, all.res[[i]]$res$p.stack, direction = "<")
    }else{
      roc.stack = roc(all.res[[i]]$res$y, all.res[[i]]$res$p.stack, direction = "<")
      roc.ind.all = lapply(grep("p.in", names(all.res[[i]]$res)), function(x){
        roc(all.res[[i]]$res$y, all.res[[i]]$res[[x]], direction = "<")
      })
      roc.ind = roc.ind.all[[order(sapply(roc.ind.all, "[[", "auc"), decreasing = T)[1]]]
      list(roc.stack,
           roc.ind)
    }
  })
  auc.ind = lapply(seq_along(all.res), function(i){
    if(i<4){
      sapply(all.res[[i]]$ind, function(x){
        auc.stack = auc(x$data$y, x$data$p.stack, direction = "<")
        auc.in = max(sapply(grep("p.in|V1", colnames(x$data)), function(idx){
          auc(x$data$y, x$data[,idx], direction = "<")
        }))
        c( auc.stack, auc.in )
      })
    }else if(i==4){
      sapply(all.res[[i]]$ind, function(x){
        auc(x$data$y, x$data$p.stack, direction = "<")
      })
    }
  })
  list(merge = auc.merge, ind = auc.ind)
}

plot.auc = function(all.auc){
  all.length = c(length(all.auc$merge[[1]][[1]]$sensitivities), length(all.auc$merge[[2]][[1]]$sensitivities),
                 length(all.auc$merge[[3]][[1]]$sensitivities), length(all.auc$merge[[4]]$sensitivities),
                 length(all.auc$merge[[5]]$sensitivities))
  all.sen = c(all.auc$merge[[1]][[1]]$sensitivities, all.auc$merge[[2]][[1]]$sensitivities, all.auc$merge[[3]][[1]]$sensitivities,
              all.auc$merge[[4]]$sensitivities, all.auc$merge[[5]]$sensitivities)
  all.spe = c(all.auc$merge[[1]][[1]]$specificities, all.auc$merge[[2]][[1]]$specificities, all.auc$merge[[3]][[1]]$specificities,
              all.auc$merge[[4]]$specificities, all.auc$merge[[5]]$specificities)
  plt.df = data.frame(sensitivity = all.sen, specificity = all.spe,
                      Method = factor(rep(c("PEM-ENet", "PEM-SVM", "PEM-RF", "PDEM", "GLMER"), all.length),
                                      levels = c("GLMER", "PEM-ENet", "PEM-SVM", "PEM-RF", "PDEM"),
                                      ordered = T))
  ggplot(plt.df, aes(x = 1-specificity, y = sensitivity, color = Method)) + geom_line() + scale_color_discrete()+
    theme_bw() + geom_abline(slope = 1, intercept = 0, linetype = 2) + coord_fixed() + xlab("1 - Specificity") + ylab("Sensitivity") +
    theme(axis.text = element_text(size = 16),
          axis.title = (element_text(size = 18)), plot.title = element_text(size = 20),
          legend.text = element_text(size = 16), legend.title = element_text(size = 18))
}

table.auc = function(all.auc){
  c( auc(all.auc$merge[[1]][[1]]), auc(all.auc$merge[[2]][[1]]), auc(all.auc$merge[[3]][[1]]),
     auc(all.auc$merge[[4]]), auc(all.auc$merge[[5]]))
}

table.acc = function(all.res){
  c(sapply(all.res[1:4], function(x){
    get.optim.accuracy(x$res$y, x$res$p.stack)
  }), get.optim.accuracy(all.res[[5]]$y, all.res[[5]]$p.glmer))
}

table.sen = function(all.res){
  cbind(sapply(all.res[1:4], function(x){
    get.optim.ss(x$res$y, x$res$p.stack)
  }), get.optim.ss(all.res[[5]]$y, all.res[[5]]$p.glmer))
}

get.optim.accuracy = function( y, p.pred ){
  y = as.logical(y)
  p.prop = mean(y==1)
  roc.res = roc(y, p.pred)
  opt.threshold = roc.res$thresholds[order((1 - roc.res$sensitivities)^2 + (1 - roc.res$specificities)^2)[1]]
  y.pred = 1*(p.pred>=opt.threshold)
  mean(y==y.pred)
}

get.optim.ss = function( y, p.pred ){
  y = as.logical(y)
  roc.res = roc(y, p.pred)
  opt.idx = order((1 - roc.res$sensitivities)^2 + (1 - roc.res$specificities)^2)[1]
  c(roc.res$sensitivities[opt.idx], roc.res$specificities[opt.idx])
}

get.ind.metric = function(all.res, emotion, metric = c("MSE", "AUC")){
  # browser()
  combined_diff = sapply(all.res, function(res){
    if(metric == "AUC"){
      # need to check which ind works better
      metric_all = sapply(res$ind, function(x){
        ens_res = auc(x$data$y, x$data$p.stack, direction = "<")
        ind_res_all = apply( x$data %>% select(contains("p.in") | contains("V1")), 2, function(z){
          auc(x$data$y, z, direction = "<")
        })
        c( ens_res,
           max(ind_res_all) )
      })
    }else{
      metric_all = sapply(res$ind, function(x){
        ens_res = mean((as.numeric(x$data$y) - 1 - x$data$p.stack)^2)
        ind_res_all = apply( x$data %>% select(contains("p.in") | contains("V1")), 2, function(z){
          mean((as.numeric(x$data$y) - 1 - z)^2)
        })
        c( ens_res, min(ind_res_all) )
      })
    }
    metric_all[1,] - metric_all[2,]
  })
  data.frame( ID = sapply(all.res[[1]]$ind, function(x) x$data$ID[1]),
              Diff = c(combined_diff),
              Emotion = emotion,
              Method = rep(c("ENet", "SVM", "RF"), each = nrow(combined_diff)) )
}

get.ind = function(all.res, emotion, metric = c("MSE", "AUC")){
  # browser()
  all_metrics = lapply(all.res, function(res){
    if(metric == "AUC"){
      # need to check which ind works better
      metric_all = sapply(res$ind, function(x){
        ens_res = auc(x$data$y, x$data$p.stack, direction = "<")
        ind_res_all = apply( x$data %>% select(contains("p.in") | contains("V1")), 2, function(z){
          auc(x$data$y, z, direction = "<")
        })
        c( ens_res,
           max(ind_res_all) )
      })
    }else{
      metric_all = lapply(res$ind, function(x){
        ens_res = mean((as.numeric(x$data$y) - 1 - x$data$p.stack)^2)
        ind_res_all = apply( x$data %>% select(contains("p.in") | contains("V1")), 2, function(z){
          mean((as.numeric(x$data$y) - 1 - z)^2)
        })
        c( ens_res, min(ind_res_all) )
      })
    }
  })
  all_metrics = do.call(cbind, all_metrics)
  data.frame( ID = sapply(all.res[[1]]$ind, function(x) x$data$ID[1]),
              ens =  all_metrics[1,],
              idio = all_metrics[2,],
              Emotion = emotion,
              Method = rep(c("ENet", "SVM", "RF"), each = ncol(all_metrics)/3) )
}

