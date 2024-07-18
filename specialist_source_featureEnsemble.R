enet.wrapper = function(X, y){
  if(sum(y==1)<=3){
    y.tmp = sample(c(1,1,rep(0, length(y)-2)))
    glm(y~., family = "binomial", data = data.frame(y = y.tmp, X))
  }
  glm(y~., family = "binomial", data = data.frame(y = y, X))
}

svm.wrapper = function(X, y){
  data_fit = data.frame(y = y, X)
  svm( y ~ ., data = data_fit, type = "C-classification", probability = T)
}

rf.wrapper = function(X, y){
  data_fit = data.frame(y = y, X)
  randomForest::randomForest(y ~ ., data = data_fit)
}

specialist.wrapper.bothfeature = function(data.ind, outcome, formula.lst, learner.lst, other.models, lambda = 0, nfolds = 10){
  # browser()
  y = as.factor(data.ind[[outcome]])
  folds = createFolds(y, k = nfolds)
  
  L_tot = length(learner.lst)*length(formula.lst)
  
  # create
  p.stacking = do.call(rbind, lapply(folds, function(fold){
    in.models = do.call(c, lapply(formula.lst, function(formula){
      X = model.matrix(formula, data = data.ind)[,-1]
      lapply(learner.lst, function(learner){
        learner(X[-fold,,drop=F], y[-fold])
      })
    }))
    all.models = c(in.models, other.models)
    
    sapply(all.models, function(ind.model){
      if(class(ind.model)[1] == "glm"){
        predict(ind.model, data.ind[fold,,drop = F], type = "response")
      }else if(class(ind.model)[1] == "svm.formula"){
        attributes(predict(ind.model, data.ind[fold,,drop = F], probability = T))$probabilities[,2]
      }else if(class(ind.model)[1] == "randomForest.formula"){
        predict(ind.model, data.ind[fold,,drop = F], type = "prob")[,2]
      }
    })
  }))
  
  y.obs = y[unlist(folds)]
  # cv.res = cv.glmnet(x = p.stacking, y = y.obs, family = "binomial", alpha = 1)
  # glm.res = glmnet(x = p.stacking, y = y.obs, family = "binomial", alpha = 1, lambda = cv.res$lambda.min)
  # c(glm.res$a0, as.numeric(glm.res$beta))
  # nnls(A = p.stacking, b = as.numeric(y.obs)-1)$x
  w.simplex( p.stacking, as.numeric(y.obs) - 1, L_tot, lambda = lambda )
}


penalized.stack.wrapper.bothfeature = function(data.ls, outcome, formula.lst, learner.lst, n.cpus = 4){
  library(snowfall)
  sfInit(parallel = T, cpus = n.cpus)
  # browser()
  # L is the number of algorithms in consideration
  L = length(learner.lst)
  # M is the number of formula
  M = length(formula.lst)
  L_tot = L*M
  
  all.models = do.call(c, lapply(data.ls, function(data.ind){
    # loop over formulas
    do.call(c, lapply(formula.lst, function(formula){
      # get rid of the intercept
      X = model.matrix(formula, data = data.ind)[,-1]
      y = as.factor(data.ind[[outcome]])
      model.res = lapply(learner.lst, function(learner){
        learner(X, y)
      })
      model.res
    }))
  }))
  
  # specialist training
  sfExport("data.ls", "outcome", "formula.lst", "learner.lst", "all.models", "L_tot")
  sfSource("specialist_source.R")
  sfSource("specialist_source_both.R")
  sfLibrary(caret)
  sfLibrary(alabama)
  sfLibrary(e1071)
  sfLibrary(randomForest)
  specialist.auc.all = sfLapply(seq_along(data.ls), function(idx){
    print(idx)
    # X = model.matrix(formula, data = data.ls[[idx]])[,-1]
    y = as.factor(data.ls[[idx]][[outcome]])
    
    folds = createFolds(y, k = 10)
    remove.idx = (idx-1)*L_tot + (1:L_tot)
    
    res.stack = lapply(folds, function(fold){
      print(fold)
      in.models = do.call(c, lapply(formula.lst, function(formula){
        X = model.matrix(formula, data = data.ls[[idx]])[,-1]
        lapply(learner.lst, function(learner){
          learner(X[-fold,,drop=F], y[-fold])
        })
      }))
      w.spe = specialist.wrapper.bothfeature(data.ls[[idx]], outcome, formula.lst, learner.lst, 
                                 all.models[-remove.idx], nfolds = 10)
      
      all.pred = sapply(c(in.models, all.models[-remove.idx]), function(ind.model){
        if(class(ind.model)[1] == "glm"){
          predict(ind.model, data.ls[[idx]][fold,,drop = F], type = "response")
        }else if(class(ind.model)[1] == "svm.formula"){
          attributes(predict(ind.model, data.ls[[idx]][fold,,drop = F], probability = T))$probabilities[,2]
        }else if(class(ind.model)[1] == "randomForest.formula"){
          predict(ind.model, data.ls[[idx]][fold,,drop = F], type = "prob")[,2]
        }
      })
      if(is.null(dim(all.pred))){
        p.in = matrix(all.pred[1:L_tot], nrow = 1)
      }else{
        p.in = all.pred[,1:L_tot,drop = F]
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

analysis.wrapper.both = function(outcome, formula.pca, data.use, n.cpus = 4, cv = c("regular", "ts")){
  # browser()
  df.ls = pca.preproc(data.use, outcome, formula.pca)
  df.combined = lapply(seq_along(data.use), function(i){
    data.use[[i]] = cbind(data.use[[i]], df.ls[[i]] %>% select(starts_with("Comp")))
  })
  names(df.combined) = names(df.ls)
  
  formula.lst = list( as.formula(sprintf("%s ~ %s", outcome, as.character(formula.pca)[2])),
                      as.formula(paste(outcome,paste(grep("Comp", colnames(df.ls[[1]]), value = T), collapse = "+"),
                                       sep = "~")))
  
  data.ls = df.combined[sapply(df.combined, function(x) sum(x[[outcome]]))>4]
  
  if(cv == "regular"){
    print("glmer")
    glmer.res = glmer.wrapper(df.ls[sapply(df.ls, function(x) sum(x[[outcome]]))>4], outcome)
    print("enet")
    enet.res = penalized.stack.wrapper.bothfeature(data.ls, outcome, formula.lst, c(enet.wrapper), n.cpus = n.cpus)
    print("svm")
    svm.res = penalized.stack.wrapper.bothfeature(data.ls, outcome, formula.lst, c(svm.wrapper), n.cpus = n.cpus)
    print("rf")
    rf.res = penalized.stack.wrapper.bothfeature(data.ls, outcome, formula.lst, c(rf.wrapper), n.cpus = n.cpus)
    print("DE")
    de.res = penalized.stack.wrapper.bothfeature(data.ls, outcome, formula.lst, c(enet.wrapper, svm.wrapper, rf.wrapper), n.cpus = n.cpus)
  }else if(cv == "ts"){
    glmer.res = NA
    enet.res = penalized.stack.tsCV.wrapper(data.ls, formula.data, outcome, c(enet.wrapper))
    svm.res = penalized.stack.tsCV.wrapper(data.ls, formula.data, outcome, c(svm.wrapper))
    rf.res = penalized.stack.tsCV.wrapper(data.ls, formula.data, outcome, c(rf.wrapper))
    de.res = penalized.stack.tsCV.wrapper(data.ls, formula.data, outcome, c(enet.wrapper, svm.wrapper, rf.wrapper))
  }
  
  list(enet.res, svm.res, rf.res, de.res, glmer.res)
}
