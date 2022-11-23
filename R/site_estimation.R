# Site Estimation -------------------------------------------------------------
# Description: Estimate CATE using only data in each site
# Date: Oct, 2022
# -----------------------------------------------------------------------------
# Input:
# X - n * p matrix: binary covariates
# A - n * 1 vector: treatment assignments
# Y - n * 1 vector: outcomes
# threshold - integer: the occurrence of xn cannot be shared if below threshold 
# ps.model - "LR" or "SL": fit propensity score model by logistic regression or superlearner
# or.model - "LR" or "SL": fit outcome regression model by logistic regression or superlearner
# ite.model - "tree" or "lmtree": fit CATE by a pruned tree or a linear model tree
# -----------------------------------------------------------------------------
# Output: tree/linear model tree learned from different methods
# DR_lmtree - lmtree object: fitted linear model tree based on DR pseudo outcomes
# muhat_lmtree  - lmtree object: fitted linear model tree based on mu1 - mu0
# IPW_lmtree    - lmtree object: fitted linear model tree based on mu1 - mu0 weighted by 1/ps
# xn            - 2^p * 1 vector: the number of occurrence
# test_rank_x   - p * 1 vector: indicator of full rank of covariate matrix
# -----------------------------------------------------------------------------

site_estimation <- function(X, A, Y, ps.model="LR", or.model="SL", threshold=0, ite.model="lmtree"){
  n <- dim(X)[1]
  p <- dim(X)[2]
  
  X_full <- as.matrix(expand.grid(rep(list(0:1), p)))
  # count the occurrence of X
  xn <- as.vector(table(rbind(X, X_full) %*% 2^c(0:(p-1))) - 1)
  xn[xn < threshold] <- 0 # only share cell with enough occurrence
  
  X <- data.frame(X)
  X[,1:p] <- lapply(X[,1:p], factor) # binary covariates
  ## check if all columns of X have different values to avoid rank deficient
  test_rank_x <- sapply(droplevels(X), nlevels) > 1 
  ##############################################################################

  ## estimate nuisance functions
  print("Estimating nuisance functions...")
  sl_lib = c("SL.ranger", "SL.glmnet", "SL.rpartPrune", "SL.mean", "SL.glm.interaction")
  if (ps.model == "SL"){
    pi_fit <- SuperLearner(Y =A, X = X[, test_rank_x], family = binomial(),
                            SL.library = sl_lib)
    pi_est <- predict(pi_fit, X[,test_rank_x], onlySL = TRUE)$pred
  } else if (ps.model == "LR"){
    df_ps <- data.frame(cbind(X[,test_rank_x], A))
    colnames(df_ps) = c(colnames(X[,test_rank_x]),"A")
    pi_fit <- glm(A ~ ., data = df_ps, family = "binomial")
    pi_est <- predict(pi_fit, X[,test_rank_x], type = "response")
  }
  if (or.model == "SL"){
    mu1_fit <- SuperLearner(Y =Y[A==1], X = X[A==1, test_rank_x], family = binomial(),
                            SL.library = sl_lib)
    mu0_fit <- SuperLearner(Y =Y[A==0], X = X[A==0, test_rank_x], family = binomial(),
                            SL.library = sl_lib)
    mu1_est <- predict(mu1_fit, X[,test_rank_x], onlySL = TRUE)$pred
    mu0_est <- predict(mu0_fit, X[,test_rank_x], onlySL = TRUE)$pred
  } else if (or.model == "LR"){
    df_mu <- data.frame(cbind(X[, test_rank_x], Y))
    colnames(df_mu) = c(colnames(X[,test_rank_x]),"Y")
    mu1_fit <- glm(Y ~ ., data = df_mu[A==1,], family = "binomial")
    mu0_fit <- glm(Y ~ ., data = df_mu[A==0,], family = "binomial")
    mu1_est <- predict(mu1_fit, X[,test_rank_x], type = "response")
    mu0_est <- predict(mu0_fit, X[,test_rank_x], type = "response")
  }
  print("Estimating nuisance functions is finished!")
  
  ## construct pseudo outcomes
  DR <- ((A-pi_est)/(pi_est*(1-pi_est)))*(Y-A*mu1_est-(1-A)*mu0_est) + mu1_est-mu0_est
  df_DR <- data.frame(cbind(X,DR))
  colnames(df_DR) = c(colnames(X),"DR")
  ## fit a local linear tree or a tree based on pseudo outcomes
  if (ite.model == "tree"){
    DR_lmtree_fit <- rpart(DR ~ ., method="anova", data=df_DR[,c(test_rank_x, TRUE)], 
                               control = rpart.control(cp = 0.0001))
    DR_lmtree_fit<- prune(DR_lmtree_fit,
                              cp= DR_lmtree_fit$cptable[which.min(DR_lmtree_fit$cptable[,"xerror"]),"CP"])
  } else if (ite.model =="lmtree"){
    DR_lmtree_fit <- lmtree(DR ~ . | .,data = df_DR[,c(test_rank_x, TRUE)])
  }
  print("Doubly robust pseudo-outcome regression is finished!")
  
  ##############################################################################
  ############# Simply using muhat_1 - muhat_0  ############# 
  muhat_est <- A*(Y - mu0_est)+(1-A)*(mu1_est - Y)
  # fit a local linear tree based on muhat
  df_muhat = data.frame(cbind(X, muhat_est))
  colnames(df_muhat) = c(colnames(X),"muhat_est")
  if (ite.model == "tree"){
    muhat_lmtree_fit <- rpart(muhat_est ~ ., method="anova", data=df_muhat[,c(test_rank_x, TRUE)], 
                              control = rpart.control(cp = 0.0001))
    muhat_lmtree_fit<- prune(muhat_lmtree_fit,
                             cp= muhat_lmtree_fit$cptable[which.min(muhat_lmtree_fit$cptable[,"xerror"]),"CP"])
  } else if (ite.model =="lmtree"){
    muhat_lmtree_fit <- lmtree(muhat_est ~ . | .,data = df_muhat[,c(test_rank_x, TRUE)])
  }
  print("muhat regression is finished!")
  
  ############# IPW muhat_1 - muhat_0  ############# 
  IPW.pseudo <- ((A-pi_est)/(pi_est*(1-pi_est)))*Y
  df_IPW = data.frame(cbind(X, IPW.pseudo))
  colnames(df_IPW) = c(colnames(X),"IPW.pseudo")
  if (ite.model == "tree"){
    IPW_lmtree_fit <- rpart(IPW.pseudo ~ ., method="anova", data=df_IPW[,c(test_rank_x, TRUE)], 
                            control = rpart.control(cp = 0.0001))
    IPW_lmtree_fit<- prune(IPW_lmtree_fit,
                           cp= IPW_lmtree_fit$cptable[which.min(IPW_lmtree_fit$cptable[,"xerror"]),"CP"])
  } else if (ite.model =="lmtree"){
    IPW_lmtree_fit <- lmtree(IPW.pseudo ~ . | .,data = df_IPW[,c(test_rank_x, TRUE)])
  }
  print("IPW regression is finished!")
  
  res <- list(DR_lmtree = DR_lmtree_fit, 
              muhat_lmtree = muhat_lmtree_fit,
              IPW_lmtree = IPW_lmtree_fit,
              xn = xn,
              test_rank_x = test_rank_x)
  return(res)
}
