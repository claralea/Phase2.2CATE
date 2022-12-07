# Site Estimation -------------------------------------------------------------
# Description: Estimate CATE using only data in each site
# Date: Oct, 2022
# -----------------------------------------------------------------------------
# Input:
# X - n * p matrix: binary covariates
# A - n * 1 vector: treatment assignments
# Y - n * 1 vector: outcomes
# ps.model - "LR" or "SL": fit propensity score model by logistic regression or superlearner
# or.model - "LR" or "SL": fit outcome regression model by logistic regression or superlearner
# ite.model - "tree" or "lmtree" or "SL": fit CATE by a pruned tree or a linear model tree or superlearner
# -----------------------------------------------------------------------------
# Output: tree/linear model tree learned from different methods
# DR_fit - lmtree object: fitted linear model tree based on DR pseudo outcomes
# test_rank_x   - p * 1 vector: indicator of full rank of covariate matrix
# -----------------------------------------------------------------------------

site_estimation <- function(X, A, Y, ps.model="LR", or.model="SL", ite.model="lmtree"){
  n <- dim(X)[1]
  p <- dim(X)[2]
  
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
    DR_fit <- rpart(DR ~ ., method="anova", data=df_DR[,c(test_rank_x, TRUE)], 
                    control = rpart.control(cp = 0.0001))
    DR_fit<- prune(DR_fit,
                   cp= DR_fit$cptable[which.min(DR_fit$cptable[,"xerror"]),"CP"])
  } else if (ite.model =="lmtree"){
    DR_fit <- lmtree(DR ~ . | .,data = df_DR[,c(test_rank_x, TRUE)])
  } else if (ite.model =="SL"){
    DR_fit <- SuperLearner(Y =DR, X = X[, test_rank_x], family = gaussian(),
                           SL.library = sl_lib)
  }
  print("Doubly robust pseudo-outcome regression is finished!")
  
  res <- list(DR_fit = DR_fit, 
              test_rank_x = test_rank_x)
  return(res)
}

deiden <- function(object, test_rank_x){
  # remove potentially identifiable pts info and reduce output size
  object$data = object$data[1,] 
  object$data[1,] = c(0, rep(factor(1), sum(test_rank_x)))
  object$fitted <- NULL
  rm(list=ls(envir = attr(object$terms, ".Environment")),
     envir = attr(object$terms,".Environment"))
  object$terms <- NULL
  object$info$terms <- NULL
  object$node$info[! names(object$node$info) %in% "object"] <- NULL
  object$node$info$object[c("residuals", "df.residual", "fitted.values","effects","contrasts")] <- NULL
  return(object)
}