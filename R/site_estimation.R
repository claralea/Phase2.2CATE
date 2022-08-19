# Site Estimation -------------------------------------------------------------
# Description: Estimate CATE using only data in each site
# Date: June, 2022
# -----------------------------------------------------------------------------
# Input:
# X - n * p matrix: binary covariates
# A - n * 1 vector: treatment assignments
# Y - n * 1 vector: outcomes
# threshold - integer: the occurrence of xn cannot be shared if below threshold
# -----------------------------------------------------------------------------
# Output: tree/linear model tree learned from different methods
# pseudo_lmtree - lmtree object: fitted linear model tree based on pseudo outcomes
# muhat_lmtree  - lmtree object: fitted linear model tree based on mu1 - mu0
# IPW_lmtree    - lmtree object: fitted linear model tree based on mu1 - mu0 weighted by 1/ps
# xn            - 2^p * 1 vector: the number of occurrence
# test_rank_x   - p * 1 vector: indicator of full rank of covariate matrix
# -----------------------------------------------------------------------------
#' Estimate CATE using only data in each site
#'
#' @param X matrix containing all binary covariates
#' @param A matrix containing binary vaccine information
#' @param Y matrix containing binary outcome
#' @param threshold integer to only count covariates with enough occurences
#' @import SuperLearner
#' @import rpart
#' @import rattle
#' @import partykit
#' @import reshape2
#' @import ggpubr

site_estimation <- function(X, A, Y, threshold=0){
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
  ####### Doubly robust pseudo outcome regression following Kennedy, 2020 ######

  ## estimate nuisance functions
  print("Estimating nuisance functions...")
  pi_fit <- glm(A ~ ., data = data.frame(cbind(X, A)), family = "binomial")
  mu1_fit <- SuperLearner(Y =Y[A==1], X = X[A==1, test_rank_x], family = binomial(),
                          SL.library = c("SL.mean", "SL.glmnet", "SL.rpartPrune"))
  mu0_fit <- SuperLearner(Y =Y[A==0], X = X[A==0, test_rank_x], family = binomial(),
                          SL.library = c("SL.mean", "SL.glmnet", "SL.rpartPrune"))
  print("Estimating nuisance functions is finished!")

  ## construct pseudo outcomes
  pi_est <- predict(pi_fit, X[,test_rank_x], type = "response")
  mu1_est <- predict(mu1_fit, X[,test_rank_x], onlySL = TRUE)$pred
  mu0_est <- predict(mu0_fit, X[,test_rank_x], onlySL = TRUE)$pred
  pseudo <- ((A-pi_est)/(pi_est*(1-pi_est)))*(Y-A*mu1_est-(1-A)*mu0_est) + mu1_est-mu0_est
  df_pseudo <- data.frame(cbind(X,pseudo))
  colnames(df_pseudo) = c(colnames(X),"pseudo")

  ## fit a local linear tree based on pseudo outcomes
  pseudo_lmtree_fit <- lmtree(pseudo ~ . | .,data = df_pseudo[,c(test_rank_x, TRUE)])
  print("Doubly robust pseudo-outcome regression is finished!")

  ##############################################################################
  ############# Simply using muhat_1 - muhat_0  #############
  # muhat_est <- mu1_est - mu0_est
  muhat_est <- A*(Y - mu0_est)+(1-A)*(mu1_est - Y)

  # fit a local linear tree based on muhat
  df_muhat = data.frame(cbind(X, muhat_est))
  colnames(df_muhat) = c(colnames(X),"muhat_est")
  muhat_lmtree_fit <- lmtree(muhat_est ~ . | .,data = df_muhat[,c(test_rank_x, TRUE)])
  print("muhat regression is finished!")

  ############# IPW muhat_1 - muhat_0  #############
  IPW.pseudo <- ((A-pi_est)/(pi_est*(1-pi_est)))*Y
  df_IPW = data.frame(cbind(X, IPW.pseudo))
  colnames(df_IPW) = c(colnames(X),"IPW.pseudo")
  IPW_lmtree_fit <- lmtree(IPW.pseudo ~ . | .,data = df_IPW[,c(test_rank_x, TRUE)])
  print("IPW regression is finished!")

  res <- list(pseudo_lmtree = pseudo_lmtree_fit,
              muhat_lmtree = muhat_lmtree_fit,
              IPW_lmtree = IPW_lmtree_fit,
              xn = xn,
              test_rank_x = test_rank_x)
  return(res)

}
