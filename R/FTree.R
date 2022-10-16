# Federated Combination -------------------------------------------------------
# Description: Build a master linear model tree based on estimates from k sites
# Date: 2022
# -----------------------------------------------------------------------------
# Input:
# p:              the number of binary covariates
# res_all_select: predictions of models fitted by each site
# -----------------------------------------------------------------------------
# Output: Estimated CATE
# res_federated: predictions of linear regression trees after federated combination
# -----------------------------------------------------------------------------

FTree <- function(p, res_all_select, prune=NULL, ite.model="lmtree"){
  X_full <- as.matrix(expand.grid(rep(list(0:1), p)))
  X_full_df <- data.frame(X_full)
  X_full_df[,1:p] <- lapply(X_full_df[,1:p], factor)
  
  federated_combo <- res_all_select[res_all_select$xn > 0,] %>%
    group_by(x) %>%
    summarise(pseudo_lmtree = sum(pseudo_lmtree * xn) / sum(xn),
              muhat_lmtree = sum(muhat_lmtree * xn) / sum(xn),
              IPW_lmtree = sum(IPW_lmtree * xn) / sum(xn))
  federated_combo <- cbind(federated_combo, X_full_df[federated_combo$x,])
  federated_combo$x <- NULL
  federated_combo <- melt(federated_combo,
                          id.var = paste("Var", c(1:p), sep=""),
                          variable.name = "Method")
  res_federated <- NULL
  fedtree <- NULL
  for (m in levels(federated_combo$Method)){
    df_fed <- federated_combo[federated_combo$Method == m,]
    df_fed$Method <- NULL
    if (ite.model == "tree"){
      federated_tree_fit <- rpart(value ~ ., method="anova", data=df_fed,
                                  control = rpart.control(cp = 0.0001))
      federated_tree_fit<- prune(federated_tree_fit,
                             cp= federated_tree_fit$cptable[which.min(federated_tree_fit$cptable[,"xerror"]),"CP"])
      federated_tree_est <- predict(federated_tree_fit, newdata = X_full_df)
      res_federated <- rbind(res_federated,
                             cbind(rep(m, 2^p), federated_tree_est))
      fedtree[[m]] <- federated_tree_fit
    } else if (ite.model == "lmtree"){
      federated_lmtree_fit <- lmtree(value ~ . | .,data = df_fed,prune=prune)
      federated_lmtree_est <- suppressWarnings(predict(federated_lmtree_fit, newdata = X_full_df))
      res_federated <- rbind(res_federated,
                             cbind(rep(m, 2^p), federated_lmtree_est))
      fedtree[[m]] <- federated_lmtree_fit
    }
  }
  res_federated <- data.frame(res_federated)
  colnames(res_federated) <- c("Method", "Est")
  res_federated$Est <- as.numeric(res_federated$Est)
  return(list(fedtree.predict=res_federated, fedtree=fedtree))
}

MultiFTree <- function(p, res_all_select, prune = NULL){ 
  X_full <- as.matrix(expand.grid(rep(list(0:1), p)))
  X_full_df <- data.frame(X_full)
  X_full_df[,1:p] <- lapply(X_full_df[,1:p], factor)
  
  federated_combo <- res_all_select[res_all_select$xn > 0,] %>% 
    group_by(x, outcome) %>%
    summarise(pseudo_lmtree = sum(pseudo_lmtree * xn) / sum(xn), 
              muhat_lmtree = sum(muhat_lmtree * xn) / sum(xn),
              IPW_lmtree = sum(IPW_lmtree * xn) / sum(xn))
  federated_combo <- cbind(federated_combo, X_full_df[federated_combo$x,])
  federated_combo$x <- NULL
  federated_combo2 <- melt(federated_combo, 
                           id.var = c(paste("Var", c(1:p), sep=""), "outcome"), 
                           variable.name = "Method")
  federated_combo2$outcome <- as.factor(federated_combo2$outcome)
  
  fedtree.predict <- NULL
  fedtree <- NULL
  for (m in levels(federated_combo2$Method)){
    df_fed <- federated_combo2[federated_combo2$Method == m, ]
    df_fed$Method <- NULL
  
    formula <- as.formula(paste("value ~ ",
                                paste(paste("outcome * Var", c(1:p), sep=""),collapse="+"),
                                "|", paste(paste("Var", c(1:p), sep=""),collapse="+"),
                                sep=""))
    
    federated_lmtree_fit <- lmtree(formula,data = df_fed, prune=prune)
    fedtree[[m]] <- federated_lmtree_fit
    
    federated_lmtree_est <- NULL
    for (outcome in levels(df_fed$outcome)){
      df_pred <- cbind(X_full_df, outcome=rep(outcome, 2^p))
      federated_lmtree_est_o <- suppressWarnings(predict(federated_lmtree_fit, newdata = df_pred))
      federated_lmtree_est <- cbind(federated_lmtree_est, federated_lmtree_est_o)
    }
    fedtree.predict <- rbind(fedtree.predict,
                           cbind(rep(m, 2^p), federated_lmtree_est))
  }
  fedtree.predict <- data.frame(fedtree.predict)
  colnames(fedtree.predict) <- c("Method", levels(df_fed$outcome))
  return(list(fedtree.predict=fedtree.predict, fedtree=fedtree))
}