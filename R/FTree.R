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

#' @import SuperLearner
#' @import rpart
#' @import rattle
#' @import partykit
#' @import reshape2
#' @import ggpubr
#'
FTree <- function(p, res_all_select){
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
  for (m in levels(federated_combo$Method)){
    federated_lmtree_fit <- lmtree(value ~ . | .,data = federated_combo[federated_combo$Method == m, - (p+1)])
    federated_lmtree_est <- suppressWarnings(predict(federated_lmtree_fit, newdata = X_full_df))
    res_federated <- rbind(res_federated,
                           cbind(rep(m, 2^p), federated_lmtree_est))
  }
  res_federated <- data.frame(res_federated)
  colnames(res_federated) <- c("Method", "Est")
  res_federated$Est <- as.numeric(res_federated$Est)
  return(res_federated)
}
