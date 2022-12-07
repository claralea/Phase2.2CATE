#' run CATE analysis
#'
#' @param input.path string with the path to the folder containing phase 2.2 data
#' @param output.path string with the path to the folder to save the results
#' @param siteid string
#' @return dataframe containing the results for all outcomes of interest
#' @import SuperLearner data.table dplyr tidyr partykit rpart IsingFit bnlearn
#' @importFrom stats  as.formula binomial cov gaussian glm na.omit predict
#' @importFrom utils write.csv
#' @export



# input.path = dir.repo = dir.data = output.path = "../"
# siteid="VA"
# load("sysdata.rda")

run_analysis = function(input.path, output.path, siteid="BI", threshold=0){
  if (siteid=="VA"){
    threshold <- 11
  } else{
    threshold <- 5 #occurrence of xn cannot be shared if below threshold
  }
  "Reading the data..."
  data.input = create.table(input.path)

  # count the empirical distribution of X
  print("Estimating distribution of X...")
  X = as.matrix(rbind(data.input$X.train, data.input$X.test))

  #### Ising Model using "IsingFit" package ####
  Ising.model <- IsingFit(X, family='binomial', plot=FALSE, progressbar=TRUE)

  #### Bayesian Network using "bnlearn" package ####
  X_df <- data.frame(X)
  X_df[,1:ncol(X_df)] <- lapply(X_df[,1:ncol(X_df)], factor)
  structure.bic <- hc(data.frame(X_df), score = "ebic")
  bn.mod.bic <- bn.fit(structure.bic, data = X_df)
  structure.bde <- hc(data.frame(X_df), score = "bde")
  bn.mod.bde <- bn.fit(structure.bde, data = X_df)

  #### count nx > threshold ####
  p <- dim(X)[2]
  X_full <- as.matrix(expand.grid(rep(list(0:1), p)))
  xn <- as.vector(table(rbind(X, X_full) %*% 2^c(0:(p-1))) - 1)
  xn[xn < threshold] <- 0 # only share cell with enough occurrence
  xn <- as.integer(xn)
  res_all = list(Ising.graph=Ising.model$weiadj,
                 Ising.external.field=Ising.model$thresholds,
                 BayesNet.bic=bn.mod.bic,
                 BayesNet.bde=bn.mod.bde,
                 xn=xn,
                 mean.X = colMeans(X))

  X = as.matrix(data.input$X.train)
  A = as.matrix(data.input$A.train)

  for(out in colnames(data.input$Y.train)){
    Y = as.matrix(data.frame(data.input$Y.train)[, out])
    print(paste0("Running analysis for outcome: ", out))

    res_y = site_estimation(X, A, Y, or.model="LR")

    res_all[[out]] = deiden(res_y$DR_fit, res_y$test_rank_x)
  }
  res_all[["test_rank_x"]] = res_y$test_rank_x

  print(paste0("Saving the results in ", output.path))
  save(res_all, file = paste0(output.path, '/CATE_results_', siteid,".rda"))
}
