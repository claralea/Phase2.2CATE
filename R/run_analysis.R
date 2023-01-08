#' run CATE analysis
#'
#' @param input.path string with the path to the folder containing phase 2.2 data
#' @param output.path string with the path to the folder to save the results
#' @param siteid string
#' @return dataframe containing the results for all outcomes of interest
#' @import SuperLearner data.table dplyr tidyr partykit rpart
#' @importFrom stats  as.formula binomial cov gaussian glm na.omit predict
#' @importFrom utils write.csv
#' @export



# input.path = dir.repo = dir.data = output.path = "../"
# siteid="VA"
load("sysdata.rda") # Harrison's mapping from concept_code to PheCode
mapping.old=icd10.phecode.map # rename icd10.phecode.map stored in sysdata.rda
mapping = mapping.old #In case where site has their own mapping such as VA. By default, use Harrison's mapping
variant.dist = fread(paste0(input.path, '/4CE_Face_Variant_Distr.csv')) # variant distribution from 2020-04 to 2022-12
variant.dist = variant.dist %>% dplyr::group_by(week) %>% dplyr::slice_max(freq) # dominant variant
colnames(variant.dist) = c('month', 'variant', 'freq')

run_analysis = function(input.path, output.path, siteid="VA"){
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
  X_df <- data.frame(X)
  X_df[,1:ncol(X_df)] <- lapply(X_df[,1:ncol(X_df)], factor)
  rank_x <- sapply(droplevels(X_df), nlevels) > 1 
  
  #### Ising Model using "IsingFit" package #### 
  Ising.model <- IsingFit(X[,rank_x], family='binomial', plot=FALSE, progressbar=TRUE) 
  
  #### Bayesian Network using "bnlearn" package ####
  structure.bic <- hc(data.frame(X_df[,rank_x]), score = "bic") 
  bn.mod.bic <- bn.fit(structure.bic, data = X_df[,rank_x])
  structure.bde <- hc(data.frame(X_df[,rank_x]), score = "bde") 
  bn.mod.bde <- bn.fit(structure.bde, data = X_df[,rank_x])
  
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
                 mean.X = colMeans(X),
                 rank.x = rank_x,
                 site.size = dim(X)[1])
  
  X = as.matrix(data.input$X.train)
  A = as.matrix(data.input$A.train)
  
  for(out in colnames(data.input$Y.train)){
    Y = as.matrix(data.frame(data.input$Y.train)[, out])
    print(paste0("Running analysis for outcome: ", out))
    
    res_y = site_estimation3(X, A, Y)
    
    res_all[[out]] = list(mu1.fit = res_y$mu1_est,
                          mu0.fit = res_y$mu0_est)
  }
  res_all[["train_rank_x"]] = res_y$test_rank_x
  
  print(paste0("Saving the results in ", output.path))
  save(res_all, file = paste0(output.path, '/CATE_results_', siteid,".rda"))
}