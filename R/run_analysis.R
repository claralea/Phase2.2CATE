#' run CATE analysis
#'
#' @param input.path string with the path to the folder containing phase 2.2 data
#' @param output.path string with the path to the folder to save the results
#' @param siteid string
#' @return dataframe containing the results for all outcomes of interest
#' @export

# setwd("/Users/")

library(SuperLearner)
library(ranger)
library(data.table)
library(dplyr)
library(tidyr)
library(partykit)
library(rpart)

input.path = dir.repo = dir.data = output.path = "../" 
siteid="VA"
if (siteid=="VA"){
  threshold <- 11
} else{
  threshold <- 0 #occurrence of xn cannot be shared if below threshold 
  
}
source("data_input.R")
source("site_estimation.R")
load("sysdata.rda")

run_analysis = function(input.path, output.path, siteid="BI",threshold=0){
  
  "Reading the data..."
  data.input = create.table(input.path)
  
  # count the empirical distribution of X
  X = as.matrix(rbind(data.input$X.train, data.input$X.test))
  p <- dim(X)[2]
  X_full <- as.matrix(expand.grid(rep(list(0:1), p)))
  xn <- as.vector(table(rbind(X, X_full) %*% 2^c(0:(p-1))) - 1)
  xn[xn < threshold] <- 0 # only share cell with enough occurrence
  xn <- as.integer(xn)
  
  X = as.matrix(data.input$X.train)
  A = as.matrix(data.input$A.train)
  
  X_full_df <- data.frame(X_full)
  X_full_df[,1:p] <- lapply(X_full_df[,1:p], factor)
  colnames(X_full_df) = colnames(X)
  X_full_eval <- X_full_df[xn > 0,]
  
  res_all = list(xn=xn)
  
  for(out in colnames(data.input$Y.train)){
    Y = as.matrix(data.frame(data.input$Y.train)[, out])
    print(paste0("Running analysis for outcome: ", out))
    
    res_y = site_estimation(X, A, Y)
    
    ITE_pred <- suppressWarnings(predict(res_y$DR_fit, newdata = X_full_eval[, res_y$test_rank_x]))
    res_all[[out]] = ITE_pred
  }
  
  print(paste0("Saving the results in ", output.path))
  save(res_all, file = paste0(output.path, '/CATE_results_', siteid,".rda"))
}
