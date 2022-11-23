#' run CATE analysis
#'
#' @param input.path string with the path to the folder containing phase 2.2 data
#' @param output.path string with the path to the folder to save the results
#' @param siteid string
#' @return dataframe containing the results for all outcomes of interest
#' @export
run_analysis = function(input.path, output.path, siteid){
  
  "Reading the data..."
  data.input = create.table(input.path)
  X = as.matrix(data.input$X)
  A = as.matrix(data.input$A)
  
  set.seed(123)
  train_ind <- sample(1:nrow(A), size=floor(0.2 * nrow(A)))
  X <- X[train_ind,]
  A <- A[train_ind]
  
  p=ncol(X)
  X_full <- as.matrix(expand.grid(rep(list(0:1), p)))
  X_full_df <- data.frame(X_full)
  X_full_df[,1:p] <- lapply(X_full_df[,1:p], factor)
  colnames(X_full_df) = colnames(X)
  
  res_all = NULL
  
  for(i in seq(1:12)){
    Y = as.matrix(data.input$Y.all[, ..i])
    Y <- Y[train_ind]
    print(paste0("Running analysis for outcome: ", colnames(data.input$Y.all[, ..i])))
    
    res_y = site_estimation(X, A, Y)
    
    pseudo_lmtree <- suppressWarnings(predict(res_y$DR_lmtree, newdata = X_full_df[, res_y$test_rank_x]))
    muhat_lmtree <- suppressWarnings(predict(res_y$muhat_lmtree, newdata = X_full_df[, res_y$test_rank_x]))
    IPW_lmtree <- suppressWarnings(predict(res_y$IPW_lmtree, newdata = X_full_df[,res_y$test_rank_x]))
    
    res_y_df = cbind(data.frame(pseudo_lmtree, muhat_lmtree, IPW_lmtree, xn=res_y$xn),
                     outcome = rep(colnames(data.input$Y.all[, ..i])),
                     siteid = siteid )
    
    res_all = rbind(res_all, res_y_df)
    
  }
  
  print(paste0("Saving the results in ", output.path))
  
  write.csv(res_all, file = paste0(output.path, '/CATE_results_', siteid, '.csv'), row.names = F)
  
}
