# install packages if not already installed
if (!require("tidyverse")) install.packages("tidyverse")
if (!require("caret")) install.packages("caret")
if (!require("purrr")) install.packages("purrr")
if (!require("broom")) install.packages("broom")
if (!require("Metrics")) install.packages("Metrics")

# load necessary libraries
library(tidyverse)
library(caret)
library(purrr)
library(broom)
library(Metrics)
library(doParallel)
# function to perform rolling window linear regression on each column in df.Y using df.X
library(rollRegres)
library(data.table)
one_eval=function(train.X,train.Y,test.X,test.Y,features){
  if (length(features) == 0 || (length(features)==1 && (features == 'EMPTY'))) {
    train.predictions <- rep(0, length(train.Y))
    test.predictions <- rep(0, length(test.Y))
  } else {
    model <- lm(train.Y ~ .-1, data = data.frame(train.X))
    
    # Make predictions
    train.predictions <- predict(model, newdata = data.frame(train.X))
    test.predictions <- predict(model, newdata = data.frame(test.X))
  }
  
  # Compute in-sample and out-sample RMSE and MAE
  in_sample_rmse <- sqrt(mean((train.Y - train.predictions)^2))
  in_sample_mae <- mean(abs(train.Y - train.predictions))
  
  
  out_sample_resid <- test.Y - test.predictions
  
  results_ <- tibble(Time = t + window_size, 
                     Outcome = names(df.Y)[i], 
                     In_Sample_RMSE = in_sample_rmse, 
                     In_Sample_MAE = in_sample_mae, 
                     Out_Sample_RE = out_sample_resid,
                     Out_Sample_Y = test.Y,
                     Out_Sample_Yhat=test.predictions
  )
  return(results_)
}

regress_and_evaluate_rolling <- function(df.X, df.Y, window_size, list_of_features) {
  T_obs <- nrow(df.X)
  
  if (T_obs != nrow(df.Y) || T_obs <= window_size || length(list_of_features) != (T_obs - window_size+1)) {
    stop("Check dimensions of data and window size")
  }
  
  # Initialize results data frame
  results <- tibble(
    Time = integer(),
    Outcome = character(),
    In_Sample_RMSE = double(),
    In_Sample_MAE = double(),
    Out_Sample_RE = double(),
    Out_Sample_Y = double(),
    Out_Sample_Yhat = double()
  )
  
  # Iterate over the windows
  for (t in seq_len(T_obs - window_size)) {
    # print(t)
    features <- list_of_features[[t]]
    
    # Prepare the rolling window train and test datasets
    if (length(features) == 0 || (length(features)==1 && (features == 'EMPTY'))) {
      train.X <- 0
      test.X <-0 
      
    }else{
      train.X <- df.X[t:(t + window_size - 1), features, drop = FALSE]
      test.X <- df.X[(t + window_size), features, drop = FALSE]
      
    }
    # Iterate over the columns of df.Y
    all_y_result=foreach(i= seq_len(ncol(df.Y)),.packages=c('tidyverse','data.table'),.combine='rbind' ) %dopar% {
      train.Y <- df.Y[t:(t + window_size - 1), i]
      test.Y <- df.Y[(t + window_size), i]
      if (length(features) == 0 || (length(features)==1 && (features == 'EMPTY'))) {
        train.predictions <- rep(0, length(train.Y))
        test.predictions <- rep(0, length(test.Y))
      } else {
        model <- lm(train.Y ~ .-1, data = data.frame(train.X))
        
        # Make predictions
        train.predictions <- predict(model, newdata = data.frame(train.X))
        test.predictions <- predict(model, newdata = data.frame(test.X))
      }
      
      # Compute in-sample and out-sample RMSE and MAE
      in_sample_rmse <- sqrt(mean((train.Y - train.predictions)^2))
      in_sample_mae <- abs(mean(train.Y - train.predictions))
      
      
      out_sample_resid <- test.Y - test.predictions
      
      results_ <- tibble(Time = t + window_size, 
                         Outcome = names(df.Y)[i], 
                         In_Sample_RMSE = in_sample_rmse, 
                         In_Sample_MAE = in_sample_mae, 
                         Out_Sample_RE = out_sample_resid,
                         Out_Sample_Y = test.Y,
                         Out_Sample_Yhat=test.predictions
      )
      return(results_)
    }
    results=results %>% add_row(all_y_result)
    
  }
  
  return(results)
}

DS_location = 'data/DS_data.csv'
HL_location = 'data/HML_data.csv'
DS_returns=read.csv(DS_location,row.names = 'yyyymm')
t_orig=nrow(DS_returns)
HL_factors = read.csv(HL_location,row.names = 'yyyymm')

eval_methods = c("B_OLS","N_LASSO",'P_LASSO','B_POSI','P_POSI','EMPTY','FF5')

var_names = colnames(HL_factors)

ptm = proc.time() 
window_size = 60
saved_selection_results_loc="step2/output/variable_selection_details_list_min_se_HL_60.RData"
saved_selection_results_loc="step2/output/variable_selection_details_list_min_se_HL_60_wFF3prior.RData"

variable_selection_result_list <- readRDS(saved_selection_results_loc)


# Initialize an empty list to store results for all FWER thresholds
month_names = rownames(DS_returns)[1:length(variable_selection_result_list)]
all_results_list = list()
i_window=1
for (i_window in 1:length(variable_selection_result_list)) {
  method = 'P_POSI'
  temp = variable_selection_result_list[[i_window]]
  temp = temp[temp$rho_inv.N.p_1<1,]
  temp$Starting_month = month_names[i_window]
  
  # Append result to the list
  all_results_list[[i_window]] = temp
}


# Combine all results into one big data.table
final_result_dt = rbindlist(all_results_list)
row.names(final_result_dt)

final_result_dt[,'Selected_by_0.05':=1.0*(final_result_dt$rho_inv.N.p_1<=0.05)]
final_result_dt[,'Selected_by_0.01':=1.0*(final_result_dt$rho_inv.N.p_1<=0.01)]
final_result_dt[,'Selected_by_0.1':=1.0*(final_result_dt$rho_inv.N.p_1<=0.1)]

count_table=final_result_dt[,.('Count_FWER_0.05'=sum(Selected_by_0.05),
                   'Count_FWER_0.01'=sum(Selected_by_0.01),
                   'Count_FWER_0.1'=sum(Selected_by_0.1)),by=c('factor')]

fwrite(final_result_dt,'step3/output/selection_detail_windowSize_60_wFF3prior.csv')
fwrite(count_table,'step3/output/selection_countl_windowSize_60_wFF3prior.csv')
length(unique(final_result_dt$Starting_month))

print(final_result_dt)
cl <- makeCluster(12)
registerDoParallel(cl)
ptm = proc.time() 
window_size = 60

saved_selection_results_loc=paste("step2/output/variable_selection_result_list_min_se_HL_",window_size,".RData",sep='')

variable_selection_result_list <- readRDS(saved_selection_results_loc)
for(FWER_threshold in FWER_threshold_grid){
variable_selection_result_list_restructured = list()
for (i_window in 1:length(variable_selection_result_list)){
    for (method in c('P_POSI')){ # eval_methods
    
    temp=variable_selection_result_list[[i_window]]
    
    removed_columns=c('FWER','METHOD')
    feature_names =setdiff( colnames(variable_selection_result_list[[1]]), removed_columns)
    
    
    
    selected_var_names=var_names[!is.na(temp[METHOD==method & FWER==FWER_threshold,.SD,.SDcols=feature_names])]
    
    if(i_window==1){
        variable_selection_result_list_restructured[[method]]=list()
    }
    variable_selection_result_list_restructured[[method]][[i_window]]=selected_var_names
    }
    variable_selection_result_list_restructured[['EMPTY']][[i_window]]='EMPTY'
    variable_selection_result_list_restructured[['FF5']][[i_window]]=c("SMB","HML","RMW","CMA" )
}
evaluation_list=list()
for (method_ in eval_methods) {
    evaluation_list[[method_]]=regress_and_evaluate_rolling (df.X=HL_factors, df.Y=DS_returns,
                                                            window_size=window_size , list_of_features=variable_selection_result_list_restructured[[method_]])
    time_elapsed = round( proc.time() - ptm,2)

    print(paste('Done:',time_elapsed[3],'s'))

}
saved_selection_results_loc=paste("step3/output/selected_model_performance_list_min_se_HL_",window_size,"FWER_",FWER_threshold,".RData",sep='')

saveRDS(evaluation_list, file=saved_selection_results_loc)
}


stopCluster(cl)
time_elapsed = round( proc.time() - ptm,2)

print(paste('Step 3 completed with ',time_elapsed[3],'s'))

