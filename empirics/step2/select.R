select_for_all_methods.common_factor_model=function(OLS_log_t_pval_vec,this_PoSI_output,methods_,FWER_threshold,d,N){
  each_window_variable_selection_ = data.frame(matrix(data=NA,nrow = length(methods_),ncol = dim(this_PoSI_output$posi_log_pval_matrix)[1]))
  var_names_features  = row.names(this_PoSI_output$posi_log_pval_matrix)
  colnames(each_window_variable_selection_)=var_names_features
  
  N_OLS = var_names_features[which(OLS_log_t_pval_vec<=log(FWER_threshold))]
  each_window_variable_selection_[1,N_OLS]=1
  B_OLS = var_names_features[which(OLS_log_t_pval_vec+log(d)+log(N)<=log(FWER_threshold))]
  each_window_variable_selection_[2,B_OLS]=1
  
  
  naive_t_p_val=this_PoSI_output$naive_t_p_val_matrix
  
  
  LASSO_min_p =exp(apply(as.matrix(naive_t_p_val), 1, FUN = min) )
  N_LASSO=var_names_features[which(LASSO_min_p<=FWER_threshold)]
  B_LASSO=var_names_features[which(LASSO_min_p*d*N<=FWER_threshold)]
  CB_LASSO=var_names_features[which(LASSO_min_p*sum(this_PoSI_output$posi_log_pval_matrix!=Inf)<=FWER_threshold)]
  each_window_variable_selection_[3,N_LASSO]=1
  each_window_variable_selection_[4,B_LASSO]=1
  each_window_variable_selection_[5,CB_LASSO]=1
  
  
  P_LASSO_df=panel_unordered(posi_log_pval_matrix = naive_t_p_val)
  P_LASSO_df$factor=row.names(P_LASSO_df)
  
  P_LASSO = P_LASSO_df$factor[P_LASSO_df$rho_inv.N.p_1<=FWER_threshold]
  each_window_variable_selection_[6,P_LASSO]=1
  
  
  P_POSI_df=panel_unordered(posi_log_pval_matrix = as.matrix(this_PoSI_output$posi_log_pval_matrix))
  conditional.bonf = sum(this_PoSI_output$posi_log_pval_matrix!=Inf)
  
  P_POSI_df$factor=row.names(P_POSI_df)
  P_POSI = P_POSI_df$factor[P_POSI_df$rho_inv.N.p_1<=FWER_threshold]
  PoSI_min_log_p = apply(as.matrix(this_PoSI_output$posi_log_pval_matrix), 1, FUN = min)
  
  B_POSI = var_names_features[which(PoSI_min_log_p+log(d)+log(N)<=log(FWER_threshold))]
  CB_POSI = var_names_features[which(PoSI_min_log_p+log(sum(this_PoSI_output$posi_log_pval_matrix!=Inf))<=log(FWER_threshold))]
  each_window_variable_selection_[7,B_POSI]=1
  each_window_variable_selection_[8,CB_POSI]=1
  each_window_variable_selection_[9,P_POSI]=1
  each_window_variable_selection_$FWER = FWER_threshold
  each_window_variable_selection_$METHOD =   methods_
  
  return(as.data.table(each_window_variable_selection_))
}

report_selection_details=function(this_PoSI_output,d,N){
  
  
  P_POSI_df=panel_unordered(posi_log_pval_matrix = as.matrix(this_PoSI_output$posi_log_pval_matrix))
  conditional.bonf = sum(this_PoSI_output$posi_log_pval_matrix!=Inf)
  P_POSI_df$condi.bonf = conditional.bonf
  P_POSI_df$factor=row.names(P_POSI_df)

  return(as.data.table(P_POSI_df))
}

panel_nested_reject = function(posi_log_pval_matrix,significance_level){
  # Calculate Y_check    
  Y_check_matrix=-posi_log_pval_matrix
  Y_check_matrix[is.na(Y_check_matrix)]=-Inf
  
  
  Simes_rejection_table=as.data.frame(matrix(data=F,ncol = 2,nrow = nrow(posi_log_pval_matrix)))
  row.names(Simes_rejection_table)=row.names(posi_log_pval_matrix)
  colnames(Simes_rejection_table)=c('q_check','N_check')
  
  # N_check_k iteration
  
  for(i_row in nrow(Simes_rejection_table):1){
    this_pval_row=Y_check_matrix[i_row,]
    meaningful_ind=!(this_pval_row==-Inf)
    meaningful_ind[is.na(meaningful_ind)]=F
    meaningful_ind=which(meaningful_ind) # script_K_i
    
    if(i_row==nrow(Simes_rejection_table)){
      
      script_K_check = meaningful_ind
    }else{
      
      script_K_check = union(meaningful_ind,script_K_check)
    }
    N_check = sum(Y_check_matrix[i_row:nrow(Y_check_matrix),script_K_check]!=-Inf)
    
    Simes_rejection_table[i_row,2]=N_check
    
  }
  
  # Z_check iteration
  Z_check_cumulative_switch=F
  for(i_row in nrow(Simes_rejection_table):1){
    this_pval_row=Y_check_matrix[i_row,]
    meaningful_ind=!(this_pval_row==-Inf)
    meaningful_ind[is.na(meaningful_ind)]=F
    meaningful_ind=which(meaningful_ind) # script_K_i
    
    if(sum(meaningful_ind)==0){
      Simes_rejection_table[i_row,1]=1
      this_row_contrib_to_Z_check=0
      next
    }
    
    if(!Z_check_cumulative_switch){
      if(i_row==nrow(Simes_rejection_table)){
        denominator=Simes_rejection_table[1,2]
      }else{
        denominator=Simes_rejection_table[1,2]-Simes_rejection_table[i_row+1,2]  
      }
      this_row_contrib_to_Z_check=sum(this_pval_row[meaningful_ind]/denominator)
      Z_check_cumulative=this_row_contrib_to_Z_check
      Z_check_cumulative_switch=T
      
    }else{
      denominator=Simes_rejection_table[1,2]-Simes_rejection_table[i_row+1,2]
      this_row_contrib_to_Z_check=sum(this_pval_row[meaningful_ind]/denominator)
      Z_check_cumulative=this_row_contrib_to_Z_check+Z_check_cumulative
      
      
    }
    q_check = exp(-Z_check_cumulative)
    
    
    Simes_rejection_table[i_row,1]=q_check
    
    
  }
  
  q_check_threshold=Simes_rejection_table$N_check/(nrow(Y_check_matrix)*ncol(Y_check_matrix))*significance_level
  # RHS=(1:nrow(Y_check_matrix))/nrow(Y_check_matrix)*significance_level
  Rejected_indices=which(Simes_rejection_table$q_check<=q_check_threshold)
  Simes_rejection_table$q_check_threshold=q_check_threshold
  Simes_rejection_table$Reject=factor(
    ifelse(Simes_rejection_table$q_check<=q_check_threshold,'Reject','NotReject'), levels = c("Reject", "NotReject"))
  return(Simes_rejection_table)
}


library(data.table)
library(corpcor)
library(glmnet)
library(pracma)
library(stringr)
function_file = 'functions.R'
source(file = function_file)
library(doParallel)  # will load parallel, foreach, and iterators
DS_location = 'data/DS_data.csv'
HL_location = 'data/HML_data.csv'
DS_returns=read.csv(DS_location,row.names = 'yyyymm')
t_orig=nrow(DS_returns)
HL_factors = read.csv(HL_location,row.names = 'yyyymm')

methods = c('N_OLS','B_OLS','N_LASSO','B_POSI','P_POSI')
debug_=F
cl <- makeCluster(12)
registerDoParallel(cl)

ptm = proc.time() 

window_size = 60
regression_results=c()
# read in the data from each file
file_list=list()
variable_selection_result_details = list()

sliding_window = seq_len(nrow(DS_returns) - window_size + 1)
# specified results location
save_location = "Regression_results_list_min_se_HL_Win_"

saved_results_location="step1/outuput/Regression_results_list_min_se_HL_Win_60_wFF3prior.RData"
curr_rds=readRDS(saved_results_location)
for (temp in curr_rds){
  if(is.null(temp)){next}
  # break
  OLS_log_t_pval_matrix=temp$OLS
  this_PoSI_output=temp$LASSO
  i_window=temp$i_window
  OLS_log_t_pval_vec=apply(X=OLS_log_t_pval_matrix, MARGIN=2, FUN=min) 
  
  variable_selection_result_details[[i_window]]=  report_selection_details(this_PoSI_output,d,N)

  time_elapsed = round( proc.time() - ptm,2)
  print(paste('Window #',i_window,' for width',window_size,':',time_elapsed[3],'s'))
}
selection_output = step2
to_save_results_location="step2/outuput/Regression_results_list_min_se_HL_Win_60_wFF3prior.RData"
saveRDS(variable_selection_result_details, 
        file=to_save_results_location)


FWER_threshold_grid=c(0.1,0.05,0.01)
sliding_window = seq_len(nrow(DS_returns) - window_size + 1)

saved_results_location='step2/output/Regression_results_list_min_se_HL_Win_60.RData'
curr_rds=readRDS(saved_results_location)
for (temp in curr_rds){
  if(is.null(temp)){next}
  OLS_log_t_pval_matrix=temp$OLS
  this_PoSI_output=temp$LASSO
  i_window=temp$i_window
  OLS_log_t_pval_vec=apply(X=OLS_log_t_pval_matrix, MARGIN=2, FUN=min) 
  matrix_of_sums <- foreach(FWER_threshold =FWER_threshold_grid,.packages='data.table') %dopar% {
    select_for_all_methods.common_factor_model(OLS_log_t_pval_vec,this_PoSI_output, methods,FWER_threshold,d = ncol(HL_factors),N=ncol(DS_returns))
  }

  time_elapsed = round( proc.time() - ptm,2)
  print(paste('Window #',i_window,' for width',window_size,':',time_elapsed[3],'s'))
}

saveRDS(variable_selection_result_list, file="step2/output/variable_selection_result_list_min_se_HL_60.RData")

stopCluster(cl)
print('Step 2 variable selection completed.')
