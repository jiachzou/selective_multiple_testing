
library(data.table)
library(corpcor)
library(glmnet)
library(pracma)
library(stringr)
library(doParallel)  # will load parallel, foreach, and iterators
library(selectiveInference)
library(lmtest)
library(nlme)
library(zeallot)
library(corpcor)
library(stringr)

library(pracma)
library(nleqslv)
library(stats)
library(MASS)

simulate_beta <- function(s, d, N) {
  if (N %% s != 0) {
    stop("N is not divisible by s")
  }
  
  beta <- matrix(0, d, N) # Create a matrix filled with zeros of size d x N
  
  # units_to_sample <- N
  for (i in 1:s) {
    num_selected <- round(N * (s - i + 1) / s)
    sampled_units <- sample(1:N, num_selected, replace = FALSE)
    
    beta[i, sampled_units] <- runif(num_selected,-1,1)
    
  }
  
  return(beta)
}

generate_sample <- function(d, N, T_obs, s_true, noise_std, strong_factor_uniform_bound, xsec_cov, phi_0){
  # Initialize the full_beta matrix
  full_beta <- matrix(0, nrow = d, ncol = N)
  
  # Split the 1:N sequence into s groups
  splitted <- split(1:N, cut(1:N, s_true))
 
  # Iterate through the groups
  for(i_factor in seq_along(splitted)){
    # Note that : i_factor only impacts i_factor / s_true of the N units
    # Concatenate the last i_factor groups
    curr_units <- unlist(splitted[(length(splitted) - i_factor + 1):length(splitted)])
    # Generate random beta values
    full_beta[i_factor, curr_units] <- runif(length(curr_units), min = -strong_factor_uniform_bound, max = strong_factor_uniform_bound)
    
  }
  
  # Define the noise covariance matrix
  noise_cov <- diag(noise_std, N)
  noise_cov[lower.tri(noise_cov)] <- xsec_cov
  noise_cov <- t(noise_cov)
  noise_cov[lower.tri(noise_cov)] <- xsec_cov
  
  # Add a small positive number to the diagonal to ensure positive semi-definiteness
  jitter <- 1e-6
  noise_cov <- noise_cov + diag(jitter, N)
  
  # Initialize the noise matrix
  noises <- matrix(0, nrow = T_obs, ncol = N)
  
  # Generate noise for each time period
  for(t in 1:T_obs){
    noise_t <- MASS::mvrnorm(n = 1, mu = rep(0, N), Sigma = noise_cov)
    if(t == 1){
      noises[t, ] <- noise_t
    } else {
      noises[t, ] <- phi_0 * noises[(t - 1), ] + noise_t
    }
  }
  
  # Generate all covariates
  all_covariates <- matrix(rnorm(T_obs * (d),mean = 1), nrow = T_obs, ncol = (d))
  # all_covariates <- cbind(all_covariates,(all_covariates[,1:s])^3)
  normalized_all_covariates <- apply(all_covariates, 2, function(x) x / sd(x))
  
  # Generate all responses
  all_response <- normalized_all_covariates %*% full_beta + noises
  # Calculate the signal to noise ratio (SNR)
  SNR_list <- apply(full_beta, 2, function(beta) sqrt(sum(beta^2)) / (sqrt(sum(noises^2)) * sqrt(T_obs)))
  SNR <- mean(SNR_list)
  colnames(all_response)=paste('Y',1:N,sep='')
  colnames(all_covariates)=paste('X',1:d,sep='')
  return(list(all_response = all_response, all_covariates = all_covariates, SNR = SNR))
}
panel_unordered = function(posi_log_pval_matrix){
  
  # posi_log_pval_matrix=this_PoSI_output$posi_log_pval_matrix
  
  # posi_log_pval_matrix is J by N
  
  posi_log_pval_matrix[is.na(posi_log_pval_matrix)]=Inf
  HD_cardi=sum(posi_log_pval_matrix!=Inf)
  ############################################################################################
  
  time_series_Bonf_rejection_table=as.data.frame(matrix(data=F,ncol = 4,nrow = nrow(posi_log_pval_matrix)))
  row.names(time_series_Bonf_rejection_table)=row.names(posi_log_pval_matrix)
  colnames(time_series_Bonf_rejection_table)=c('rho_inv.N.p_1','rho_inv.N','N','p_1')
  K_set = rowSums((posi_log_pval_matrix!=Inf)*1.0)
  M_set = colSums((posi_log_pval_matrix!=Inf)*1.0)
  N_vec = vector(length = nrow(posi_log_pval_matrix),mode = 'numeric')
  for(i_row in 1:nrow(time_series_Bonf_rejection_table)){
    this_pval_row=posi_log_pval_matrix[i_row,]
    meaningful_ind=!(this_pval_row==Inf)
    # meaningful_ind[is.na(meaningful_ind)]=F
    N_vec[i_row]=sum(M_set[meaningful_ind])
  }
  rho_inv=sum(K_set[N_vec>0]/N_vec[N_vec>0])
  rho = 1/rho_inv
  for(i_row in 1:nrow(time_series_Bonf_rejection_table)){
    this_pval_row=posi_log_pval_matrix[i_row,]
    meaningful_ind=!(this_pval_row==Inf)
    # meaningful_ind[is.na(meaningful_ind)]=F
    if(sum(meaningful_ind)>0){
      my_df=N_vec[i_row]
      p_1=min(this_pval_row[meaningful_ind],na.rm = T)
      
      bonf_level=exp(p_1)*(my_df)*rho_inv
      # rejection_char=ifelse(bonf_level<significance_level ,'Reject','NotReject')
      time_series_Bonf_rejection_table[i_row,]=c(bonf_level,my_df*rho_inv,my_df,p_1)
      
      
    }else{
      my_df=0
      time_series_Bonf_rejection_table[i_row,]=c(NA,my_df*rho_inv,my_df,NA)
      
      
    }
    
  }
  time_series_Bonf_rejection_table$p_1=exp(as.numeric(time_series_Bonf_rejection_table$p_1))
  
  time_series_Bonf_rejection_table =time_series_Bonf_rejection_table[order(time_series_Bonf_rejection_table$rho_inv.N.p_1),]
  time_series_Bonf_rejection_table$rho = rho
  time_series_Bonf_rejection_table[is.na(time_series_Bonf_rejection_table)]=Inf
  
  time_series_Bonf_rejection_table$HD_cardi=HD_cardi
  return(time_series_Bonf_rejection_table)
}


posi_lognorm_pval_enforce_dimension = function(X, Y, penalty_omega_inv=NA){
  
  
  #######################################
  # Enumerating thru the test portfolios #
  #######################################
  
  t_orig = nrow(X)
  t_obs = t_orig
  N = ncol(Y)
  J = ncol(X)
  
  all_portfolio_names=colnames(Y)
  
  posi_log_pval_matrix=as.data.frame(matrix(data=Inf,ncol = N,nrow = J))
  colnames(posi_log_pval_matrix)=all_portfolio_names
  rownames(posi_log_pval_matrix)=colnames(X)
  
  
  naive_t_p_val=posi_log_pval_matrix
  
  
  beta_hat_table=as.data.frame(matrix(data=NA,ncol = length(all_portfolio_names),nrow = length(colnames(X))))
  colnames(beta_hat_table)=all_portfolio_names
  rownames(beta_hat_table)=colnames(X)
  residual_matrix=matrix(ncol=ncol(Y),nrow=t_orig)
  for(i_col in 1:length(all_portfolio_names)){
    col = all_portfolio_names[i_col]
    
    y = Y[,col]

    X_pseudo_inv = pinv(X)
    beta_hat = X_pseudo_inv %*% y
    residual_from_X_pseudo_inv = y-X%*%beta_hat
    sigma_hat = sd(residual_from_X_pseudo_inv)
    # K.max = as.integer(J/3)
    
    K.max = 50
    the_path =glmnet(X,y,intercept = F,dfmax = K.max)
    id = max(which(the_path$df<=K.max))
    optimal_lambda=the_path$lambda[id]
    lasso_capped =glmnet(X,y,intercept = F,lambda = optimal_lambda    )
    lasso_Lee_etal=lasso_capped
    optimal_beta = as.vector(lasso_Lee_etal$beta)
    
    M_set=optimal_beta!=0
    M_cardi=sum(M_set)
    if(M_cardi==0){
      p_raw_vec=rep(Inf,M_cardi)
      next
    }
    y_hat=predict(lasso_Lee_etal,X)
    residuals_lasso=y-y_hat
    
    M_ind=which(M_set)
    
    sigma_squared=sum(residuals_lasso^2)/max(1,t_obs-M_cardi)
    
    
    beta_hat_LASSO=optimal_beta[M_set]
    X_M=X[,M_ind]
    X_M=as.matrix(X_M)
    s_M=sign(beta_hat_LASSO)
    
    X_min_M=X[,which(optimal_beta==0)]
    X_M_Pinv=pinv(X_M)   # the same as pinv(t(X_M)%*%X_M)%*%t(X_M)
    
    P_M = X_M %*%X_M_Pinv
    
    I_min_P_M=diag(dim(P_M)[1])-P_M
    
    y_shifted=y
    # y_shifted=y
    
    # beta_hat_one_step=X_M_Pinv%*%y
    beta_hat_one_step=X_M_Pinv%*%y_shifted
    beta_hat_table[M_ind,i_col]=beta_hat_one_step
    Fisher_Inv=pinv(t(X_M)%*%X_M)
    covariance_matrix=sigma_squared*Fisher_Inv
    
    dig_var=diag(covariance_matrix)
    dig_sd=sqrt(dig_var)
    OLS_POST_LASSO=lm(y~X_M-1)
    if(dim(X_M)[2]>=dim(X_M)[1]){
      print('Warning! The number of LASSO active set is still larger than OLS')
    }
    
    naive_t_p_val[M_ind,i_col]=log(coef(summary(OLS_POST_LASSO))[,4])
    
    
    # Using Lee, Sun, Sun and Tyalor
    omega_inv_M=penalty_omega_inv[M_ind]
    omega_inv_minus_M=penalty_omega_inv[-M_ind]
    
    component1=t(X_min_M)%*%I_min_P_M
    component2=t(X_min_M)%*%t(X_M_Pinv)%*%(s_M*omega_inv_M)
    if(M_cardi==1){
      A_mat=rbind(1/optimal_lambda*component1,
                  -1/optimal_lambda*component1,
                  -s_M%*%X_M_Pinv)
      
      
      b_vec=rbind(omega_inv_minus_M-component2,
                  omega_inv_minus_M+component2,
                  -optimal_lambda*s_M%*%Fisher_Inv%*%(s_M*omega_inv_M))
    }else{
      A_mat=rbind(1/optimal_lambda*component1,
                  -1/optimal_lambda*component1,
                  -diag(s_M)%*%X_M_Pinv)
      
      
      b_vec=rbind(omega_inv_minus_M-component2,
                  omega_inv_minus_M+component2,
                  -optimal_lambda*diag(s_M)%*%Fisher_Inv%*%(s_M*omega_inv_M))
    }
    
    
    positive_gap = b_vec-A_mat%*%y_shifted
    b_vec[positive_gap<0]=b_vec[positive_gap<0]-min(positive_gap[positive_gap<0])
    p_raw_vec=rep(Inf,M_cardi)
    for(i_of_M in 1:M_cardi){
      this_eta =X_M_Pinv[i_of_M,]
      this_xi = this_eta/sum(this_eta^2)
      xi_eta_outer=this_xi%*%t(this_eta)
      
      this_z=(diag(length(this_eta))-xi_eta_outer)%*%y_shifted
      numerators=b_vec-A_mat%*%this_z
      num_keep_ind=which(numerators>0)
      
      numerators=numerators[num_keep_ind]
      
      denominators=A_mat%*%this_xi
      denominators=denominators[num_keep_ind]
      ratio = numerators/denominators
      
      
      negative_ones = denominators<0
      positive_ones = denominators>0
      # print(paste('Ratio chosen from ',length(ratio)))
      V_minus=max(ratio[negative_ones & (beta_hat_one_step[i_of_M]-ratio>1e-32)],na.rm = T)
      
      V_plus=min(ratio[positive_ones& (ratio-beta_hat_one_step[i_of_M]>1e-32)],na.rm = T)
      
      
      # print(paste(dig_sd[i_of_M],V_minus,V_plus))
      if(is.na(V_minus) | is.na(V_plus)| is.na(dig_sd[i_of_M])){
        next
      }
      epsilon=1e-8
      if(V_plus-V_minus<=epsilon || dig_sd[i_of_M]<0){
        print('bizzare')
        
        next
        
      }else{
        This_Df=max(1,nrow(X_M)-M_cardi)
        if(beta_hat_one_step[i_of_M]>0){
          if(-V_minus/dig_sd[i_of_M]+V_plus/dig_sd[i_of_M]<=1e-3){
            next
          }
          
          right_tail = log(ptruncnorm(q=-beta_hat_one_step[i_of_M], a=-V_plus-epsilon, b=-V_minus+epsilon, mean = 0, sd = dig_sd[i_of_M]))
          left_tail = log(ptruncnorm(q=-beta_hat_one_step[i_of_M], a=V_minus-epsilon, b=V_plus+epsilon, mean = 0, sd = dig_sd[i_of_M]))
          
          
        }else{
          if(V_plus/dig_sd[i_of_M]-V_minus/dig_sd[i_of_M]<=1e-3){
            next
          }
          right_tail = log(ptruncnorm(q=beta_hat_one_step[i_of_M], a=V_minus-epsilon, b=V_plus+epsilon, mean = 0, sd = dig_sd[i_of_M]))
          left_tail = log(ptruncnorm(q=beta_hat_one_step[i_of_M], a=-V_plus-epsilon, b=-V_minus+epsilon, mean = 0, sd = dig_sd[i_of_M]))
          
        }
        
      }
      if(is.na(right_tail) | is.na(left_tail)){
        p_raw=Inf
        next
      } 
      
      if(right_tail==-Inf & left_tail==-Inf){
        p_raw=Inf
      } else if(abs(right_tail-left_tail)>10){
        p_raw=max(right_tail,left_tail)
      }else{
        p_raw=log(exp(right_tail)+exp(left_tail))
      }
      
      p_raw_vec[i_of_M]=p_raw
    }
    posi_log_pval_matrix[M_ind,i_col]=p_raw_vec
    
  }
  
  
  return(list(posi_log_pval_matrix=posi_log_pval_matrix,naive_t_p_val_matrix=naive_t_p_val))
}

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
  
  each_window_variable_selection_[3,N_LASSO]=1
  each_window_variable_selection_[4,B_LASSO]=1
  
  
  
  
  P_POSI_df=panel_unordered(posi_log_pval_matrix = as.matrix(this_PoSI_output$posi_log_pval_matrix))
  P_POSI_df=P_POSI_df[P_POSI_df$p_1!=Inf,]
  
  
  conditional.bonf = sum(this_PoSI_output$posi_log_pval_matrix!=Inf)
  
  P_POSI_df$factor=row.names(P_POSI_df)
  P_POSI = P_POSI_df$factor[P_POSI_df$rho_inv.N.p_1<=FWER_threshold]
  PoSI_min_log_p = apply(as.matrix(this_PoSI_output$posi_log_pval_matrix), 1, FUN = min)
  
  B_POSI = var_names_features[which(PoSI_min_log_p+log(d)+log(N)<=log(FWER_threshold))]
  
  each_window_variable_selection_[5,B_POSI]=1
  each_window_variable_selection_[6,P_POSI]=1
  
  
  each_window_variable_selection_$FWER = FWER_threshold
  each_window_variable_selection_$METHOD =   methods_
  conditional.bonf = sum(this_PoSI_output$posi_log_pval_matrix!=Inf)
  # each_window_variable_selection_$HD=conditional.bonf
  
  n_features_p_posi_more_powerful_ =sum(P_POSI_df$rho_inv.N<conditional.bonf)
  each_window_variable_selection_$Num.rho_inv.N.smaller.than.HD=n_features_p_posi_more_powerful_
  return(as.data.table(each_window_variable_selection_))
  
}

fit_single_window=function(X,Y,i_window){
  
  
  # Fitting OLS
  OLS_log_t_pval_mat = matrix(nrow = ncol(Y),ncol=ncol(X))
  
  for (i_asset in 1:ncol(Y)){
    for (j_factor in 1:ncol(X)){
      curr_LM = lm(Y[,i_asset]~X[,j_factor]-1)  
      OLS_log_t_pval_mat[i_asset,j_factor]=log(coef(summary(curr_LM))[,4])
    }
  }
  
  this_PoSI_output = posi_lognorm_pval_enforce_dimension(X=X,
                                                         Y=Y,
                                                         penalty_omega_inv = penalty_omega_inv)
  print(paste(i_window,'done'))
  return(list(OLS=OLS_log_t_pval_mat,
              LASSO = this_PoSI_output,
              i_window=i_window))
}

save_results_to_csv <- function(output_folder) {
  file_date <- format(Sys.Date(), "%Y%m%d")
  file_name <- paste(output_folder,"grid_search_", file_date, ".csv", sep = "")
  
  # Check if the file already exists, if yes, add an underscore and an enumeration
  enumeration <- 0
  while (file.exists(file_name)) {
    enumeration <- enumeration + 1
    file_name <- paste(output_folder, "grid_search_", file_date, "_", enumeration, ".csv", sep = "")
  }
  
  
  return(file_name)
}

cl <- makeCluster(12)
registerDoParallel(cl)

eval_methods = c('N_OLS','B_OLS','N_LASSO','B_LASSO','B_POSI','P_POSI')

library(tibble)

ptm=proc.time()
set.seed(1)
N_sim=200
FWER_threshold_grid=c(0.1,0.05,0.01)
param_grid <- expand.grid(d = c(100),
                          N = c(120),
                          T_obs = c(300),
                          s = c(10),
                          noise_std = c(2),
                          strong_factor_uniform_bound = c(0.5),
                          xsec_cov =c(0,1),
                          phi_0 = c(0,0.2))
print(nrow(param_grid))


output_folder <- "simulations/outputs/"

output_loc <-save_results_to_csv(output_folder)
average_result_list=list()

# Loop through all rows of the parameter grid
for (i in 1:nrow(param_grid)) {
  # Extract the current parameters
  params <- param_grid[i, ]
  print(paste('Current param.',params))
  # Run the model with these parameters and store the results

  # Assign the values of the parameters to variables with their corresponding names
  list2env(as.list(params), envir = .GlobalEnv)
  
  # Run the model with these parameters
  
  penalty_omega_inv=rep(1,d)
  i_sim=1
  matrix_of_sums <- foreach(i_sim =1:N_sim,.packages = c('data.table','MASS','pracma','glmnet','truncnorm','expm','flare')) %dopar% {
  
    training_sample = 1:as.integer(0.5*T_obs)
    
    one_sample <- generate_sample(d, N, T_obs, s, noise_std, 
                                  strong_factor_uniform_bound, xsec_cov, phi_0)
    
    penalty_omega_inv=rep(1,d)
    temp=fit_single_window(X = one_sample$all_covariates[training_sample,], Y = one_sample$all_response[training_sample,],i_sim)
    OLS_log_t_pval_matrix=temp$OLS
    this_PoSI_output=temp$LASSO
    

    OLS_log_t_pval_vec=apply(X=OLS_log_t_pval_matrix, MARGIN=2, FUN=min) 
    results_list=list()
    for(i_FWER_threshold in 1:length(FWER_threshold_grid)){
      
      FWER_threshold=FWER_threshold_grid[i_FWER_threshold]
      selected_result_ = select_for_all_methods.common_factor_model(OLS_log_t_pval_vec,this_PoSI_output, eval_methods,FWER_threshold,d = d,N=N)
      mat_ = as.matrix(1.0*!is.na(selected_result_))[,1:d]
      
      # Calculate true positives, false positives and false negatives
      tp <- rowSums(mat_[, 1:s])
      fp <- rowSums(mat_[, (s+1):ncol(mat_)])
      fn <- s - tp
      
      # Calculate precision, recall and F1 score
      precision <- tp / (tp + fp)
      recall <- tp / s
      f1 <- 2 * (precision * recall) / (precision + recall)
      var_names_ = colnames(one_sample$all_covariates)
      
      X_train <- one_sample$all_covariates[training_sample,]
      Y_train <- one_sample$all_response[training_sample,]
      X_test <- one_sample$all_covariates[-training_sample,]
      Y_test <- one_sample$all_response[-training_sample,]
      
      unique_methods_=unique(selected_result_$METHOD)
      OOS_R_squared_vec=rep(0,length(unique_methods_))
      OOS_RMSE_vec=rep(0,length(unique_methods_))
      INS_RMSE_vec=rep(0,length(unique_methods_))
      INS_R_squared_vec=OOS_R_squared_vec
      for(i_method in 1:length(unique_methods_)){
        method_=unique_methods_[i_method]
        selected_var_names_=var_names_[!is.na(selected_result_[METHOD==method_,1:d])]
        selected_features_idx <- match(selected_var_names_, var_names_)
        # Predict Y values for the test set using matrix algebra
        if(length(selected_features_idx)==0){
          Y_pred=Y_test*0
          Y_fit=Y_train*0
        }else{
          betas <- pinv(X_train[, selected_features_idx,drop=FALSE]) %*% Y_train
          Y_pred <- X_test[, selected_features_idx,drop=FALSE] %*% betas  
          Y_fit <-X_train[, selected_features_idx,drop=FALSE] %*% betas  
        }
        
        residuals_ =Y_test - Y_pred
        # Calculate R-squared
        RSS <- sum((Y_test - Y_pred)^2)
        TSS <- sum((Y_test - colMeans(Y_test))^2)
        OOS_R_squared_vec[i_method] <- 1 - RSS / TSS
        OOS_RMSE_vec[i_method]=sqrt(RSS/(N*nrow(X_test)))
        
        INS_RSS=sum((Y_train - Y_fit)^2)
        INS_TSS= sum((Y_train - colMeans(Y_fit))^2)
        
        INS_RMSE_vec[i_method]=sqrt(INS_RSS/(N*nrow(X_train)))
        INS_R_squared_vec[i_method]=1 - INS_RSS / INS_TSS
      }
      
      
      
      # selected_result_$Num.rho_inv.N.smaller.than.HD
      # Create dataframe
      results <- data.table(
        False_Selection = fp,
        True_Selection = tp,
        Total_Selection = rowSums(mat_),
        F1_Score = f1,
        SNR=one_sample$SNR,
        INS_R_squared=INS_R_squared_vec,
        INS_RMSE=INS_RMSE_vec,
        OOS_R_squared=OOS_R_squared_vec,
        OOS_RMSE=OOS_RMSE_vec,
        FWER_threshold=FWER_threshold,
        Num.rho_inv.N.smaller.than.HD=selected_result_$Num.rho_inv.N.smaller.than.HD
      )
      results_list[[length(results_list) + 1]] <-results
    }

    return(as.matrix(rbindlist(results_list)))
  }

  
  # Now calculate element-wise averages
  avg_mat <- Reduce("+", matrix_of_sums) / length(matrix_of_sums)
  
  # Convert back to data.table if needed
  avg_dt <- as.data.table(avg_mat)
  avg_dt$method=rep(eval_methods,length(FWER_threshold_grid))
  avg_dt$N_sim=N_sim
  avg_dt$d=d
  avg_dt$N=N
  avg_dt$T_obs=T_obs
  avg_dt$s=s
  avg_dt$strong_factor_uniform_bound=strong_factor_uniform_bound
  avg_dt$noise_std=noise_std
  avg_dt$phi_0=phi_0
  avg_dt$xsec_cov=xsec_cov
  
  time_elapsed = round( proc.time() - ptm,2)
  
  print(paste('Simulation:',time_elapsed[3],'s -- all'))
  average_result_list[[length(average_result_list) + 1]] <-avg_dt
  to_write_ = rbindlist(average_result_list)
  fwrite(to_write_,output_loc)
}


time_elapsed = round( proc.time() - ptm,2)

print(paste('Model fit:',time_elapsed[3],'s -- all'))

stopCluster(cl)

to_write_ = rbindlist(average_result_list)

fwrite(to_write_,output_loc)

print('Simulations completed done')
print(paste('Total time elapsed:',time_elapsed[3],'s -- all'))
