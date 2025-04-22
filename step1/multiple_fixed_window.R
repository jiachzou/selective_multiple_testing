
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
library(stringr)

library(pracma)
library(nleqslv)
library(stats)
library(TruncatedNormal)

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
  # X=as.matrix(HL_factors)
  # Y=as.matrix(DS_returns)
  residual_matrix=matrix(ncol=ncol(Y),nrow=t_orig)
  # i_col=1
  for(i_col in 1:length(all_portfolio_names)){
    col = all_portfolio_names[i_col]
    # break  
    
    y = Y[,col]
    y = 5*y
    k_fold = 4
    
    df_guarantee=nrow(X)-4
    # penalty_omega_inv=rep(1,ncol(X))
    cvob0=cv.glmnet(X, y,nfolds = k_fold,penalty.factor = penalty_omega_inv,intercept = F,standardize=F,thresh = 1e-16,parallel =T,pmax=df_guarantee)#,lambda=exp((32:-32)/4)*log(ncol(X))/sqrt(t_obs))
    lasso_screening=glmnet(X,y,lambda = cvob0$lambda,penalty.factor = penalty_omega_inv,intercept = F,standardize = F,thresh = 1e-16,parallel =T,pmax=df_guarantee)
    non_zeros = colSums(lasso_screening$beta!=0)
    rank_defining_lambda_boolean = non_zeros<df_guarantee
    lambda_ = cvob0$lambda[rank_defining_lambda_boolean]
    
    optimal_lambda =  cvob0$lambda.1se 
    optimal_lasso=glmnet(X,y,lambda = optimal_lambda,penalty.factor = penalty_omega_inv,intercept = F)
    optimal_beta=as.vector(optimal_lasso$beta)
    M_set=optimal_beta!=0
    M_cardi=sum(M_set)
    
    if(! (optimal_lambda %in% lambda_)){
      diff_in_lambda_ = lambda_-optimal_lambda
      closest_idx = which(diff_in_lambda_==min(diff_in_lambda_[diff_in_lambda_>0]))
      optimal_lambda=lambda_[closest_idx]
    }
    
    if(M_cardi==0){
      p_raw_vec=rep(Inf,M_cardi)
      next
    }
    y_hat=predict(optimal_lasso,X)
    residuals_lasso=y-y_hat
    
    M_ind=which(M_set)    
    sigma_squared=sum(residuals_lasso^2)/max(1,t_obs-M_cardi)
    
    beta_hat_LASSO=optimal_beta[M_set]
    X_M=X[,M_ind]
    X_M=as.matrix(X_M)
    s_M=sign(beta_hat_LASSO)
    
    X_min_M=X[,which(optimal_beta==0)]
    X_M_Pinv=pinv(X_M) # pseudo-inverse of X_M
    
    P_M = X_M %*%X_M_Pinv
    
    I_min_P_M=diag(dim(P_M)[1])-P_M
    
    y_shifted=y
    
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
    
    
    # Using bounds similar to Lee, Jason D., Dennis L. Sun, Yuekai Sun, and Jonathan E. Taylor. "Exact post-selection inference, with application to the lasso." (2016)
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
      V_minus=max(ratio[negative_ones & (beta_hat_one_step[i_of_M]-ratio>1e-32)],na.rm = T)
      
      V_plus=min(ratio[positive_ones& (ratio-beta_hat_one_step[i_of_M]>1e-32)],na.rm = T)
      
      if(is.na(V_minus) | is.na(V_plus)| is.na(dig_sd[i_of_M])){
        next
      }
      epsilon=1e-16
      if(V_plus-V_minus<=epsilon || dig_sd[i_of_M]<0){
        print('bizzare')
        
        next
        
      }else{
        This_Df=max(1,nrow(X_M)-M_cardi)
        if(beta_hat_one_step[i_of_M]>0){
          if(-V_minus/dig_sd[i_of_M]+V_plus/dig_sd[i_of_M]<=1e-3){
            next
          }
          
          right_tail=ptmvnorm(q=-beta_hat_one_step[i_of_M],mu=0,sigma = dig_sd[i_of_M],
                              lb =-V_plus-epsilon,ub = -V_minus+epsilon,type = 'qmc',
                              log = T,B=20000)
          left_tail=ptmvnorm(q=-beta_hat_one_step[i_of_M],mu=0,sigma = dig_sd[i_of_M],
                             lb =V_minus-epsilon,ub = V_plus+epsilon,type = 'qmc',
                             log = T,B=20000)

        }else{
          if(V_plus/dig_sd[i_of_M]-V_minus/dig_sd[i_of_M]<=1e-3){
            next
          }
          right_tail=ptmvnorm(q=beta_hat_one_step[i_of_M],mu=0,sigma = dig_sd[i_of_M],
                              lb =V_minus-epsilon,ub = V_plus+epsilon,type = 'qmc',
                              log = T,B=20000)
          left_tail=ptmvnorm(q=beta_hat_one_step[i_of_M],mu=0,sigma = dig_sd[i_of_M],
                             lb =-V_plus-epsilon,ub = -V_minus+epsilon,type = 'qmc',
                             log = T,B=20000)

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


fit_single_window=function(X,Y,i_window){
  
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


training_ratio = 1/2

set.seed(123)

# Set the file locations for DS and HML data
DS_location = 'data/DS_data.csv' # response variables: double-sorted portfolio returns from Kenneth French's website
HL_location = 'data/HML_data.csv' # factor variables: high-minus-low risk factors from Hou et al's website
DS_returns=read.csv(DS_location,row.names = 'yyyymm')
t_orig=nrow(DS_returns)
HL_factors = read.csv(HL_location,row.names = 'yyyymm')
factor_choice = 'HL wFF3'

# factor_choice corresponds to the factor priors or using PCs
# - wFF3: use FF3 factors as priors
# - wFF5: use FF5 factors as priors
if(factor_choice %in% c('HL','HL wFF3','HL wFF5')){
  all_factors=HL_factors
} else if(factor_choice=='HL PC'){
  all_factors=HL_PC_factors
} else if(factor_choice=='HL + HL PC'){
  all_factors=cbind(HL_factors,HL_PC_factors)
} 

column_normalize=function(x){x/sd(x)}
calc_pca <- function(i,window_size,covariates) {
  start_row <- i
  end_row <- start_row + window_size - 1
  curr_window_X = covariates[start_row:end_row, ]
  curr_window_X_normalized=apply(curr_window_X,2,column_normalize)
  # Fit PCA using the `prcomp` function
  pca<- prcomp(curr_window_X_normalized, scale = F,rank.	=20,center = F)
  
  return(list(pc=pca$x,rotation=pca$rotation,start_row=start_row,window_size=window_size))
}



all_factors =HL_factors
ff_portfolios=DS_returns

N = ncol(ff_portfolios)
d = ncol(all_factors)

if (factor_choice == 'HL wFF3'){
  penalty_omega_inv = rep(1,d)
  # penalty_omega_inv[(d-4):(d-2)]=0
  penalty_omega_inv[(d-3):(d-2)]=0
  penalty_omega_inv=penalty_omega_inv/sum(penalty_omega_inv)*d
} else if (factor_choice == 'HL wFF5'){
  penalty_omega_inv = rep(1,d)
  penalty_omega_inv[(d-3):(d)]=0
  penalty_omega_inv=penalty_omega_inv/sum(penalty_omega_inv)*d
} else {
  penalty_omega_inv = rep(1,d)
}

print(paste('Factor choice is',factor_choice))
# Creating a 12-core cluster
cl <- makeCluster(12)
registerDoParallel(cl)
# Creating a 60-monthrolling window
window_size=60
sliding_window = seq_len(nrow(ff_portfolios) - window_size + 1)
methods = c('N_OLS','B_OLS','N_LASSO','B_POSI','P_POSI')
ptm=proc.time()
# specifies where to save the results
save_location = "~/Regression_results_list_min_se_HL_Win_"
all_fitted_models <- foreach(i =sliding_window,.packages=c('pracma','glmnet','TruncatedNormal','corpcor'),.inorder = T) %dopar% {
fit_single_window(X = as.matrix(all_factors[i:(i+window_size-1),]),
                    Y = as.matrix(ff_portfolios[i:(i+window_size-1),]),i_window = i)
}
time_elapsed = round( proc.time() - ptm,2)

print(paste('Model fit:',time_elapsed[3],'s -- window',window_size))
saveRDS(all_fitted_models, file=paste(,window_size,".RData",sep=''))

time_elapsed = round( proc.time() - ptm,2)
print(paste('Model fit:',time_elapsed[3],'s -- all'))
stopCluster(cl)
print('Step 1 p-value calculations completed.')