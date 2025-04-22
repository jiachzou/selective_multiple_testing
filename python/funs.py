
# Created by Jiacheng Zou @ Sept 2, 2022

import pandas as pd
import numpy as np

def panel_unordered(log_pval_matrix):
    '''
    # --------------------------------------------------
    # Input: a matrix of log p values. dimension J by N.
    # -- J: number of features. N: number of units.
    # --------------------------------------------------
    # Output: a pd.DataFrame of sorted multiple testing evidence.
    # -- rho_inv.N.p_1:     each feature's aggregated p values, adjusted for panel multiplicity
    # -- rho_inv.N:         each feature's panel multiplicity, adjusted for rho
    # -- N:                 each feature's panel multiplicity, adjusted for rho
    # -- p_1:               each feature's smallest p value
    # -- rho:               panel cohesive coefficient
    # --------------------------------------------------
    '''
    if type(log_pval_matrix)!=pd.core.frame.DataFrame:
        log_pval_matrix= pd.DataFrame(log_pval_matrix)

    log_pval_matrix = log_pval_matrix.fillna(float('inf'))
    
    result_df = pd.DataFrame(False, index=log_pval_matrix.index, 
                                                    columns=['rho_inv.N.p_1','rho_inv.N','N','p_1'])

    K_set = (log_pval_matrix != float('inf')).sum(axis=1)
    M_set = (log_pval_matrix != float('inf')).sum(axis=0)
    N_vec = []

    for i_row in range(log_pval_matrix.shape[0]):
        this_pval_row = log_pval_matrix.iloc[i_row,:]
        meaningful_ind = ~(this_pval_row == float('inf'))
        N_vec.append(M_set[meaningful_ind].sum())
    N_vec = np.array(N_vec)
    rho_inv = (K_set[N_vec>0] / N_vec[N_vec>0]).sum()
    rho = 1 / rho_inv
    
    for i_row in range(log_pval_matrix.shape[0]):
        this_pval_row = log_pval_matrix.iloc[i_row,:]
        meaningful_ind = ~(this_pval_row == float('inf'))

        if meaningful_ind.sum() > 0:
            my_df = N_vec[i_row]
            p_1 = this_pval_row[meaningful_ind].min()

            bonf_level = np.exp(p_1) * my_df * rho_inv
            result_df.iloc[i_row,:] = [bonf_level,my_df*rho_inv,my_df,p_1]
        else:
            result_df.iloc[i_row,:] = [np.nan,0,0,np.nan]

    result_df['p_1'] = np.exp(result_df['p_1'].astype(float))
    result_df = result_df.sort_values('rho_inv.N.p_1')
    result_df['rho'] = rho

    return result_df

def panel_unordered_singleFWER(log_pval_matrix, FWER=0.05):
    '''
    # --------------------------------------------------
    # Input: 
    # - log_pval_matrix :: a matrix of log p values. dimension J by N
    # -- J: number of features. N: number of units.
    # - FWER:: the FWER threshold.
    # --------------------------------------------------
    # Output: 
    # - selected_features:: the names of the selected features, subject to the FWER specified.
    # - result_df:: a pd.DataFrame of sorted multiple testing evidence. See more in panel_unordered's documentation
    # --------------------------------------------------
    '''
    if type(log_pval_matrix)!=pd.core.frame.DataFrame:
        log_pval_matrix= pd.DataFrame(log_pval_matrix)

    log_pval_matrix = log_pval_matrix.fillna(float('inf'))
    
    result_df = pd.DataFrame(False, index=log_pval_matrix.index, 
                                                    columns=['rho_inv.N.p_1','rho_inv.N','N','p_1'])

    K_set = (log_pval_matrix != float('inf')).sum(axis=1)
    M_set = (log_pval_matrix != float('inf')).sum(axis=0)
    N_vec = []

    for i_row in range(log_pval_matrix.shape[0]):
        this_pval_row = log_pval_matrix.iloc[i_row,:]
        meaningful_ind = ~(this_pval_row == float('inf'))
        N_vec.append(M_set[meaningful_ind].sum())
    N_vec = np.array(N_vec)
    rho_inv = (K_set[N_vec>0] / N_vec[N_vec>0]).sum()
    rho = 1 / rho_inv
    
    for i_row in range(log_pval_matrix.shape[0]):
        this_pval_row = log_pval_matrix.iloc[i_row,:]
        meaningful_ind = ~(this_pval_row == float('inf'))

        if meaningful_ind.sum() > 0:
            my_df = N_vec[i_row]
            p_1 = this_pval_row[meaningful_ind].min()

            bonf_level = np.exp(p_1) * my_df * rho_inv
            result_df.iloc[i_row,:] = [bonf_level,my_df*rho_inv,my_df,p_1]
        else:
            result_df.iloc[i_row,:] = [np.nan,0,0,np.nan]

    result_df['p_1'] = np.exp(result_df['p_1'].astype(float))
    result_df = result_df.sort_values('rho_inv.N.p_1')
    result_df['rho'] = rho
    selected_features =np.sort(result_df.index[result_df['rho_inv.N.p_1']<=FWER]).tolist()

    return selected_features,result_df
