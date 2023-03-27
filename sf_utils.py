#define helper functions
import pandas as pd
import numpy as np
from scipy.stats import zscore
import igraph as ig
import scipy
import upsetplot
import scanpy as sc

import adjustText
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec


from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
import shap
import pickle
from math import comb

import cassiopeia as cas
from cassiopeia.critique.critique_utilities import get_outgroup
import itertools

def get_informative_triplets(static_tree):
    '''
    Create informative triplets from a casseopiea tree. An informative triplet is one with defined out group.
    Inputs: casseopiea tree object
    Outputs: list of leaf triplets along with outgroup
    '''
    
    #omits all triplets with 3 root connected leaves
    root_connected_leaves = list(filter(lambda x:static_tree.parent(x) == static_tree.root, static_tree.leaves))
    non_root_connected_leaves = list(filter(lambda x:static_tree.parent(x) != static_tree.root, static_tree.leaves))
    
    print("Omitting all triplets with 3 root connected leaves")
    triplet_list = []
    triplet_log = []

    #possible triplets: 1 root + 2 non-root; 2 root + 1 non-root; 3 non-root

    # 1 root + 2 non-root
    triplets = list(itertools.combinations(non_root_connected_leaves, 2))
    for i in itertools.product(root_connected_leaves, triplets):
        triplet_curr = [i[0],i[1][0],i[1][1]]
        triplet_curr = np.sort(triplet_curr)
        triplet_curr_str = "/".join(triplet_curr)

        #find out grp
        out_grp = get_outgroup(static_tree, tuple(triplet_curr))

        #if out group is not null add to list of triplets
        if(out_grp!="None"):
            triplet_list.append([triplet_curr, out_grp])
    print("Done case 1: 1 root + 2 non-root leaves")

    # 2 root + 1 non-root
    triplets = list(itertools.combinations(root_connected_leaves, 2))
    for i in itertools.product(non_root_connected_leaves, triplets):
        triplet_curr = [i[0],i[1][0],i[1][1]]
        triplet_curr = np.sort(triplet_curr)
        triplet_curr_str = "/".join(triplet_curr)

        #find out grp
        out_grp = get_outgroup(static_tree, tuple(triplet_curr))

        #if out group is not null add to list of triplets
        if(out_grp!="None"):
            triplet_list.append([triplet_curr, out_grp])
    print("Done case 2: 2 root + 1 non-root leaves")

    # 3 non-root      
    for i in itertools.combinations(non_root_connected_leaves, 3):
        triplet_curr = i
        triplet_curr = np.sort(triplet_curr)
        triplet_curr_str = "/".join(triplet_curr)

        #find out grp
        out_grp = get_outgroup(static_tree, tuple(triplet_curr))

        #if out group is not null add to list of triplets
        if(out_grp!="None"):
            triplet_list.append([triplet_curr, out_grp])
    print("Done case 3: 3 non-root leaves")
    
    return(triplet_list)


def compare_triplet_accuracy(triplet_list, new_tree):
    '''
    Compare accuracy of outgroup calling between an existing tree and a new tree. Takes output of 'get_information_triplets' as input.
    Inputs: leaf triplet-outgroup list from 'get_information_triplets'; Casseopiea object of tree to test.
    Outputs: percentage of correct outgroup matches between new and existing tree
    '''
    match_list = []
    for i in triplet_list:
        ter_outgrp = get_outgroup(new_tree, i[0])

        if(ter_outgrp!='None'):
            match_list.append([i[1] == ter_outgrp])
    match_list = np.array(match_list).flatten()

    return(sum(match_list)/len(match_list))



### NON-TREE FUNCTIONS ###
def shap_loader(file_addr, mode='rb'):
    with open(file_addr, mode) as f:
        return(pickle.load(f))

def get_top_n_features(shap_object, num_ft):
    return(np.array(shap_object.feature_names)[np.argsort(np.abs(shap_object.values).mean(axis=0))[::-1]][:num_ft])

def corr_shap(shap_obj_curr, df, top_vals, return_fig, figsize=(4,10), ax_user=None,color_corr=False,
palette=sns.color_palette('Spectral', n_colors=101).as_hex()):
    '''
    arguments in order:
    shap object,
    cell x ft matrix of feature values,
    number of values to plot, 
    boolean for returning figure,
    outputs: shap-ft correlation table, figure (optional)
    '''
    
    df_shap = shap_obj_curr.values
    top_idx = np.argsort(np.abs(df_shap).sum(axis=0))[::-1][:top_vals]
    df_shap = df_shap[:,top_idx]
    df = df.iloc[:,top_idx].copy()
    #import matplotlib as plt
    # Make a copy of the input data
    shap_v = pd.DataFrame(df_shap)
    feature_list = df.columns
    shap_v.columns = feature_list
    df_v = df.copy().reset_index().drop('cell.bc',axis=1)
    
    # Determine the correlation in order to plot with different colors
    corr_list = list()
    for i in feature_list:
        b = np.corrcoef(shap_v[i],df_v[i])[1][0]
        corr_list.append(b)
    corr_df = pd.concat([pd.Series(feature_list),pd.Series(corr_list)],axis=1).fillna(0)
    # Make a data frame. Column 1 is the feature, and Column 2 is the correlation coefficient
    corr_df.columns  = ['Variable','Corr']
    corr_df['Sign'] = np.where(corr_df['Corr']>0,sns.color_palette('PiYG').as_hex()[-1],sns.color_palette('PiYG').as_hex()[0])
    
    if(color_corr):
        cols = np.array(list(palette))
        corr_df['Sign'] = cols[np.floor(50*(1 + corr_df['Corr']/np.max(np.abs(corr_df['Corr'])))).astype(int)]

    # Plot it
    shap_abs = np.abs(shap_v)
    k=pd.DataFrame(shap_abs.mean()).reset_index()
    k.columns = ['Variable','SHAP_abs']
    k2 = k.merge(corr_df,left_on = 'Variable',right_on='Variable',how='inner')
    k2 = k2.sort_values(by='SHAP_abs',ascending = True)
    colorlist = k2['Sign']
    
    fig = plt.figure(figsize=figsize)
    if(ax_user is not None):
        ax=ax_user
    else:
        ax=plt.axes()
    # k2.plot.barh(x='Variable',y='SHAP_abs',color = colorlist, figsize=(5,6),legend=False, ax=ax)
    plt.barh(k2['Variable'].values, k2['SHAP_abs'].values,color=k2['Sign'], )
    ax.set_xlabel("Absolute SHAP value")
    
    if(return_fig):
        return(fig, corr_df)
    
    return(corr_df)

#simulator for homoplasy in celltagging data
def find_homoplasy(n_cells, moi, barcode_abundance, ct_min=2, ct_max=25, n_iters=100):
    print("Simulating celltag data")
    duplication_rate_list = []
    
    for iteration_curr in range(n_iters):
        if(iteration_curr%10==0):
            print(f"Iteration: {iteration_curr}")
    
        #simulate n_cells
        cells = np.random.poisson(moi, n_cells)

        #filter out cells outside ct_min and ct_max constraints
        filter_1 = cells >= ct_min
        filter_2 = cells <= ct_max

        filter_final = filter_1 & filter_2
        cells = cells[filter_final]
        # print(len(cells))

        #assign tags to each cell
        celltag_sigs = {}
        seen_lens = set()

        # print("Generating CellTag signatures")
        for i in range(len(cells)):
            #simulate a celltag signature based on barcode abundance
            celltag_sig_curr = np.sort(np.random.choice(barcode_abundance.index, size=cells[i], p=barcode_abundance.iloc[:,0]))
            celltag_sig_curr = "/".join(celltag_sig_curr)

            #add sorted signature to dict - based on length of signature
            if(cells[i] in celltag_sigs.keys()):
                celltag_sigs[cells[i]].append(celltag_sig_curr)
            else:
                celltag_sigs[cells[i]] = []
                celltag_sigs[cells[i]].append(celltag_sig_curr)

        #check duplication rates
        # print("Checking for duplicates")
        net_dup_pairs = 0
        for i in celltag_sigs.keys():
            df_curr = pd.DataFrame(celltag_sigs[i])
            df_count = df_curr.value_counts()
            df_count = df_count[df_count > 1].copy()

            if(len(df_count) > 0):
                dup_pairs = sum([comb(x,2) for x in df_count.values])
                net_dup_pairs += dup_pairs
                # print(f"Duplicate pairs ({i} celltags/cell): {dup_pairs}")
        duplication_rate_list.append(net_dup_pairs/comb(len(cells),2))
    
    print("Finished!")
    return(duplication_rate_list)


def plot_trajectories_large(adata,emb_name, rep_clones, de_clones, clone_table, rep_col,
                            de_col,assay=['RNA'], print_d3=False, ret_fig=False):
    
    #define cell groups in each trajectory
    rep_d3 = clone_table[(clone_table['clone_id'].isin(rep_clones)) & (clone_table['day'] == 'D3') & (clone_table['assay'].isin(assay))]['cell.bc']
    de_d3 = clone_table[(clone_table['clone_id'].isin(de_clones)) & (clone_table['day'] == 'D3') & (clone_table['assay'].isin(assay))]['cell.bc']

    rep_d12 = clone_table[(clone_table['clone_id'].isin(rep_clones)) & (clone_table['day'] == 'D12') & (clone_table['assay'].isin(assay))]['cell.bc']
    de_d12 = clone_table[(clone_table['clone_id'].isin(de_clones)) & (clone_table['day'] == 'D12') & (clone_table['assay'].isin(assay))]['cell.bc']

    rep_d21 = clone_table[(clone_table['clone_id'].isin(rep_clones)) & (clone_table['day'] == 'D21') & (clone_table['assay'].isin(assay))]['cell.bc']
    de_d21 = clone_table[(clone_table['clone_id'].isin(de_clones)) & (clone_table['day'] == 'D21') & (clone_table['assay'].isin(assay))]['cell.bc']
    

    #get coordinates for kdeplots
    de_ctr_d3_all = pd.DataFrame(adata[adata.obs_names.isin(de_d3),].obsm[emb_name])
    rep_ctr_d3_all = pd.DataFrame(adata[adata.obs_names.isin(rep_d3),].obsm[emb_name])

    de_ctr_d12_all = pd.DataFrame(adata[adata.obs_names.isin(de_d12),].obsm[emb_name])
    rep_ctr_d12_all = pd.DataFrame(adata[adata.obs_names.isin(rep_d12),].obsm[emb_name])

    de_ctr_d21_all = pd.DataFrame(adata[adata.obs_names.isin(de_d21),].obsm[emb_name])
    rep_ctr_d21_all = pd.DataFrame(adata[adata.obs_names.isin(rep_d21),].obsm[emb_name])
    
    if(print_d3):
        print(f'number of reprogramming cells on d3: {len(rep_ctr_d3_all)}')
        print(f'number of dead-end cells on d3: {len(de_ctr_d3_all)}')
    


    
    #init figure
    fig = plt.figure(figsize=(12,8))
    gs = GridSpec(2,3)
    ax0 = plt.subplot(gs[0])
    ax1 = plt.subplot(gs[1])
    ax2 = plt.subplot(gs[2])
    
    ax3 = plt.subplot(gs[3])
    ax4 = plt.subplot(gs[4])
    ax5 = plt.subplot(gs[5])
    
    #plot reprogramming
    #contour plots of reprogramming atac clones
    sc.pl.embedding(adata, s=20, ax=ax0, show=False,basis=emb_name, frameon=False)
    sns.kdeplot(rep_ctr_d3_all[0], rep_ctr_d3_all[1], ax=ax0, color=rep_col)
    sns.kdeplot(rep_ctr_d3_all[0], rep_ctr_d3_all[1], ax=ax0, color=rep_col, shade=True, alpha=0.6)
    ax0.set_title('reprogramming d3')

    sc.pl.embedding(adata, s=20, ax=ax1, show=False,basis=emb_name, frameon=False)
    sns.kdeplot(rep_ctr_d12_all[0], rep_ctr_d12_all[1], ax=ax1, color=rep_col)
    sns.kdeplot(rep_ctr_d12_all[0], rep_ctr_d12_all[1], ax=ax1, color=rep_col, shade=True, alpha=0.6)
    ax1.set_title('reprogramming d12')

    sc.pl.embedding(adata, s=20, ax=ax2, show=False,basis=emb_name, frameon=False)
    sns.kdeplot(rep_ctr_d21_all[0], rep_ctr_d21_all[1], ax=ax2, color=rep_col)
    sns.kdeplot(rep_ctr_d21_all[0], rep_ctr_d21_all[1], ax=ax2, color=rep_col, shade=True, alpha=0.6)
    ax2.set_title('reprogramming d21')
    
    #plot dead-end
    if(len(de_ctr_d3_all)>5):
        sc.pl.embedding(adata, s=20, ax=ax3, show=False,basis=emb_name, frameon=False)
        sns.kdeplot(de_ctr_d3_all[0], de_ctr_d3_all[1], ax=ax3, color=de_col)
        sns.kdeplot(de_ctr_d3_all[0], de_ctr_d3_all[1], ax=ax3, color=de_col, shade=True, alpha=0.5)
        ax3.set_title('dead end d3')

    sc.pl.embedding(adata, s=20, ax=ax4, show=False,basis=emb_name, frameon=False)
    sns.kdeplot(de_ctr_d12_all[0], de_ctr_d12_all[1], ax=ax4, color=de_col)
    sns.kdeplot(de_ctr_d12_all[0], de_ctr_d12_all[1], ax=ax4, color=de_col, shade=True, alpha=0.6)
    ax4.set_title('dead end d12')

    sc.pl.embedding(adata, s=20, ax=ax5, show=False,basis=emb_name, frameon=False)
    sns.kdeplot(de_ctr_d21_all[0], de_ctr_d21_all[1], ax=ax5, color=de_col)
    sns.kdeplot(de_ctr_d21_all[0], de_ctr_d21_all[1], ax=ax5, color=de_col, shade=True, alpha=0.6)
    ax5.set_title('dead end d21')
    
    #set same y and x lims
    # y_limits = ax2.get_ylim()
    # x_limits = ax2.get_xlim()
    # plt.setp(ax0, xlim=x_limits, ylim=y_limits)
    # plt.setp(ax1, xlim=x_limits, ylim=y_limits)
    # plt.setp(ax2, xlim=x_limits, ylim=y_limits)
    # plt.setp(ax3, xlim=x_limits, ylim=y_limits)
    # plt.setp(ax4, xlim=x_limits, ylim=y_limits)
    # plt.setp(ax5, xlim=x_limits, ylim=y_limits)
    
    plt.suptitle(f'{assay} large clones trajectories')
    
    if(ret_fig):
        return(fig)


def get_clone_cell_embed(adata_obj, clone_table, clone_weight = 1):
    '''
    Creates clone cell embedding, returns combined object with new dims and cluster assignments
    clone_table: columns: first column: cell.bc, second column: clone_id, (opt) third column and beyond: metadata to add to final object
    '''


    all_clones = clone_table['clone_id'].unique()
    all_obs = pd.concat((pd.Series(adata_obj.obs_names), pd.Series(all_clones)))
    cell_cols = adata_obj.obs_names

    #create coembed object
    adata_obj_coembed = sc.AnnData(np.zeros((len(all_obs),100)))

    print("Imputing connectivities matrix")
    #create new connectivities based on edge list
    adata_obj_connectivities = adata_obj.obsp['connectivities']

    if(clone_weight >= 1):
        adata_obj_connectivities = adata_obj_connectivities/clone_weight

    #create clone-cell connectivities

    if(clone_weight < 1):
        imputed_connectivities, new_rows, new_cols = table_to_spmtx(clone_table['clone_id'],
                                                                             clone_table['cell.bc'],clone_weight*np.ones(len(clone_table)))
    else:
        imputed_connectivities, new_rows, new_cols = table_to_spmtx(clone_table['clone_id'],
                                                                             clone_table['cell.bc'],np.ones(len(clone_table)))

    #add missing cells (columns)
    missing_cells = np.array([*set.difference(set(cell_cols), set(new_cols))])
    new_cols = np.hstack((new_cols,missing_cells))
    imputed_connectivities = scipy.sparse.hstack((imputed_connectivities,
                                                  scipy.sparse.csr_matrix(np.zeros((len(all_clones),len(missing_cells)))))).tocsr()

    # reorder 
    idx = []
    for i in cell_cols:
        new_idx = np.where(i==new_cols)
        idx.append(new_idx[0][0])

    assert((new_cols[idx] == cell_cols).all())

    imputed_connectivities = imputed_connectivities.tocsr()[:,idx]
    new_cols = cell_cols

    conn1 = scipy.sparse.vstack((adata_obj_connectivities, imputed_connectivities))
    conn2 = scipy.sparse.vstack((imputed_connectivities.transpose(), scipy.sparse.csr_matrix(np.zeros((len(new_rows),len(new_rows))))))
    final_connectivities = scipy.sparse.hstack((conn1,conn2)).tocsr()

    #assert connectivities mtx properties
    assert((final_connectivities - final_connectivities.T).sum() == 0)
    assert((final_connectivities.diagonal()==0).all())

    #add imputed connectivities to coembed object
    adata_obj_coembed.uns['neighbors'] = dict()
    adata_obj_coembed.uns['neighbors']['params'] = adata_obj.uns['neighbors']['params']
    adata_obj_coembed.obsp['connectivities'] = final_connectivities

    all_obs = np.hstack((np.array(new_cols),new_rows))
    adata_obj_coembed.obs_names = all_obs

    
    return(adata_obj_coembed)

#quantile norm. a dataframe, columns are samples, rows are features
def enrich_fn(cells_1, cells_2, data_mtx,col1_id="grp1", col2_id="grp2", sort_by="delta", fdr=0.05, return_sig = False):
    '''
    This function performs differential enrichment of each feature in data_mtx(cell x ft) across 2 lists of cells, cells_1 and cells_2.
    '''
    
    mtx1 = data_mtx[data_mtx.index.isin(cells_1)]
    mtx2 = data_mtx[data_mtx.index.isin(cells_2)]
    
    mtx1_mean = mtx1.mean(axis=0)
    mtx2_mean = mtx2.mean(axis=0)

    #find mean difference
    mean_diff = mtx1_mean - mtx2_mean

    #find M and A for MA plots (add  to prevent 0s in log)
    val1 = (1 + mtx1.sum(axis=0))/len(mtx1)
    val2 = (1 + mtx2.sum(axis=0))/len(mtx2)
    logfc_diff = np.log2(val1/val2)
    
    #log of 1 + (sum of means/2)
    logfcA = np.log2(1 + ((mtx1_mean + mtx2_mean)/2))
    
    score_table = pd.DataFrame([mtx1_mean,mtx2_mean, mean_diff, logfc_diff, logfcA], index = [col1_id, col2_id, "delta","log2fc","A"]).T

    #p-val and p-adj
    p_list = []
    for i in score_table.index:
        t, p = mannwhitneyu(mtx1[i].values, mtx2[i].values, alternative='two-sided')
        p_list.append(p)

    score_table['p-val'] = p_list

    bh_vals = multipletests(score_table['p-val'], method='fdr_bh', alpha=fdr)
    score_table['p-adj'] = bh_vals[1]
    score_table['p-adj-log'] = -1*np.log10(score_table['p-adj'])
    score_table['is_significant'] = bh_vals[0]


    if(return_sig):
        return(score_table[score_table['is_significant'] == True].sort_values(sort_by, ascending= False))
    return(score_table.sort_values(sort_by, ascending= False))

def quantile_normalize(df):
    """
    input: dataframe with numerical columns
    output: dataframe with quantile normalized values
    """
    df_sorted = pd.DataFrame(np.sort(df.values,
                                     axis=0), 
                             index=df.index, 
                             columns=df.columns)
    df_mean = df_sorted.mean(axis=1)
    df_mean.index = np.arange(1, len(df_mean) + 1)
    df_qn =df.rank(method="min").stack().astype(int).map(df_mean).unstack()
    return(df_qn)

def assign_fate(clone_mtx, fate_col = 'day', fate_key='d5', cell_type_key = 'cell_type2'):
    '''
    Function to assign fate based on  aper clone basis. For each clone, fate is defined as the most common cell type on fate day.
    Inputs:
    clone_mtx: table of a single clone (subsetted from clone table dataframe)
    fate_col: string; column to use for identifying fate cells (default is 'day')
    fate_key: string; value in fate_col that defines fate cells (default is 'd5')
    cell_type_key: string; column to use for cell type (default 'cell_type2')
    '''
    clone_mtx_curr = clone_mtx[clone_mtx[fate_col] == fate_key].copy()
    if(len(clone_mtx_curr) == 0):
        clone_mtx['fate'] = 'no_fate_cells'
        clone_mtx['fate_pct'] = 0
    else:
        clone_mtx['fate'] = clone_mtx_curr[cell_type_key].value_counts(dropna=False).sort_index().idxmax()
        clone_mtx['fate_pct'] = 100*clone_mtx_curr[cell_type_key].value_counts(dropna=False).sort_index().max()/len(clone_mtx_curr)
    return clone_mtx


def plotVariability(feature_mtx, label_pts, label_x, label_y):
    '''
    ArchR inspired fucntion to plot a standard deviation rank plot for a set of features across all cells in a dataset. Features with higher varibility could be biologically more interesting (not always true)
    '''
    return('tbd')
    

#function to create 2 groups of cells, cells belonging to a given fate and cells not belonging to a given fate
def filter_clones(clone_table, fate = ["Mono"], threshold = None, min_cells = None, fate_col = 'day', fate_key='d5', cell_type_key = 'cell_type2'):
    '''itertate through clones, for each state-fate clone, define its fate group based on thresholds,
    if fate grp same as user defined fate, add to fate_table, else add to nonfate_table'''

    if((threshold is None) and (min_cells is None)):
        print("Specify at least threshold or min_cells")
        return
        
    #initialize output tables
    fate_table = pd.DataFrame()
    nonfate_table = pd.DataFrame()
    clones_sf = pd.DataFrame()
    
    for i in clone_table['clone.id'].unique():
        clone_curr = clone_table[clone_table['clone.id'] == i]
        
        sf = clone_curr[fate_col].values
        assay_type = clone_curr['assay'].values
        
        #if clone is state-fate
        if(len(set(sf)) == 2):
            
            #get fates
            fates = clone_curr[clone_curr[fate_col] == fate_key][cell_type_key]
            
            #if both min_cells and threshold specified, use both. If just one of them is specified, use that
            if(min_cells is not None):
                if(threshold is not None):
                    if(((fates.isin(fate)).sum() >= min_cells) | ((fates.isin(fate)).sum()/len(fates) >= threshold)):
                        fate_table = fate_table.append(clone_curr)
                    else:
                        nonfate_table = nonfate_table.append(clone_curr)
                else:
                    if(((fates.isin(fate)).sum() >= min_cells)):
                        fate_table = fate_table.append(clone_curr)
                    else:
                        nonfate_table = nonfate_table.append(clone_curr)
            else:
                if((fates.isin(fate)).sum()/len(fates) >= threshold):
                    fate_table = fate_table.append(clone_curr)
                else:
                    nonfate_table = nonfate_table.append(clone_curr)               
                
            clones_sf = clones_sf.append(clone_curr)
    return fate_table, nonfate_table, clones_sf


def merge_nn(nn_graph, all_cells, cell_list):
    '''
    This function merges a given cell list with its n nearest neighbors
    '''
    nn_cell_set = set(cell_list)
    
    for i in cell_list:
        if i in all_cells:
            cells_curr = all_cells[nn_graph[np.where(all_cells == i)[0][0],:].nonzero()[1]]
            nn_cell_set = nn_cell_set.union(cells_curr)
    return nn_cell_set

#this function tf-idf normalizes (signac implementation) and z-scores a feature x cell matrix (not cell x feature)
def atac_norm(data, tf_idf = True, z_score = True, scale_factor=10000):
    if(tf_idf):
        tf = scale_factor*data/data.sum(axis=0)
        idf = len(data.columns)/data.sum(axis=1)
        data = np.log(1+tf.multiply(idf,axis=0))
    
    if(z_score):
        data_norm = pd.DataFrame(zscore(data,axis=1), index = data.index, columns = data.columns)
        return(data_norm)
    return(data)


def naive_atac_rna_pairing(clone, seed = 100, state_day = 'd2'):
    '''Takes a state fate clone with both RNA and ATAC states and creates pairs naively (random pairing)'''
    
    np.random.seed(seed)
    state_clone = clone[clone['day'] == 'd2'].copy(deep=True)
    atac_cells = state_clone[state_clone['assay'] == 'atac']['cell.barcode'].values
    rna_cells = state_clone[state_clone['assay'] == 'rna']['cell.barcode'].values
    
    atac_len = len(atac_cells)
    rna_len = len(rna_cells)
    
    #pair cells
    if(atac_len == rna_len):
        pair_list = np.vstack((atac_cells, rna_cells))
    if(atac_len < rna_len):
        pair_list = np.hstack((np.vstack((atac_cells, rna_cells[:atac_len])), np.vstack((np.random.choice(atac_cells, size=rna_len-atac_len), rna_cells[atac_len:]))))
    if(rna_len < atac_len):
        pair_list = np.hstack((np.vstack((atac_cells[:rna_len], rna_cells)), np.vstack((atac_cells[rna_len:], np.random.choice(rna_cells, size=atac_len-rna_len)))))
    
    fate_curr = clone['fate'].values[0]
    fate_curr = np.array([fate_curr]*pair_list.shape[1]).reshape(1,-1)  
    return(np.vstack((pair_list, fate_curr)))


def table_to_spmtx(row_data, col_data, count_data):
    '''
    Convert a table to csr sparse matrix
    '''

    cb_u = list(np.sort(np.unique(row_data)))
    celltag_u = list(np.sort(np.unique(col_data)))

    data = count_data.tolist()
    row = pd.Categorical(row_data, categories= cb_u).codes
    col = pd.Categorical(col_data, categories = celltag_u).codes
    celltag_mat = scipy.sparse.csr_matrix((data, (row, col)), shape=(len(cb_u), len(celltag_u)))

    cells = np.array(cb_u)
    celltags = np.array(celltag_u)

    return(celltag_mat, cells, celltags)



