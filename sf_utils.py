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
def find_homoplasy(n_cells, moi, barcode_abundance, ct_min=2, ct_max=25, return_celltag_sigs=False):
    '''
    Function to simulate rate of two unrelated cells getting the same celltag signature
    n_cells: number of starting cells
    moi: expected MOI of transduction
    barcode_abundance: a dataframe with barcode names as index and relative abundance as first column
    ct_min: minimum celltags/cell threshold
    ct_max: maximum celltags/cell threshold
    return_celltag_sigs: whether to return celltag signature dictionary
    '''

    print("Simulating celltag data")

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

    print("Generating CellTag signatures")
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
    print("Checking for duplicates")
    net_dup_pairs = 0
    for i in celltag_sigs.keys():
        df_curr = pd.DataFrame(celltag_sigs[i])
        df_count = df_curr.value_counts()
        df_count = df_count[df_count > 1].copy()

        if(len(df_count) > 0):
            dup_pairs = sum([comb(x,2) for x in df_count.values])
            net_dup_pairs += dup_pairs
            print(f"Duplicate pairs ({i} celltags/cell): {dup_pairs}")
    if(net_dup_pairs > 0):
        print(f"Total duplication rate: {net_dup_pairs/comb(len(cells),2)}\n")
    else:
        print("No duplicates found\n")


    print("Finished!")
    if(return_celltag_sigs):
        return(celltag_sigs)


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


def jaccard_similarities(mat):
    #finds jaccard similarity between all pairs of columns
    cols_sum = mat.getnnz(axis=0)
    ab = mat.T * mat

    # for rows
    aa = np.repeat(cols_sum, ab.getnnz(axis=0))
    # for columns
    bb = cols_sum[ab.indices]

    similarities = ab.copy()
    similarities.data /= (aa + bb - ab.data)

    return(similarities)

def call_clones(jac_mat_low, cells, jac_th=0.6, return_graph = True):
    '''
    Performs CellTag clone calling and returns a clone table
    '''
    jac_mat = jac_mat_low > jac_th

    edge_df = pd.DataFrame(np.array(jac_mat.nonzero()).T)
    vertex_df = pd.DataFrame(cells)
    vertex_df['cb'] = vertex_df[0]
    vertex_df[0] = vertex_df.index
    
    #create clone graph
    g = ig.Graph.DataFrame(edges=edge_df, directed=False, vertices=vertex_df)
    g.vs.select(_degree = 0).delete()
    
    #call clones
    clones = g.components()

    
    #evaluate cols for clone table
    edge_den = [[g.induced_subgraph(i).density()]*len(i) for i in clones]
    clones_bc = [g.induced_subgraph(i).vs['cb'] for i in clones]
    clone_id = [[j+1]*len(i) for j,i in enumerate(clones_bc)]
    
    #create and return clone table
    clone_table = pd.concat([pd.DataFrame((i,j,k), index = ['clone.id','cell.bc','edge.den']).T for i,j,k in zip(clone_id,clones_bc,edge_den)])

    if(return_graph):
        return(g, clone_table)
        
    return(clone_table)

def ident_sparse_clones(clone_info, n_largest = 10, density_th = 0.2, plot=False, **kwargs):
    '''
        Function for sparse clone identification. Takes clone meta data as input and returns sparse clone metadata as output.
        returns clone metadata for sparse clones and optionally plots them on QC plot
    '''
    clone_info_subset = clone_info.nlargest(n_largest, columns='size')
    clone_info_subset = clone_info_subset[clone_info_subset['edge.den'] < density_th].copy()

    if(len(clone_info_subset) == 0):
        print("No sparse clones found!")
        if(plot):
            ax=plot_size_by_den(clone_info, **kwargs)
            return(None, ax)
        return(None)

    if(plot):
        ax=plot_size_by_den(clone_info, red_clones=clone_info_subset['clone.id'].values, **kwargs)
        return(clone_info_subset, ax)
    return(clone_info_subset)

def fix_sparse_clones(clone_graph, sparse_ids = None):
    '''Function to fix sparse clones and rebuild clone table. This takes in the original clone graph and list of sparse clone IDs and returns a new clone table.'''

    if(isinstance(sparse_ids, type(None))):
        print("No sparse clones found")
        return()

    new_clones = []
    edge_density = []
    og_clones = clone_graph.components()

    #define variables for cell number snaity check
    init_cells = sum([len(i) for i in og_clones])
    wasted_cells = 0
    new_cells = 0

    for i in range(len(og_clones)):

        #if sparse clone is encountered
        if(i in sparse_ids-1):
            idx_curr = i

            #extract clone
            sub_gr = clone_graph.induced_subgraph(og_clones[i])

            #extract maximum clique from clone and add to new_clones
            while sub_gr.ecount() > 0:
                new_idx = sub_gr.largest_cliques()[0]
                new_gr = sub_gr.induced_subgraph(new_idx)
                new_cb = new_gr.vs['cb']

                #add new cb to list
                new_clones.append(new_cb)

                #assert edge density is 1 and add to density list
                try:
                    assert(new_gr.density() == 1)
                except AssertionError:
                    print("Maximum clique density is less than 1, graph subsetting is incorrect")
                    return()

                edge_density.append([new_gr.density()]*len(new_cb))

                sub_gr.delete_vertices(new_idx)

            wasted_cells = wasted_cells + sub_gr.vcount()

        #add to clone list as normal if not sparse ID
        else:
            new_clones.append(clone_graph.induced_subgraph(og_clones[i]).vs['cb'])
            edge_density.append([clone_graph.induced_subgraph(og_clones[i]).density()]*len(og_clones[i]))

    #cell number sanity check
    new_cells = sum([len(i) for i in new_clones])

    try:
        assert(new_cells == (init_cells - wasted_cells))
    except AssertionError:
        print("Cell number sanity check failed!")
        return()

    clone_id = [[j+1]*len(i) for j,i in enumerate(new_clones)]
            
    #create and return clone table
    clone_table = pd.concat([pd.DataFrame((i,j,k), index = ['clone.id','cell.bc','edge.den']).T for i,j,k in zip(clone_id,new_clones,edge_density)])
    return(clone_table)




def filter_celltag_table(celltag_mat_met, cells_met, celltags_met):
    '''
    Filters out low confidence 1ctpc signatures from the metric filtered matrix
    '''

    #create 1cptc table
    ct1_cellidx = (celltag_mat_met.sum(axis=1) == 1).nonzero()[0]
    ct1_mtx = celltag_mat_met[ct1_cellidx,]
    ct1_cells = cells_met[ct1_cellidx]

    ct1_tagidx = (ct1_mtx.sum(axis=0) > 0).nonzero()[1]
    ct1_mtx = ct1_mtx[:,ct1_tagidx]
    ct1_celltags = celltags_met[ct1_tagidx]
    assert((ct1_mtx.sum(axis=1) == 1).all())

    #create >1ctpc table
    ctmany_cellidx = (celltag_mat_met.sum(axis=1) > 1).nonzero()[0]
    ctmany_mtx = celltag_mat_met[ctmany_cellidx,]
    ctmany_cells = cells_met[ctmany_cellidx]

    ctmany_tagidx = (ctmany_mtx.sum(axis=0) > 0).nonzero()[1]
    ctmany_mtx = ctmany_mtx[:,(ctmany_mtx.sum(axis=0) > 0).nonzero()[1]]
    ctmany_celltags = celltags_met[ctmany_tagidx]
    assert((ctmany_mtx.sum(axis=1) > 1).all())

    #get all 1cptc tags
    ct1 = set(ct1_celltags)

    #get all >1ctpc tags
    ct_many = set(ctmany_celltags)

    #set of 1cptc tags that aren't >1ctpc
    ct_diff = ct1.difference(ct_many)

    #subset 1ctpc mtx to contain only high confidence celltag sigs
    ct_diff_idx = np.isin(ct1_celltags, list(ct_diff)).nonzero()[0]
    ct1_mtx_fil = ct1_mtx[:,ct_diff_idx]
    ct_diff_idx_2 = (ct1_mtx_fil.sum(axis=1) > 0).nonzero()[0]
    ct1_mtx_fil = ct1_mtx_fil[ct_diff_idx_2,:]

    #create final matrix
    ct1_df = pd.DataFrame.sparse.from_spmatrix(ct1_mtx_fil, index=ct1_cells[ct_diff_idx_2], columns=ct1_celltags[ct_diff_idx]).T
    ctmany_df = pd.DataFrame.sparse.from_spmatrix(ctmany_mtx, index=ctmany_cells, columns=ctmany_celltags).T

    ct_final = ctmany_df.merge(ct1_df, how="outer", right_index=True, left_index=True).fillna(0)
    return(ct_final.T)


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


def get_clone_celltag_mtx(clones, celltag_mat_met, cells_met, celltags_met, sig_type="core"):
    '''
    Takes as input clone table, celltag metric filtered matrix, cells and tags. Returns clone x celltag matrix,cells, tags
    '''
    clone_id = np.empty([0,])
    clone_tag = np.empty([0,])

    for i,j in clones.groupby('clone.id'):

        #get clone_mtx
        clone_mtx_curr = celltag_mat_met[np.isin(cells_met,j['cell.bc']),]

        if(sig_type=="core"):
            ct_sig = celltags_met[(clone_mtx_curr.sum(axis=0) > 1).nonzero()[1]]

        elif(sig_type=="union"):
            ct_sig = celltags_met[(clone_mtx_curr.sum(axis=0) > 0).nonzero()[1]]

        clone_tag = np.hstack([clone_tag, ct_sig])
        clone_id = np.hstack([clone_id, np.ones_like(ct_sig, dtype=int)*i])


    return(table_to_spmtx(clone_id, clone_tag, np.ones_like(clone_id, dtype=int)))



###### PLOTTING FUNCTIONS FOR CLONES ######

def plot_sig_by_size(clone_meta, sig="core", ax=None, palette = 'viridis_r', s=100, edgecolor='black', **kwargs):
    '''
        Function for plotting signature size by clone size colored by edge density
        clone_meta: clone metadata; expected to have a 'size' column and either a "ctpc_core" or "ctpc_union" column
        sig: either "core" or "union"

        returns axes
    '''

    if(ax==None):
        ax=plt.gca()

    if(sig=="core"):
        y_curr = "ctpc_core"
    elif(sig=="union"):
        y_curr = "ctpc_union"
    else:
        print("Sig can only be core or union")
        return()

    #plot scatterplot
    data = clone_meta.sort_values('edge.den', ascending=False)
    plot_curr = plt.scatter(x=data["size"], y=data[y_curr],c=data["edge.den"], cmap=palette, s=s, edgecolor=edgecolor, **kwargs)
    plt.xlabel("size")
    plt.ylabel(y_curr)
    cbar = plt.colorbar(label="edge density")
    # cbar.ax.yaxis.tick_left()
    # cbar.ax.set_ylabel('edge density', rotation=270)
    # plt.tight_layout()

    return(plot_curr)

def plot_size_by_den(clone_meta, red_clones=None, green_clones=None, ax=None, **kwargs):
    '''
        Function for plotting clone size by edge density scatterplot.
        clone_meta: clone metadata: expected to have a 'size', 'edge.den' and 'clone.id' (if red/green clones passed) column
        red_clones: list/ numpy array of clone IDs to be plotted in red (and labelled)
        green_clones: list/ numpy array of clone IDs to be plotted in green (and labelled)
    ''' 
    if(ax==None):
        ax=plt.gca()

    ax = sns.scatterplot(x = "edge.den", y="size", data=clone_meta,color='royalblue', **kwargs)

    #select clone IDs to label on the plot
    if([x for x in (red_clones,green_clones) if x is not None]):
        clone_ids_to_label = np.concatenate([x for x in (red_clones,green_clones) if x is not None])
        clone_meta_sub = clone_meta[clone_meta['clone.id'].isin(clone_ids_to_label)].copy()
        texts = [plt.text(x['edge.den'], x['size'], "clone "+str(int(x['clone.id'])), ha='center', va='center', fontsize=12) for _,x in clone_meta_sub.iterrows()]
        adjustText.adjust_text(texts, expand_points=(1.5, 1.5), arrowprops=dict(arrowstyle="-", color='k', lw=0.5))

    #color cells
    if(red_clones is not None):
        ax=sns.scatterplot(x = "edge.den", y="size", data=clone_meta[clone_meta['clone.id'].isin(red_clones)], color='red', **kwargs)
    if(green_clones is not None):
        ax=sns.scatterplot(x = "edge.den", y="size", data=clone_meta[clone_meta['clone.id'].isin(green_clones)], color='limegreen', **kwargs)

    return(ax)


def plot_clone_upset(clone_table, fname, days=['D3','D12','D21'], title="", **kwargs):
    '''
        This function creates an upset plot from clone table, based on the "day" label. Expects "day" and "clone.id" in clone table. Saves
        output file as a pdf

        clone_table: clone calling table
        fname: file name for pdf
        days: list of days to include in upSet plot
        title: title of the upSet plot
    '''
    idx = [["is_"+x, x] for x in days]
    idx = [item for sublist in idx for item in sublist]
    idx = ['clone.id'] + idx
    clone_upset = pd.DataFrame(columns=idx)
    
    for i,clone in clone_table.groupby("clone.id"):
        list_curr = []
        list_curr.append(i)
        
        for day_str in days:
            day_curr = clone[clone['day'] == day_str]
            
            if(len(day_curr) > 0):
                list_curr.append(True)
            else:
                list_curr.append(False)
            list_curr.append(len(day_curr))
        
        # print(list_curr, idx)
        clone_upset = pd.concat((clone_upset,pd.DataFrame(list_curr, index = idx).T))
        
    clone_upset = clone_upset.set_index(["is_"+x for x in days])
    
    # return(clone_upset)
    upset = upsetplot.UpSet(clone_upset, subset_size='count', intersection_plot_elements=2, element_size=57, sort_by="cardinality",
                  show_counts = True, totals_plot_elements=1, orientation="horizontal", **kwargs)

    for day_str in days:
        upset.add_catplot(value=day_str, kind='strip', size=3, elements=3)
    
    upset.plot()
    plt.title(title)
    plt.savefig(fname)

