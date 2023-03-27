#define helper functions
import pandas as pd
import numpy as np
import igraph as ig
import scipy

import adjustText
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec


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

#def plot_clone_upset(clone_table, fname, days=['D3','D12','D21'], title="", **kwargs):
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

