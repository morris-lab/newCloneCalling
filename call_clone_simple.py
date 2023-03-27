import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
import clone_calling_fn as cc


from scipy import io

#PARAMS
TRIPLET_TH = 1
STARCODE_TH = 2
BIN_TH = 1
METRIC_LOW = 1
METRIC_HIGH = 25

#Create list of files to import for celltag matrix
KEYS = ['sample1_RNA','sample2_RNA']

ct_reads_list = []

#import and filter celltag table from each file and add to list
for KEY_CURR in KEYS:
    
    print(f"processing: {KEY_CURR}")
    bam_fil = pd.read_csv(f"../celltag_reads/{KEY_CURR}/{KEY_CURR}_bam_parse.txt", sep="\t")
    print(f"Total filtered CellTag Reads: {len(bam_fil)}")
    

    #create UMI counts per CB-celltag pair
    bam_fil['concat'] = bam_fil['Cell.BC'] + "." + bam_fil['Cell.Tag']+ "." + bam_fil['UMI']
    
    #filter triplets on read counts
    reads, counts = np.unique(bam_fil['concat'].values, return_counts=True)
#     plt.hist(counts,bins=100)
#     plt.tight_layout()
    bam_umi = pd.DataFrame(reads[(counts > TRIPLET_TH)])
    
    seq_sat = 100*(1 - len(reads[(counts == 1)])/len(bam_fil))
    print("CellTag Sequencing saturation: ", str(seq_sat))
    
    bam_umi['Cell.BC'] = bam_umi[0].apply(lambda x:str(x).split(".")[0])
    bam_umi['Cell.Tag'] = bam_umi[0].apply(lambda x:str(x).split(".")[1])
    bam_umi.drop(columns=0,inplace=True)
    del bam_fil
    
    #starcode collapse
    (bam_umi['Cell.BC'].apply(lambda x: x[:-2]) + bam_umi['Cell.Tag']).to_csv("collapsing.txt",
                                                                             sep='\t',
                                                                             index=False,
                                                                             header=False)

    os.system("~/starcode/starcode -t 4 -d {} -s collapsing.txt > collapsing_result.txt".format(str(STARCODE_TH)))
    ct_reads_final = pd.read_csv("collapsing_result.txt", sep='\t', header=None)
    ct_reads_final['CB'] = ct_reads_final[0].apply(lambda x: x[:16] + "-1")
    ct_reads_final['celltag'] = ct_reads_final[0].apply(lambda x: x[16:])
    ct_reads_final.rename(columns={1:"count"}, inplace = True)
    ct_reads_final.drop(columns=[0], inplace = True)
    os.system('rm collapsing_result.txt')
    os.system('rm collapsing.txt')
    
    #allowlisting
    allow = pd.read_csv("~/misc_files/18N-multi-v1-allowlist.csv", header = None)
    ct_reads_final = ct_reads_final[ct_reads_final['celltag'].isin(allow[0])].copy()
    ct_reads_final['CB'] = KEY_CURR + "-" + ct_reads_final['CB']
    ct_reads_list.append(ct_reads_final)
    print()



KEY_CURR = "RNA-merged"
if(not os.path.isdir("../proc_files/{0}".format(KEY_CURR))):
    os.makedirs("../proc_files/{0}".format(KEY_CURR))



#merge all celltag tables into 1 table and save to disk
ct_reads_merged = pd.concat(ct_reads_list)
ct_reads_merged.to_csv("../proc_files/{0}/{0}_all_celltag_reads.csv".format(KEY_CURR))

#create allowlisted celltag UMI count matrix
celltag_mat, cells, celltags = cc.table_to_spmtx(ct_reads_merged['CB'],
                                                       ct_reads_merged['celltag'],
                                                       count_data=ct_reads_merged['count'])

#write allowlisted celltag UMI count matrix to file
print("Writing allowlisted matrix to file")
io.mmwrite("../proc_files/{0}/{0}_allow_ctmat.mtx".format(KEY_CURR), celltag_mat)
np.savetxt("../proc_files/{0}/{0}_allow_cells.txt".format(KEY_CURR),cells, delimiter='\t', fmt='%s')
np.savetxt("../proc_files/{0}/{0}_allow_celltags.txt".format(KEY_CURR),celltags, delimiter='\t', fmt='%s')


#binarize
celltag_mat_bin = celltag_mat > BIN_TH

#metric filter
row_fil = ((celltag_mat_bin.sum(axis=1) > METRIC_LOW) & (celltag_mat_bin.sum(axis=1) < METRIC_HIGH))
temp = celltag_mat_bin[row_fil.nonzero()[0],]
col_fil = temp.sum(axis=0) > 0
celltag_mat_met = temp[:,col_fil.nonzero()[1]].copy()
celltag_mat_met = celltag_mat_met*1

cells_met = np.array(cells)[row_fil.nonzero()[0]]
celltags_met = np.array(celltags)[col_fil.nonzero()[1]]

#write metric filtered matrix to file
print("Total cells: " + str(celltag_mat_met.shape[0]))
print("Total celltags: " + str(celltag_mat_met.shape[1]))
io.mmwrite("../proc_files/{0}/{0}_met_ctmat.mtx".format(KEY_CURR), celltag_mat_met)
np.savetxt("../proc_files/{0}/{0}_met_cells.txt".format(KEY_CURR),cells_met, delimiter='\t', fmt='%s')
np.savetxt("../proc_files/{0}/{0}_met_celltags.txt".format(KEY_CURR),celltags_met, delimiter='\t', fmt='%s')


#jacard similarity
print("Calculating Jaccard similarity")
jac_mat = cc.jaccard_similarities(celltag_mat_met.tocsc().astype(np.float64).transpose())
jac_mat.setdiag(0)              


g, clones = cc.call_clones(scipy.sparse.tril(jac_mat), cells_met)

r2 = celltag_mat_met
fig, ax = plt.subplots(2,2, figsize=(10,10))
ax[0,0].scatter(np.expand_dims(np.arange(0,r2.shape[0]),axis=1), np.array(r2.sum(axis=1)), rasterized=True)
ax[0,1].scatter(np.expand_dims(np.arange(0,r2.shape[1]),axis=1), np.array(r2.sum(axis=0)), rasterized=True)
ax[1,0].hist(np.array(r2.sum(axis=1)), bins=100, rasterized=True)
# ax[1,0].set_ylim((0,50))
ax[1,1].hist(np.array(r2.sum(axis=0).transpose()), bins=100, rasterized=True)

ax[0,0].set_title("celltags/cell")
ax[0,1].set_title("cells/celltag")

fig.suptitle("metric filtered mtx QC", fontsize=19)
fig.tight_layout()


#write files to disk
plt.savefig("./metric_filter_QC_plots.pdf")
clones.write.csv("./clones.csv")





