import numpy as np
import pandas as pd
import os
import h5py
import scipy.stats as stats
import sys

#INPUT_FILENAMES = ['c2.all.v7.0.symbols_JustK.gmt','c2.all.v7.0.symbols_JustR.gmt', 'x', 'c5.bp.v7.0.symbols_SHORT.gmt', 'c5.mf.v7.0.symbols.gmt']
# python gmts_to_adj_matrices_offdiagonal.py c2.all.v7.0.symbols_JustK.gmt x c2.all.v7.0.symbols_JustR.gmt


if (len(sys.argv)<2):
    print("USAGE: gmts_to_adj_matrices_offdiagonal.py [list all gmt files for rows] x [list all gmt files for columns]")

print(sys.argv)    
INPUT_FILENAMES = sys.argv[1:]

if ("x" not in np.array([x.lower() for x in INPUT_FILENAMES])):
    print("USAGE: gmts_to_adj_matrices_offdiagonal.py [list all gmt files for rows] x [list all gmt files for columns]") 


splitloc = np.where(np.array([x.lower() for x in INPUT_FILENAMES]) == "x")[0][0]

row_filenames = INPUT_FILENAMES[:splitloc]
col_filenames = INPUT_FILENAMES[splitloc+1:]




############ LOAD IN GMT FILES -->  "names" = pathway names, "genes" = list of genes in pathway 

path_to_pways = "../pathways_raw/"


gmt_dfs_rows = []
for f in row_filenames:
    gmt = pd.read_csv(path_to_pways + f, header=None)
    gmt["names"] = gmt[0].apply(lambda x: x.split("\t")[0])
    gmt["genes"] = gmt[0].apply(lambda x: x.split("\t")[2:])
    
    gmt_dfs_rows.append(gmt)
    
    
pathways_df_rows = pd.concat(gmt_dfs_rows).drop_duplicates(0)

gmt_dfs_columns = []
for f in col_filenames:
    gmt = pd.read_csv(path_to_pways + f, header=None)
    gmt["names"] = gmt[0].apply(lambda x: x.split("\t")[0])
    gmt["genes"] = gmt[0].apply(lambda x: x.split("\t")[2:])
    
    gmt_dfs_columns.append(gmt)
    
    
pathways_df_columns = pd.concat(gmt_dfs_columns).drop_duplicates(0)

all_genes = np.union1d(np.unique(np.hstack(pathways_df_rows["genes"].values)), np.unique(np.hstack(pathways_df_columns["genes"].values)))
print("# genes: %i"%len(all_genes))
print("# pathways: %i rows, %i columns"%(len(pathways_df_rows), len(pathways_df_columns)))



#########  

#########  

odds_ratios = np.zeros([len(pathways_df_rows), len(pathways_df_columns)])
pvals = np.zeros([len(pathways_df_rows), len(pathways_df_columns)])

for i,pway1 in pathways_df_rows.iterrows():
    if i%10==0:
        print("%i of %i"%(i, len(pathways_df_rows)))
    
    for j,pway2 in pathways_df_columns.iterrows():
        if pway1["names"] == pway2["names"]:
            odds_ratios[i,j] = odds_ratios[j,i]
            pvals[i,j] = pvals[j,i]

        else:
            mx = np.zeros([2,2])
            mx[0,0] = len(np.intersect1d(pway1["genes"],pway2["genes"]))
            mx[1,0] = len(np.setdiff1d(pway1["genes"],pway2["genes"]))
            mx[0,1] = len(np.setdiff1d(pway2["genes"],pway1["genes"]))
            mx[1,1] = len(all_genes)-len(np.union1d(pway1["genes"],pway2["genes"]))
            

            odds_ratios[i,j], pvals[i,j] = stats.fisher_exact(mx)


            
pway_subfolder = "-".join([x[:-4] for x in INPUT_FILENAMES])
if not os.path.isdir("adj_matrices/offdiagonal/" +pway_subfolder ):
    os.makedirs("adj_matrices/offdiagonal/" +pway_subfolder )
np.savetxt("adj_matrices/offdiagonal/" +pway_subfolder + "/odd_ratios.txt", odds_ratios)
np.savetxt("adj_matrices/offdiagonal/" + pway_subfolder + "/pvals.txt", pvals)
np.savetxt("adj_matrices/offdiagonal/" +pway_subfolder + "/pathway_names_rows.txt", pathways_df_rows["names"].values, fmt="%s")
np.savetxt("adj_matrices/offdiagonal/" +pway_subfolder + "/pathway_names_columns.txt", pathways_df_columns["names"].values, fmt="%s")

print("SAVED TO FOLDER:  %s" %("adj_matrices/offdiagonal/" +pway_subfolder))