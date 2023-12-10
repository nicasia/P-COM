import numpy as np
import pandas as pd
import os
import h5py
import scipy.stats as stats
import sys

#INPUT_FILENAMES = ['c2.all.v7.0.symbols_JustK.gmt']
# gmts_to_adj_matrices.py c2.all.v7.0.symbols_JustK.gmt

if len(sys.argv) < 2:
    print("USAGE:  gmts_to_adj_matrices.py [list all gmt files with spaces in between]")
    

print(sys.argv)    
INPUT_FILENAMES = sys.argv[1:]



############ LOAD IN GMT FILES -->  "names" = pathway names, "genes" = list of genes in pathway 

path_to_pways = "../pathways_raw/"


gmt_dfs = []
for f in INPUT_FILENAMES:
    gmt = pd.read_csv(path_to_pways + f, header=None)
    gmt["names"] = gmt[0].apply(lambda x: x.split("\t")[0])
    gmt["genes"] = gmt[0].apply(lambda x: x.split("\t")[2:])
    
    gmt_dfs.append(gmt)
    
    
pathways_df = pd.concat(gmt_dfs).drop_duplicates(0)
all_genes = np.unique(np.hstack(pathways_df["genes"].values))
print("# pathways: %i"%len(pathways_df))
print("# genes: %i"%len(all_genes))


#########  

odds_ratios = np.zeros([len(pathways_df), len(pathways_df)])
pvals = np.zeros([len(pathways_df), len(pathways_df)])

for i,pway1 in pathways_df.iterrows():
    
    if i%10==0:
        print("%i of %i"%(i, len(pathways_df)))
    
    for j,pway2 in pathways_df.iterrows():
        if i > j:
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
if not os.path.isdir("adj_matrices/" +pway_subfolder ):
    os.makedirs("adj_matrices/" +pway_subfolder )
np.savetxt("adj_matrices/" +pway_subfolder + "/odd_ratios.txt", odds_ratios)
np.savetxt("adj_matrices/" + pway_subfolder + "/pvals.txt", pvals)
np.savetxt("adj_matrices/" +pway_subfolder + "/pathway_names.txt", pathways_df["names"].values, fmt="%s")

print("SAVED TO FOLDER:  %s" %("adj_matrices/" +pway_subfolder))