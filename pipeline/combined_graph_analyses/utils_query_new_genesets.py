import pandas as pd
import numpy as np
import os
import networkx as nx
from scipy import stats
import copy

def reverse_dict(dic, assume_unique=False):
    newdict = {}
    for val in np.unique(list(dic.values())):
        newdict[val] = np.array([k for k,v in dic.items() if v == val])
        if assume_unique:
            newdict[val] = newdict[val][0]
    return newdict

def load_pathway_df_genes(gsets_folders, path_to_pways):
    gmt_dfs = []
    for f in gsets_folders:
        gmt = pd.read_csv(path_to_pways + f + ".gmt", header=None)
        gmt["names"] = gmt[0].apply(lambda x: x.split("\t")[0])
        gmt["genes"] = gmt[0].apply(lambda x: x.split("\t")[2:])
        gmt_dfs.append(gmt)

    pathways_df = pd.concat(gmt_dfs).drop_duplicates(0)
    pathways_df = pathways_df.reset_index(drop=True)
    all_genes = np.unique(np.hstack(pathways_df["genes"].values))
    return pathways_df, all_genes


def load_pathways_G_with_weights(path_to_adj_matrices, gsets_folders, pathways_df):
    pvals_allrows = []
    for i, foldername1 in enumerate(gsets_folders):
        print(foldername1)
        pvals_row_blocks = []

        for j, foldername2 in enumerate(gsets_folders):
            if foldername1 == foldername2:
#                 print("\t",foldername1, len(np.loadtxt(os.path.join(path_to_adj_matrices, foldername1, "pathway_names.txt"), dtype=str)))
                pvals_row_blocks.append(np.loadtxt(os.path.join(path_to_adj_matrices, foldername1, "pvals.txt")))

            # rows and cols match up with file we read in
            elif "%s--%s"%(foldername1, foldername2) in os.listdir(os.path.join(path_to_adj_matrices, "offdiagonal")):
#                 print("\t", foldername1, foldername2)
                pvals_row_blocks.append(np.loadtxt(os.path.join(path_to_adj_matrices, "offdiagonal", "%s--%s"%(foldername1, foldername2), "pvals.txt")))

            # need to transpose
            elif "%s--%s"%(foldername2, foldername1) in os.listdir(os.path.join(path_to_adj_matrices, "offdiagonal")):
#                 print("\ttransposed!", foldername1, foldername2)
                pvals_row_blocks.append(np.loadtxt(os.path.join(path_to_adj_matrices, "offdiagonal", "%s--%s"%(foldername2, foldername1), "pvals.txt")).T)

#             else:
# #                 print("MISSING MATRIX: %s vs %s"%(foldername1, foldername2))
        pvals_allrows.append(np.hstack(pvals_row_blocks))
    
    pvals = np.vstack(pvals_allrows)
    pathway_names = pathways_df["names"].values
    newzero = np.min(pvals[np.nonzero(pvals)]) / 10
    pvals[pvals == 0] = newzero
    weights =  -1 * np.log10(pvals)
    weights[np.isnan(weights)] = 0.0
    weights[np.isinf(weights)] = 0.0
    weights[weights == 0] = 0.0
    weights_df = pd.DataFrame(weights.astype(float), index=pathway_names, columns=pathway_names)
    G_with_weights = nx.from_numpy_matrix(weights_df.values)
    return pathway_names, G_with_weights

# used formula from here: https://en.wikipedia.org/wiki/Louvain_modularity, and checked that it agrees with 
# results from community_louvain.modularity(candidate_partition, query_G) - original_modularity 

def d_Q(query_G, coms, node):
    k_i = query_G.degree(node, weight="weight")
    m = query_G.size("weight")
    
    d_Qs = []
    for com_idx in np.unique(coms):
        com_pways = np.where(coms==com_idx)[0]
        sigma_in = query_G.subgraph(com_pways).size("weight")
        sigma_tot = np.sum([v for (n, v) in query_G.degree(com_pways, weight="weight")])
        i_in_subG = query_G.subgraph(np.append(com_pways, node))
        k_i_in = i_in_subG.degree(node, weight="weight")
        dq = ((sigma_in + 2*k_i_in)/(2*m) - ((sigma_tot+k_i)/(2*m))**2)  - ((sigma_in/(2*m)- (sigma_tot/(2*m))**2 -  (k_i/(2*m))**2))
        d_Qs.append(dq)
    return np.array(d_Qs)

def query_modularities(pathways_df, G_with_weights, all_genes, coms, query_pway_genes):
    query_odds_ratios = np.zeros([len(pathways_df)])
    query_pvals = np.zeros([len(pathways_df)])

    for i,pway in pathways_df.iterrows():

        mx = np.zeros([2,2])
        mx[0,0] = len(np.intersect1d(query_pway_genes,pway["genes"]))
        mx[1,0] = len(np.setdiff1d(query_pway_genes,pway["genes"]))
        mx[0,1] = len(np.setdiff1d(pway["genes"],query_pway_genes))
        mx[1,1] = len(all_genes)-len(np.union1d(query_pway_genes,pway["genes"]))

        query_odds_ratios[i], query_pvals[i] = stats.fisher_exact(mx)


    newzero = np.min(query_pvals[np.nonzero(query_pvals)]) / 10
    query_pvals[query_pvals == 0] = newzero
    query_weights =  -1 * np.log10(query_pvals)
    query_weights[np.isnan(query_weights)] = 0.0
    query_weights[np.isinf(query_weights)] = 0.0
    query_weights[query_weights == 0] = 0.0

    query_G = copy.deepcopy(G_with_weights)
    query_node_idx = len(pathways_df)
    query_G.add_node(query_node_idx)
    query_G.add_weighted_edges_from([(query_node_idx,i,x)  for i,x in enumerate(query_weights) if x!=0])

    tmp_coms = copy.deepcopy(coms)
    tmp_coms = np.append(tmp_coms, np.max(np.unique(coms))+1)
    dq_array = d_Q(query_G, tmp_coms, len(pathways_df))
    
    # normalize result so that we're showing the change compared to being in own community
    candidate_dQs = dq_array[:-1] - dq_array[-1]

    return query_weights, candidate_dQs