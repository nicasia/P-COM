import networkx as nx

# CNM:
import CNM_networkx as community_CNM

# Louvain:
import community as community_louvain

# agglomerative & spectral clustering:
from sklearn.cluster import SpectralClustering, AgglomerativeClustering
import numpy as np


######## Community detection / clustering algorithms ############

def get_labels_CNM(G, **kwargs):
    if "num" in kwargs.keys():
        if kwargs["num"]:
            frozensets = community_CNM.greedy_modularity_communities_num(G,kwargs["num"])
        else:
            frozensets = community_CNM.greedy_modularity_communities(G)
    else:
        frozensets = community_CNM.greedy_modularity_communities(G)

    community_labels = com_labels_from_frozensets(frozensets,G)
    return community_labels


def get_labels_louvain(G, **kwargs):   
    if "resolution" in kwargs.keys():
        if "weight" in kwargs.keys():
            partition = community_louvain.best_partition(G, weight=kwargs["weight"], resolution=kwargs["resolution"])
        else:
            partition = community_louvain.best_partition(G, resolution=kwargs["resolution"])
    else:
        partition = community_louvain.best_partition(G)
    return np.array(list(partition.values()))

def get_labels_agglomerative(pvals, num, **kwargs):
    
    kw = {"l": "precomputed", "linkage": "complete"}
    for key in kwargs.keys():
        kw[key] = kwargs[key]
    
    community_labels = AgglomerativeClustering(affinity=kw['affinity'], n_clusters=num, linkage=kw["linkage"]).fit_predict(pvals)
    return community_labels

def get_labels_spectral(pvals, num, **kwargs):
    # DEFINE SIMILARITY AS 1-PVAL  (pval of 0=>sim=1,   pval of 1=> sim=0)
    sim_opp_pvals = 1-pvals

    # if we want fully connected graph (for spectral clustering), set 0 similarity to 1/10th of smallest non-zero
    if 0 in sim_opp_pvals:
        sim_opp_pvals[np.where(sim_opp_pvals==0)]=np.unique(sim_opp_pvals)[1]/10
        
    community_labels = SpectralClustering(num, affinity='precomputed').fit_predict(sim_opp_pvals)
    return community_labels


######### Mapping patitions / lists of nodes to community labels ######

def com_labels_from_frozensets(frozensets, G):
    num_nodes = len(list(G.nodes))
    com_labels = np.empty(num_nodes)
    for lab, nodelist in  enumerate(frozensets):
        for elt in list(nodelist):
            com_labels[elt] = lab
    return com_labels

