import requests
from bs4 import BeautifulSoup
import numpy as np
import pandas as pd
import os

###### Dealing with dictionaries 

def reverse_dict(dic, assume_unique=False):
    newdict = {}
    for val in np.unique(list(dic.values())):
        newdict[val] = np.array([k for k,v in dic.items() if v == val])
        if assume_unique:
            newdict[val] = newdict[val][0]
    return newdict

# reverse dict but where values in the dict are lists
def reverse_dict_lists(dic, assume_unique=False):
    newdict = {}
    for val in np.unique(np.hstack(list(dic.values()))):
        newdict[val] = np.array([k for k,v in dic.items() if val in v])
        if assume_unique:
            newdict[val] = newdict[val][0]
    return newdict


#### Loading the graph

def load_graph_edges(gsets_folders, adj_mx_path, pathways):
    pvals_allrows = []
    odds_ratios_all_rows = []
    for i, foldername1 in enumerate(gsets_folders):
        print(foldername1)

        pvals_row_blocks = []
        odds_ratios_row_blocks = []

        for j, foldername2 in enumerate(gsets_folders):

            if foldername1 == foldername2:
                print("\t",foldername1, len(np.loadtxt("%s/%s/pathway_names.txt"%(adj_mx_path,foldername1), dtype=str)))
                pvals_row_blocks.append(np.loadtxt("%s/%s/pvals.txt"%(adj_mx_path,foldername1)))
                odds_ratios_row_blocks.append(np.loadtxt("%s/%s/odd_ratios.txt"%(adj_mx_path,foldername1)))

            # rows and cols match up with file we read in
            elif "%s--%s"%(foldername1, foldername2) in os.listdir("%s/offdiagonal/"%(adj_mx_path)):
                print("\t", foldername1, foldername2)
                pvals_row_blocks.append(np.loadtxt("%s/offdiagonal/%s--%s/pvals.txt"%(adj_mx_path, foldername1, foldername2)))
                odds_ratios_row_blocks.append(np.loadtxt("%s/offdiagonal/%s--%s/odd_ratios.txt"%(adj_mx_path, foldername1, foldername2)))

            # need to transpose
            elif "%s--%s"%(foldername2, foldername1) in os.listdir("%s/offdiagonal/"%adj_mx_path):
                print("\ttransposed!", foldername1, foldername2)
                pvals_row_blocks.append(np.loadtxt("%s/offdiagonal/%s--%s/pvals.txt"%(adj_mx_path, foldername2, foldername1)).T)
                odds_ratios_row_blocks.append(np.loadtxt("%s/offdiagonal/%s--%s/odd_ratios.txt"%(adj_mx_path, foldername2, foldername1)).T)

            else:
                print("MISSING MATRIX: %s vs %s"%(foldername1, foldername2))
        pvals_allrows.append(np.hstack(pvals_row_blocks))
        odds_ratios_all_rows.append(np.hstack(odds_ratios_row_blocks))


    pvals = np.vstack(pvals_allrows)
    odds_ratios = np.vstack(odds_ratios_all_rows)
    pathway_names = np.hstack([pathways[name] for name in gsets_acronyms])


    pvals[pvals == 0] = np.min(pvals[np.nonzero(pvals)]) / 10
    weights =  -1 * np.log10(pvals)
    weights[np.isnan(weights)] = 0.0
    weights[np.isinf(weights)] = 0.0
    weights[weights == 0] = 0.0
    
    return weights


#### Getting MSigDB descriptions from the web

def msigdb_to_description(pway):
    url = "http://software.broadinstitute.org/gsea/msigdb/cards/%s"%pway
    response = requests.get(url)
    soup = BeautifulSoup(response.text, "html.parser")
    td_tags = soup.find_all("td")
    
    if ("Gene Set Not Found" in td_tags[2].text):
        print("NOT FOUND")
        return pway
    else:
        return td_tags[5].text
    
###### Loading gene sets & categories

acronym_to_folder = {"KEGG": "c2.all.v7.0.symbols_JustK", "REACTOME":"c2.all.v7.0.symbols_JustR",
                  "GO_BP": "c5.bp.v7.0.symbols_SHORT", "GO_CC": "c5.cc.v7.0.symbols", "GO_MF": "c5.mf.v7.0.symbols"}
folder_to_acronym = reverse_dict(acronym_to_folder, assume_unique = True)

pway_subfolders =  'c2.all.v7.0.symbols_JustK-c2.all.v7.0.symbols_JustR-c5.bp.v7.0.symbols_SHORT-c5.mf.v7.0.symbols'

gsets_folders = pway_subfolders.split("-")
gsets_acronyms = [folder_to_acronym[x] for x in gsets_folders]



def load_curated_labels(pway_dbs=["KEGG", "REACTOME", "GO_BP", "GO_MF"]):
    true_labels = {}
    true_labels_names = {}
    true_labels_unique_names = []
    true_labels_unique = []
    max_so_far = 0 

    gmts = {}
    for pathway_group in pway_dbs:
        if pathway_group=="KEGG":
            pway_subfolder = "c2.all.v7.0.symbols_JustK"
            level = 0
        elif pathway_group=="REACTOME":
            pway_subfolder = "c2.all.v7.0.symbols_JustR"
            level = 0
        elif pathway_group=="GO_CC":
            pway_subfolder = "c5.cc.v7.0.symbols"
            level = 1
        elif pathway_group=="GO_MF":
            pway_subfolder = "c5.mf.v7.0.symbols"
            level = 1
        elif pathway_group=="GO_BP":
            pway_subfolder = "c5.bp.v7.0.symbols_SHORT" 
            level = 1

        ##############
        #Read true labels
        if level != 0:
            labels_df = pd.read_csv("../curated_labels/%s_labels_LEVEL%d.csv"%(pathway_group, level))
            toplev_names = np.loadtxt("../curated_labels/%s_categorynames_LEVEL%d.csv"%(pathway_group, level), dtype="str", delimiter="\n")
        else:
            labels_df = pd.read_csv("../curated_labels/%s_labels.csv"%(pathway_group))
            toplev_names = np.loadtxt("../curated_labels/%s_categorynames.csv"%(pathway_group), dtype="str", delimiter="\n")

        true_labels[pathway_group] = labels_df["label"].values
        true_labels_names[pathway_group] = labels_df["category"].values


        num_true_cats = len(np.unique(true_labels[pathway_group]))
        print("Number of true categories: ", num_true_cats)

        print(np.unique(true_labels[pathway_group]))

        adjusted_labels = max_so_far + labels_df["label"].values + 1
        true_labels_unique.append(adjusted_labels)
        max_so_far = np.max(adjusted_labels)
        true_labels_unique_names.append(true_labels_names)


        gmt = pd.read_csv("../../pathways_raw/%s.gmt"%pway_subfolder, header=None)
        gmt["names"] = gmt[0].apply(lambda x: x.split("\t")[0])
        gmt["urls"] = gmt[0].apply(lambda x: x.split("\t")[1])
        gmt["genes"] = gmt[0].apply(lambda x: x.split("\t")[2:])

        gmts[pathway_group] = gmt

    true_labels_unique = np.hstack(true_labels_unique)
    return gmts, true_labels_unique, true_labels, true_labels_names