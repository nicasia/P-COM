# An automatic integrative method for learninginterpretable communities of biological pathways
### Nicasia Beebe-Wang, Ayse B. Dincer, and Su-In Lee

##### Paul G. Allen School of Computer Science & Engineering, University of Washington, Seattle

Although knowledge of biological pathways is essential for interpreting results from computational biology studies, the growing number of pathway databases complicates efforts to efficiently perform pathway analysis due to high redundancies among pathways from different databases, and inconsistencies in how pathways are created and named. We introduce the PAthway Communities (PAC) framework, which reconciles pathways from different databases and reduces pathway redundancy by revealing informative groups with distinct biological functions. Uniquely applying the Louvain community detection algorithm to a network of 4,847 pathways from KEGG, REACTOME and Gene Ontology databases, we identify 35 distinct and automatically annotated communities of pathways and show that they are consistent with expert-curated pathway categories. Further, we demonstrate that our pathway community network can be queried with new gene sets to provide biological context in terms of related pathways and communities. Our approach, combined with an interpretable web tool we provide, will help computational biologists more efficiently contextualize and interpret their biological findings. 

<img align="center" src="Concept_Figure.jpg" width="50%">.

---

#### Resources

- **Community_Members.xlsx** and **Community_Members.csv** files contain the list of all pathways included in each of the 35 pathway communities we learned.  
- **Community_kmers.xlsx** and **Community_kmers.csv** files contain the list of k-mers for each pathway community, along with the number of occurence and hubness of each pathway (within the community's subgraph).   

---

## Installation Guide 

For a full list of software packages and version numbers, see the Conda environment file environment.yml. We recommend installation of the required packages using the Conda package manager, available through the Anaconda Python distribution. Anaconda is available free of charge for non-commercial use through Anaconda Inc. After installing Anaconda and cloning this repository, use the conda command to install necessary packages: conda env create -f environment.yml

#### Datasets

[Molecular Signatures Database (MSigDB)](http://www.gsea-msigdb.org/gsea/downloads.jsp#msigdb) is a collection of annotated gene sets, which include KEGG, REACTOME, GO, and Hallmark gene sets. Curated categories for KEGG pathways are available [here](https://www.genome.jp/kegg/pathway.html). Curated categories for REACTOME pathways are available from their [interactive web browser](https://reactome.org/PathwayBrowser/). GO relations are provided [here](http://geneontology.org/docs/download-ontology/) (we used the go-basic version).

Data from the METABRIC (Molecular Taxonomy of Breast Cancer International Consortium) cohort are available from the cBioPortal for cancer genomics. The specific dataset used in our study was downloaded [here](https://cbioportal-datahub.s3.amazonaws.com/brca_metabric.tar.gz), and an interactive view of clinical features is available [here](https://www.cbioportal.org/study/summary?id=brca_metabric).


#### Setting up the datasets
Our data is organized as follows. If you would like to re-run the pipeline with new gene sets, or perform a post-hoc analysis with a different dataset, you will need to update the files in the folders below accordingly: 

    ./Datasets/
        data_clinical_patient.txt
        data_clinical_sample.txt
        data_mRNA_median_Zscores.txt
    ./pathways_raw/
        [all .gmt files from pathway databases go in this folder]
    ./pipeline/curated_labels/
        [.csv files with a “pathway” and “category” column which maps each pathway to a category] 

A note on GMT files:  GMT is a format provided by MSigDB downloads, and is commonly used for pathway enrichment analysis. If you have your gene sets saved in a different format, you can create .gmt file with the following basic rules [here](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29). 

## Code 

### PAC Overview
Below, we show how to run the code to reproduce results from our paper. The two main components are the following:
 
**Generating PAC:**
1. Creating the pathway graphs from gene sets (generating adjacency matrices by calculating pairwise overlaps between pathways)
2. Learning communities of pathways (our final approach was the Louvain algorithm, but others are possible as well -- e.g., Figure 2 in our paper)

**Exploring what PAC has learned:**
- Evaluating correspondence between our communities and Hallmark pathways (Figure 3) or curated categories from the original databases (Supplementary Figures 9-13)
- Visualizing the learned community network and subgraphs for each community (Figure 4).
- Querying a new gene set against PAC’s learned communities (Figure 5)

The "pipeline" folder contains all the scripts in our pipeline used to generate PAC. The "pipeline/combined_graph_analyses" folder contains all notebooks used to analyze learned communities (e.g., Figures 3-5). Below, we describe all scripts used (for .ipynb files, a brief overview is provided at top of each jupyter notebook).

### Generating PAC

#### Step 0: Understanding pathway redundancy
**pathway_overlaps.ipynb** is an analysis notebook which measures the overlaps among pathways for the visualization in Figure 1. This provides an overview of how much redundancy exists in the gene sets, and is not part of the pipeline itself. 

#### Step 1: Generating adjacency matrices
- **gmts_to_adj_matrices.py** and **gmts_to_adj_matrices_offdiagonal.py** define adjacency matrices for the MSigDB pathways by measuring the pairwise similarities between pathways from each database and across databases, respectively. In our paper, the specific geneset files used were: KEGG (c2.all.v7.0.symbols_JustK.gmt), REACTOME (c2.all.v7.0.symbols_JustR.gmt),  GO Molecular Function (c5.mf.v7.0.symbols.gmt) and GO Biological Processes (c5.bp.v7.0.symbols_SHORT)
- **Gmts_to_adj_matrices.py** usage: 
    - Input: list of gmt files (contained in the  “pathways_raw” folder) separated by spaces. For our paper, we ran this once for each .gmt file separately since our initial evaluations were for each individual database (i.e., step 2 below) (For example: `python gmts_to_adj_matrices.py c2.all.v7.0.symbols_JustK.gmt`)
    - Output: Fisher’s exact test p-values are saved to `./pipeline/adj_matrices/[geneset names]/pvals.txt`
- **gmts_to_adj_matrices_offdiagonal.py** usage: 
    - Input: two lists of gmt files separated by “x” (contained in the  “pathways_raw” folder).  The script calculates Fisher’s exact tests between each pathway before the “x” versus pathways after the “x”. For our paper, we ran this once for pairs of .gmt file separately (e.g., `python gmts_to_adj_matrices_offdiagonal.py c2.all.v7.0.symbols_JustK.gmt x c2.all.v7.0.symbols_JustR.gmt` to compute a KEGG vs. REACTOME adjacency matrices)
    - Output: Fisher’s exact test p-values are saved to `./pipeline/adj_matrices/offdiagonal/[geneset names]--[genset names]/pvals.txt`
- If you choose to re-run these scripts, the newly produced adjacency matrices will overwrite those currently saved in this repository. If you'd later like to ensure that the results here are consistent with those presented in our paper, you may evaluate them in the `combined_graph_louvain_with_weights.ipynb` notebook. 

#### Learning communities 

##### Comparison of clustering algorithms (Figure 2 in our paper):
While our final method uses the Louvain algorithm to learn communites of related pathways, we evaluated several other alternative approaches, separately for each database, with the following code:
- **algorithm_helpers.py**  is the helper script to run various clustering and community detection algorithms. Each of the functions in this script serves as a wrapper around the algorithms such that the resulting clusters are in the same format for easy comparison. (These functions are called in the `comparison_of_algorithms_NMI.ipynb` notebook below.)
- **CNM_networkx.py** is a modified version of the CNM algorithm (which allows us to select the number of communities to generate) originally from [NetworkX]( https://networkx.github.io/documentation/stable/_modules/networkx/algorithms/community/modularity_max.html).  The original version of the algorithm is under the function `greedy_modularity_communities` - the version we use for the paper is `greedy_modularity_communities_num` which we have slightly modified to control the number of communities learned. This function is called in the `comparison_of_algorithms_NMI.ipynb` notebook below, and takes as input a NetworkX graph and returns community labels. 
- **comparison_of_algorithms_NMI.ipynb** script executes all the clustering algorithms and compares them across all pathway databases (Figure 2).  

##### Generating combined pathway network and learning communities
**combined_graph_louvain_with_weights.ipynb** defines the combined pathway network (the user may specify a list of databases to include) and applies the Louvain algorithm to learn pathway communities. In this notebook, the user may also run tests to compare their adjacency matrix and learned communities with the results presented in our paper. 

### Exploring what PAC has learned
- **combined_graph_analyses/full_graph_summary_hallmark.ipynb** Investigates the size and database distribution for each pathway community, and compares learned communities to hallmark founders **(Figure 3)**.
- **combined_graph_analyses/full_graph_summary_curated.ipynb** Evaluates the correspondence between the full community graph and curated categories for each database (both via NMI scores and visually with clustermaps).
- **combined_graph_analyses/interpret_communities_genes.ipynb** Evaluates the relevance of genes to each community.
- **combined_graph_analyses/query_new_pathway_against_communities.ipynb** Notebook providing a function for querying new pathways against the learned community (by adding the query gene set as a node to the graph and evaluating the modularity of different community assignments). We query all hallmark pathways for use in the 'full_graph_summary_curated' notebook above.
- **combined_graph_analyses/communities_kmers-and-hubness.ipynb** Identifies automatic labels for each community and also calculates hubness among members of each community subgraph.
- **combined_graph_analyses/Network-Graphs.ipynb** A notebook containing code to visualize learned networks. We visualize both a high-level view of communities, as well as sub-graphs for each individual community. **(Figure 4)** 
- **combined_graph_analyses/BreastCancerExample.ipynb** A notebook containing a biological example of how the PAC method clarifies findings of a differential expression analysis. **(Figure 5)**

