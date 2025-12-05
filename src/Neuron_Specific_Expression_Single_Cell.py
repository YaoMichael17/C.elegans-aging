
import scanpy as sc
import pandas as pd
import numpy as np

# read the dataset
adata = sc.read_h5ad("/Users/cmdb/QB Project/ad_worm_aging.h5ad")
print("Initial shape:", adata.shape)

# list of stress genes and neuron markers we care about
stress_genes = [
    "hsp-16.2", "hsp-70", "hsp-4", "gst-4", "ctl-1", "ctl-2", "atg-18", "lgg-1",
    "bec-1", "sip-1", "myo-5", "gpx-7", "lmd-3", "ama-1", "act-1", "cdc-42",
    "Mms-19", "Exo-1", "Ung-1", "Polh-1", "Ddb-1", "pms-2", "Rnf-121",
    "Rnf-113", "trul-1", "wwp-1", "rbx-1", "cul-4", "fsn-1", "uba-1",
    "hsp-3", "hsp-17", "hsp-43", "hsp-1", "hsp-60",
    "atfs-1", "kreg-1", "F53A9.6", "F22H10.2", "fbxa-60", "gpx-1", "cdr-4",
    "R12E2.13", "fipr-24", "mev-1", "fbxa-59", "fnci-1", "nanp-1",
    "F45E1.4", "Y73B6BL.14", "F42A6.5", "nse-4", "ubql-1", "xpc-1",
    "chaf-2", "nse-1", "sti-1", "jnk-1", "natc-1", "rnp-9", "hsp-16.48",
    "eif-2bdelta", "scl-1", "prdx-2", "tag-56", "gst-10", "F43C11.7",
    "tiar-2", "tiar-1", "ife-2", "tdp-1", "hsp-110", "pept-1", "kgb-1",
    "che-11", "xbp-1", "hcf-1", "nlp-7", "clk-2", "hpk-1", "skpo-1",
    "parp-1", "frh-1", "xpa-1", "mek-1", "bar-1", "pvl-1", "spy-1",
    "cep-1", "elt-3", "fkh-9", "apy-1", "crt-1", "idpc-5", "idpc-2",
    "F41B4.3", "ech-3", "clec-121", "irld-6", "sdz-35", "F22E5.6",
    "lips-11", "mpdu-1"
]

stress_genes_present = [g for g in stress_genes if g in adata.var_names]  # keep only ones actually in dataset
if not stress_genes_present:
    raise ValueError("None of the stress genes are in the dataset!")
print("Stress genes present:", stress_genes_present)

neuron_markers = ['unc-47', 'unc-25']  # markers to find neuron clusters

# make a copy to keep raw counts around
adata.raw = adata.copy()

# normalize, log transform, scale, PCA, neighbors (ChatGPT help)
sc.pp.normalize_total(adata, target_sum = 1e4)
sc.pp.log1p(adata)
sc.pp.scale(adata)
sc.tl.pca(adata, svd_solver = 'arpack')
sc.pp.neighbors(adata, n_neighbors = 15, n_pcs = 30)

# clustering cells
sc.tl.leiden(adata, resolution = 0.5)

# figure out which clusters are neurons
cluster_expr = adata[:, neuron_markers].to_df().groupby(adata.obs['leiden']).mean()
neuron_clusters = cluster_expr[(cluster_expr > 0.1).any(axis=1)].index.tolist()
neuron_adata = adata[adata.obs['leiden'].isin(neuron_clusters)].copy()  # keep only neurons

# ORDER TIMEPOINTS CORRECTLY
FINAL_DAY_ORDER = ["d1", "d3", "d5", "d8", "d11", "d15"]
available_order = [tp for tp in FINAL_DAY_ORDER if tp in neuron_adata.obs['timepoint'].unique()]
neuron_adata = neuron_adata[neuron_adata.obs['timepoint'].isin(available_order)].copy()

#ChatGPT help for DE analysis with sp
# run DE analysis between d11 vs d1 on all neuron genes
sc.tl.rank_genes_groups(
    neuron_adata,
    groupby='timepoint',
    groups=['d11'],
    reference='d1',
    method='wilcoxon'
)

# pull out only our stress genes and convert log fold change to log2
rg = neuron_adata.uns['rank_genes_groups']
all_genes = pd.DataFrame({
    'names': rg['names']['d11'],
    'logfoldchanges': np.array(rg['logfoldchanges']['d11']) / np.log(2),  # ln -> log2
    'pvals': rg['pvals']['d11'],
    'pvals_adj': rg['pvals_adj']['d11'],
    'scores': rg['scores']['d11']
})
# Above code is to fix R misreading the uploaded .csv file

# filter for only the stress genes
de_results = all_genes[all_genes['names'].isin(stress_genes_present)].reset_index(drop=True)
de_results['logfoldchanges'] = de_results['logfoldchanges'].fillna(0)  # fill any NaNs with 0

# save to CSV
de_results_file = "/Users/cmdb/QB Project/neuronal_DE_stress_d11_vs_d1_scanpy.csv"
de_results.to_csv(de_results_file, index = False)
print("Saved neuron-specific stress gene DE results to:", de_results_file)
