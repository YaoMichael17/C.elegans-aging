
import scanpy as sc
import pandas as pd

# Load the aging worm dataset
adata = sc.read_h5ad("/Users/cmdb/QB Project/ad_worm_aging.h5ad")

print("Initial shape:", adata.shape)

# Define stress and neuron markers
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
stress_genes_present = [g for g in stress_genes if g in adata.var_names]
# print("Stress genes present:", stress_genes_present)

if not stress_genes_present:
    print("no stress genes in the dataset")

neuron_markers = ['unc-47', 'unc-25']  # used to detect neuron clusters

# Store raw data for future use
adata.raw = adata.copy()

# Normalize counts and log-transform
sc.pp.normalize_total(adata, target_sum = 1e4)
sc.pp.log1p(adata)

# Scale data + run PCA + compute neighbors graph
sc.pp.scale(adata)
sc.tl.pca(adata, svd_solver = 'arpack') # Scanpy .pca function is epic 
sc.pp.neighbors(adata, n_neighbors = 15, n_pcs = 30) # Input for leidan clustering below

# Leiden clustering to identify cell clusters
sc.tl.leiden(adata, resolution = 0.5)
print("Clusters:", adata.obs['leiden'].cat.categories)

# Identify neuron clusters based on neuron markers above and their expression (Looking for high expression above 0.1 [somewhat arbitrary but seems to filter well])
cluster_expr = adata.to_df()[neuron_markers].groupby(adata.obs['leiden']).mean()
neuron_clusters = cluster_expr[(cluster_expr > 0.1).any(axis = 1)].index.tolist()  # threshold can be adjusted
print("Neuron clusters detected:", neuron_clusters)

# Subset dataset to neuron clusters only
neuron_adata = adata[adata.obs['leiden'].isin(neuron_clusters)].copy()
print("Neuron-only dataset shape:", neuron_adata.shape)

# Compute average stress gene expression per timepoint
expr = neuron_adata.raw.to_adata().to_df()[stress_genes_present].copy()
expr['timepoint'] = neuron_adata.obs['timepoint'].values
neuron_means = expr.groupby("timepoint").mean().reset_index()

# PLEASE keep timepoints order
FINAL_DAY_ORDER = ["d1", "d3", "d5", "d8", "d11", "d15"]
available_order = [tp for tp in FINAL_DAY_ORDER if tp in neuron_means['timepoint'].unique()]
neuron_means = neuron_means.set_index("timepoint").loc[available_order].reset_index()


print("Neuron stress-gene means:\n", neuron_means.head())
#print("Neuron stress-gene means:\n", neuron_adata.head())

# Save results to .csv file
outfile = "/Users/cmdb/QB Project/neuronal_GE_timepoint_averages.csv"
neuron_means.to_csv(outfile, index = False)

print("Saved neuron stress-gene expression to:", outfile)
