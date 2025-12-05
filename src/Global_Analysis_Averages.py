
import scanpy as sc
import pandas as pd
import numpy as np

# Load single cell RNA-seq dataset
adata = sc.read_h5ad("/Users/cmdb/QB Project/ad_worm_aging.h5ad")  # AnnData object with all cells

# Define stress hallmark genes of interest
stress_genes = [
    "hsp-16.2", "hsp-70", "hsp-4", "gst-4", "ctl-1", "ctl-2", "atg-18", "lgg-1",
    "bec-1", "sip-1", "myo-5", "gpx-7", "lmd-3", "ama-1", "act-1", "cdc-42",
    "Mms-19", "Exo-1", "Ung-1", "Polh-1", "Ddb-1", "pms-2", "Rnf-121",
    "Rnf-113", "trul-1", "wwp-1", "rbx-1", "cul-4", "fsn-1", "uba-1",
    "hsp-3", "hsp-17", "hsp-43", "hsp-1", "hsp-60",
    # Newly added genes
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

# Keep only genes in the dataset
stress_genes_in_data = [gene for gene in stress_genes if gene in adata.var_names]  # Filter genes to those available

# Subset AnnData to only stress genes
adata_stress = adata[:, stress_genes_in_data].copy()  # New AnnData with only stress genes

# Normalize total counts per cell to 10,000 to correct for sequencing depth differences
sc.pp.normalize_total(adata_stress, target_sum = 1e4)

# Log-transform counts to reduce skew and stabilize variance
sc.pp.log1p(adata_stress)

# Convert sparse expression matrix to dense array for Pandas
expression_data = adata_stress.X.toarray()  # Rows = cells, Columns = genes

# Create Pandas DataFrame (ChatGPT help)
expression_df = pd.DataFrame(
    data = expression_data,
    index = adata_stress.obs_names,
    columns = adata_stress.var_names
)

# Add timepoint information
expression_df['timepoint'] = adata_stress.obs['timepoint'].values

# Compute average expression per gene for each timepoint
average_expression_df = expression_df.groupby('timepoint')[stress_genes_in_data].mean().reset_index()

# Save the tidy CSV for R visualization
average_expression_df.to_csv("/Users/cmdb/QB Project/global_GE_timepoint_averages.csv", index = False)

print(average_expression_df)  # Quick check of average expression per gene at each timepoint
