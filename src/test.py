import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys

stress_genes = ["hsp-16.2", "hsp-70", "hsp-4", "gst-4", "ctl-1", "ctl-2", "atg-18", "lgg-1", "bec-1"]

adata = sc.read_h5ad("/Users/cmdb/QB Project/ad_worm_aging.h5ad")

# Identify which stress genes are present in the dataset
all_genes = list(adata.var_names)
stress_genes_in_data = [gene for gene in stress_genes if gene in all_genes]

# Print out which genes were found and which were not
print("Stress genes found in the dataset:", stress_genes_in_data)

# Subset the data to include only the stress genes
# `adata[:, list_of_genes]` filters by genes
adata_stress_genes = adata[:, stress_genes_in_data].copy()

# Access the expression data
expression_data = adata_stress_genes.X.toarray()

# Pandas DataFrame for readability
expression_df = pd.DataFrame(
    data = expression_data,
    index = adata_stress_genes.obs_names,
    columns = adata_stress_genes.var_names
)

# SHOW ALL EXPRESSION DATA
print("\nExpression data for stress genes:")
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None) # Optional: Adjusts the width of the display
# print(expression_df)

expression_df.to_csv("/Users/cmdb/QB Project/expression_df.csv", index = False)

#------------------------------------------------------------------------------------
# Erroneous
average_expression = expression_df.mean()
print("\nAverage expression per stress gene:")
print(average_expression)