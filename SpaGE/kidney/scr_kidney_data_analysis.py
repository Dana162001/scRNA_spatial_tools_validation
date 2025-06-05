# %%
import scanpy as sc
import numpy as np
import pandas as pd
import loompy
import matplotlib.pyplot as plt
import scipy.stats as st
import warnings
warnings.filterwarnings('ignore')

# %%
spage = sc.read_h5ad('/work/rwth1209/dana_projects/test_new_imputation_tools/SpaGE/kidney/spage_3000_genes_kidney.h5ad')
original = sc.read_h5ad('/work/rwth1209/projects/Merfish_QC/500_kidney_26.2/preprocessed.h5ad')
scref = sc.read_h5ad('/work/rwth1209/dana_projects/test_new_imputation_tools/enVI/data_preprocessing/kidney/scref_3000hvg_raw_arr.h5ad')
sp = sc.read_h5ad('/work/rwth1209/dana_projects/test_new_imputation_tools/enVI/data_preprocessing/kidney/withhelded_kidney_array_raw_474_genes.h5ad')
scd = scref

# %%
spage

# %%
#original.X = original.layers["counts"].copy()
original = original[spage.obs_names, :]
sc.pp.filter_cells(original, min_counts=5)
sc.pp.filter_genes(original, min_cells=5)
sc.pp.normalize_total(original, target_sum=1e4)
sc.pp.log1p(original)
original

# %%
spage = spage[spage.obs_names, :]
sc.pp.filter_cells(spage, min_counts=5)
sc.pp.filter_genes(spage, min_cells=5)
sc.pp.normalize_total(spage, target_sum=1e4)
sc.pp.log1p(spage)
spage

# %%
original

# %%
original.obsm["spatial"]

# %%
spage.uns['spatial'] = original.obsm["spatial"]

# %%
original.X.max()

# %%
spage.X.max()

# %%
scd.X.max()

# %%
genes_to_check = ['NPHS1','NPHS2','WT1','EMCN','PECAM1','COL1A1','COL3A1','AQP2','CFH','LRP2','CUBN','TACSTD2','SLC8A1','HSD11B2','CALB1',]

# %%
# for gene in genes_to_check:
#     sc.pl.embedding(original, basis="spatial", color=gene, title=f"Original {gene}",vmax="p99")
#     sc.pl.embedding(spage, basis="spatial", color=gene, title=f"Spage {gene}",vmax="p99")

# # %%
# for gene in genes_to_check:
#     sc.pl.dotplot(original, gene, groupby="tacco",title=f"Original {gene}")
#     sc.pl.dotplot(spage, gene, groupby="tacco",title=f"Spage {gene}")

# %%
spage.write_h5ad("spage_3000_genes_norm_kidney.h5ad")

# %% [markdown]
# # Calculate Spearman correlation from two matrices

# %%
from scipy.stats import spearmanr
from sklearn.metrics import mean_squared_error

def spearman_correlation(a, b):
    a = a.toarray().flatten()
    b = b.toarray().flatten()
    print(a.shape, b.shape)
    print(mean_squared_error(a, b))
    return spearmanr(a, b)

# %%
for gene in genes_to_check:
    print(gene)
    print(spearman_correlation(original[:,gene].X, spage[:,gene].X))

# %% [markdown]
# # UMAP

# %%

# %%
spage.var

# %%
spatial_genes = list(sp.var_names)
len(spatial_genes)
genes_to_keep = [gene for gene in spage.var_names if gene not in spatial_genes]
spage_filtered = spage[:, genes_to_keep]

# %%
len(spatial_genes)

# %%
len(spage_filtered)

# %%
spage.X.max()

# %%
sc.pp.neighbors(spage)

# %%
sc.tl.umap(spage)

# %%
# run as a job!!!
sc.tl.leiden(spage)

# %%
sc.pp.neighbors(spage, n_pcs=30)
sc.tl.leiden(spage, key_added="leiden_res0_25", resolution=0.25)
sc.pl.umap(
    spage,
    color=["leiden_res0_25"],
    legend_loc="on data", save='spage_leiden_res0_25.png')

sc.tl.leiden(spage, key_added="leiden_res0_5", resolution=0.5)
sc.pl.umap(
    spage,
    color=["leiden_res0_5"],
    legend_loc="on data", save='spage_leiden_res0_5.png')
sc.pl.umap(
    spage,
    color=["tacco"],
    legend_loc="right margin", save='spage_cell_type.png')

# %%
spage.write_h5ad("spage_3000_genes_norm.h5ad")

# %%
sc.tl.umap(original)

# %%
original

# %%
sc.pl.umap(
    original,
    color=["leiden_0.5"],
    legend_loc="on data", save='original_leiden_0.5.png')

# %%
sc.pl.umap(
    original,
    color=["tacco"],
    legend_loc="on data",save='original_cell_type.png')


