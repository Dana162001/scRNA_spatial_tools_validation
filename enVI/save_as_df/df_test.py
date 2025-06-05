import scenvi 
import scanpy as sc
import anndata


sc_data = sc.read_h5ad('/work/rwth1209/dana_projects/test_new_imputation_tools/enVI/data_preprocessing/scrna_ref_all_raw_array.h5ad')
st_data = sc.read_h5ad('/work/rwth1209/dana_projects/test_new_imputation_tools/enVI/data_preprocessing/withhelded_array_raw_483_genes.h5ad') 

envi_model = scenvi.ENVI(spatial_data = st_data, sc_data = sc_data)

envi_model.train()
envi_model.impute_genes()
envi_model.infer_niche()

 

st_data.obsm['envi_latent'] = envi_model.spatial_data.obsm['envi_latent']
st_data.obsm['COVET'] = envi_model.spatial_data.obsm['COVET']
st_data.obsm['COVET_SQRT'] = envi_model.spatial_data.obsm['COVET_SQRT']
#st_data.uns['COVET_genes'] =  envi_model.CovGenes
st_data.obsm['imputation'] = envi_model.spatial_data.obsm['imputation']

# sc_data.obsm['envi_latent'] = envi_model.sc_data.obsm['envi_latent']
# sc_data.obsm['COVET'] = envi_model.sc_data.obsm['COVET']
# sc_data.obsm['COVET_SQRT'] = envi_model.sc_data.obsm['COVET_SQRT']
#sc_data.uns['COVET_genes'] =  envi_model.CovGenes


imputed_adata = sc.AnnData(st_data.obsm['imputation'])
imputed_adata.obsm['spatial'] = st_data.obsm['spatial']
imputed_adata.obsm['envi_latent'] = st_data.obsm['envi_latent']
imputed_adata.obsm['COVET'] = st_data.obsm['COVET']
imputed_adata.obsm['COVET_SQRT'] = st_data.obsm['COVET_SQRT']
imputed_adata.obsm['imputation'] = st_data.obsm['imputation']
#imputed_adata.obsm['COVET_genes'] = st_data.obsm['COVET_genes']

imputed_adata.obs = st_data.obs
imputed_adata

imputed_adata.write('imputed_adata_15_genes_all.h5ad')
#sc_data.write('output_scrna_data.h5ad')