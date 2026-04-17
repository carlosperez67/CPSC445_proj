import scanpy as sc
import os
import anndata as ad

path = "samples_batch1/"
files = [f for f in os.listdir(path) if f.endswith('_cleaned.h5ad')]

adatas = []
for f in files:
    temp = sc.read_h5ad("samples_batch1/" + f)
    temp.obs['sample_id'] = f.replace('_cleaned.h5ad', '')
    adatas.append(temp)

combined = ad.concat(adatas, label="batch", join="inner")
combined.write_h5ad("samples_merged.h5ad")


adata = combined
adata.raw = adata
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

sc.pp.highly_variable_genes(
    adata,
    n_top_genes=3000,
    batch_key='sample_id',
    flavor = "seurat"
)

adata_hvg = adata[:, adata.var['highly_variable']].copy()
adata_hvg.obs['Condition'] = adata_hvg.obs['sample_id'].str[0]
adata_hvg.obs['Condition'] = adata_hvg.obs['Condition'].map({
    'K': 'Disease',
    'P': 'Control'
})
adata_hvg.obs['Condition'].value_counts()
adata_hvg.X = adata_hvg.X.tocsc()
adata_hvg.write("signal_prepared.h5ad")