
import pandas as pd

def _write_adata_with_warning_reduction(adata, path):
    
    
    if adata.raw:
        try:
            for col in adata.raw.var.columns:
                adata.raw.var[col] = pd.Categorical(adata.raw.var[col])
        except:
            pass
        
        try:
            for col in adata.raw.obs.columns:
                adata.raw.obs[col] = pd.Categorical(adata.raw.obs[col])
        except:
            pass
    
    for col in adata.var.columns:
        adata.var[col] = pd.Categorical(adata.var[col])
    
    for col in adata.obs.columns:
        adata.obs[col] = pd.Categorical(adata.obs[col])
    
    adata.var = adata.var.reset_index(drop=True)
    adata.var.index = adata.var.index.astype(str)

    adata.obs = adata.obs.reset_index()
    adata.obs.index = adata.obs.index.astype(str)
    
    adata.write_h5ad(path)