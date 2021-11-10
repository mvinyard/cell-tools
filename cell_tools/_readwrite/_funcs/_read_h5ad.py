
import anndata as a

def _read_h5ad(h5ad_path, silent=False):
    
    """Wraps anndata.read_h5ad"""
    
    adata = a.read_h5ad(h5ad_path)
    
    if not silent:
        print(adata)
        
    return adata