
# import packages # 
# --------------- # 
import umap

# local imports # 
# ------------- # 
from ..._utilities._funcs._scale import _scale_anndata


def _umap(X, n_components=2, **kwargs):
    
    """Principle Component Analysis"""
    
    umap_model = umap.UMAP(n_components, **kwargs)
    X_umap = umap_model.fit_transform(X)
    
    return umap_model, X_umap

def _umap_anndata(adata, use="X_pca", n_components=2, return_adata=False, **kwargs):
    
    
        
    adata.uns['umap'], adata.obsm['X_umap'] = _umap(adata.obsm[use], n_components, **kwargs)
    
    if return_adata:
        return adata
