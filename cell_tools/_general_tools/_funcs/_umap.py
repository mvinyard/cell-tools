
# import packages # 
# --------------- # 
import umap

# local imports # 
# ------------- # 
from ..._utilities._funcs._scale import _scale_anndata


def _umap(X, n_components=2, n_neighbors=20, **kwargs):
    
    """UMAP"""
    
    umap_model = umap.UMAP(n_components=n_components, n_neighbors=n_neighbors, **kwargs)
    X_umap = umap_model.fit_transform(X)
    umap_model.verbose = False

    
    return umap_model, X_umap

def _umap_anndata(adata, use="X_pca", n_components=2, return_adata=False, **kwargs):
    
    
        
    adata.uns['umap'], adata.obsm['X_umap'] = _umap(adata.obsm[use], n_components, n_neighbors=20, **kwargs)
    
    if return_adata:
        return adata
