
# import packages # 
# --------------- #
import licorice
from sklearn.decomposition import PCA

# local imports # 
# ------------- # 
from ..._utilities._funcs._scale import _scale_anndata


def _pca(X, n_components=50, **kwargs):
    
    """Principle Component Analysis"""
    
    pca = PCA(n_components, **kwargs)
    X_pca = pca.fit_transform(X)
    
    return pca, X_pca

def _pca_anndata(adata, n_components=50, return_adata=False, **kwargs):
    
    try:
        adata.layers['X_scaled']
    except:
        print("\n Scaling adata.X and adding the layer: adata.layers['{}']".format(licorice.font_format("X_scaled", ['BOLD'])))
        adata = _scale_anndata(adata, return_adata=True)
        
    adata.uns['pca'], adata.obsm['X_pca'] = _pca(adata.layers['X_scaled'], n_components, **kwargs)
    
    if return_adata:
        return adata