
# general tools __init__.py

__module_name__ = "__init__.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu",])

from ._funcs._Harmony._HarmonizeModule import _Harmonize as Harmonize
from ._funcs._pca import _pca_anndata as PCA
from ._funcs._umap import _umap_anndata as UMAP

# these less common, simpler APIs for PCA and UMAP can also be used directly but they are hidden to prioritize the AnnData-compatible functions
from ._funcs._umap import _umap
from ._funcs._pca import _pca
