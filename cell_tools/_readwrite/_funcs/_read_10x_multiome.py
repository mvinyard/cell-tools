
# _read_10x_multiome.py

__module_name__ = "_read_10x_multiome.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu",])


# package imports #
# --------------- #
import anndata as a
import os

# local imports #
# ------------- #
from ._FeatureBarcodeMatrixModule import _FeatureBarcodeMatrix


def _read_10x_multiome(
    path,
    adata_filename="adata.h5ad",
    files=["matrix.mtx", "barcodes.tsv", "features.tsv"],
    obs_names=["cell_bc"],
    var_names=["ENSG", "gene_id", "omic", "Chromosome", "Start", "End"],
    set_obs_index_col=False,
    add_obs_id=False,
    write=True,
    silent=False,
):

    """
    Load 10x multiome from a 10x outs feature_bc_matrix directory. 
    
    Parameters:
    -----------
    path
        path to the directory containing the 10x outs for a feature barcode matrix. 
        type: str
    
    adata_filename
        filename to be saved within the path contianing the 10x outs in the formatted, compact .h5ad format. 
        type: str
        
    files
        files contained in `path` of the 10x outs. 
        type: str
        
    obs_names
        names of the adata.obs pandas.DataFrame column.
        type: list(str)
    
    var_names
        names of the adata.var pandas.DataFrame column.
        type: list(str)
    
    write
        If loading for the first time, determine if adata should be written as .h5ad for future l
        default: True
        type: bool
    
    silent
        Silence printed details
        default: False
        type: bool
    
    Returns:
    --------
    adata
    
    Notes:
    ------
    (1) Eventually should add support for reading from .h5 or other formats of 10x outs. 
    """

    data_10x = _FeatureBarcodeMatrix(path, adata_filename, files)

    if os.path.exists(data_10x.h5ad_path):
        print("Loading previously formatted 10x adata...\n")
        adata = a.read_h5ad(data_10x.h5ad_path)
    else:
        data_10x.load(obs_names, var_names)
        adata = data_10x.to_adata(write, silent, return_adata=True)
                
    if add_obs_id:
        for key, value in add_obs_id.items():
            adata.obs[key] = value
    if set_obs_index_col:
        adata.obs = adata.obs.set_index(set_obs_index_col)
    if not silent:
        print(adata, "\n")
    return adata