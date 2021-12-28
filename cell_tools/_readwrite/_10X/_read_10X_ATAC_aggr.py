
__module_name__ = "_read_10x_ATAC_aggr.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu",])


# import packages #
# --------------- #
import episcanpy
import os
import pandas as pd


def _read_10X_aggr_atac_h5_to_adata(path, file="filtered_peak_bc_matrix.h5"):
    
    """"""
    
    filepath = os.path.join(path, "outs", file)
    adata = episcanpy.read_h5_atac(filename=filepath)
    
    return adata

def _load_library_ids(aggr_path):
    
    """
    
    Returns:
    --------
    library_ids
        type: pandas.core.frame.DataFrame
        
    (1) Read the aggregation csv
    
    (2) assumes path: outs/aggregation_csv.csv
    """
    
    aggr_csv = pd.read_csv(os.path.join(aggr_path, "outs/aggregation_csv.csv"))
    library_ids = pd.DataFrame(aggr_csv['library_id'])
    library_ids['library_idx'] = library_ids.reset_index()['index']  + 1
    library_ids['library_idx'] = library_ids['library_idx'].astype(str)
    return library_ids

def _annotate_library(adata):
        
        
    lib_ids = adata.uns['lib_ids']
    
    tmp = adata.obs.reset_index()
    tmp = tmp['index'].str.split('-', expand=True)
    tmp.columns = ['cell_barcode', 'library_idx']
    tmp['library_idx'] = tmp['library_idx'].astype(str)
    obs_df = tmp.merge(lib_ids, on='library_idx')
    obs_df.index = obs_df.index.astype(str)
    
    del tmp
    
    return obs_df

def _format_10X_atac_adata_var(adata):
    
    tmp_var = adata.var.reset_index()
    tmp_var = tmp_var.drop('name', axis=1)
    tmp_var = tmp_var.rename({'id':'coordinates'}, axis=1)
    tmp_var.index = tmp_var.index.astype(str)
    
    return adata.var

def _read_10X_ATAC_aggr(aggr_path, silent=False):
    
    """
    Read the resulting outfiles of the 10X cellranger-atac aggr function. 
    
    Parameters:
    -----------
    aggr_path
        path to outfiles of 10X cellranger-atac aggr
        type: str
    
    silent
        type: bool
        default: False
    
    Returns:
    --------
    adata
    
    Notes:
    ------
    (1) 10X cellranger-atac aggr function documentation here: 
        https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/using/aggr
    """
    
    adata = _read_10X_aggr_atac_h5_to_adata(aggr_path)
    adata.uns['lib_ids'] = _load_library_ids(aggr_path)
    adata.obs = _annotate_library(adata)
    adata.var = _format_10X_atac_adata_var(adata)
    
    if not silent:
        print(adata)
    
    return adata