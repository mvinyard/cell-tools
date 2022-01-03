import os
import pandas as pd


from ._read_h5 import _read_h5
from ._write_adata_with_warning_reduction import _write_adata_with_warning_reduction


def _load_10X_aggregation_csv(aggr_10x_path):

    """Modifies columns / library_indices for downstream integration with AnnData"""

    aggr_csv_path = os.path.join(aggr_10x_path, "aggregation_csv.csv")
    aggr_csv = pd.read_csv(aggr_csv_path).reset_index()
    aggr_csv.columns = ["library_index", "sample", "fragments", "outs"]
    aggr_csv["library_index"] = (aggr_csv["library_index"].astype(int) + 1).astype(str)

    return aggr_csv


def _split_aggr_library_cell_barcodes(adata, barcodes="barcodes"):

    """"""

    return adata.obs[barcodes].str.split("-", expand=True)[1]


def _read_10x_aggr(path, write_h5ad="aggr.10x.adata.h5ad", silent=False):

    """
    Creates AnnData object from h5 file. Annotates with library id's.


    """

    adata = _read_h5(os.path.join(path, "outs/filtered_peak_bc_matrix.h5"), silent=True)

    aggr_csv = _load_10X_aggregation_csv(os.path.join(path, "outs/"))
    adata.uns["10X_aggr_meta"] = aggr_csv[["fragments", "outs"]]
    aggr_csv = aggr_csv[["library_index", "sample"]]

    adata.obs["library_index"] = _split_aggr_library_cell_barcodes(
        adata, barcodes="barcodes"
    )

    tmp_obs = adata.obs.merge(aggr_csv, on="library_index", how="left")
    tmp_obs.index = tmp_obs.index.astype(str)
    adata.obs = tmp_obs
    del tmp_obs
        
    if write_h5ad:
        _write_adata_with_warning_reduction(adata, write_h5ad)
        
    if not silent:
        print(adata)
    
    return adata