

# _make_peak_bed_from_multiome_features.py

__module_name__ = "_make_peak_bed_from_multiome_features.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu",])


# package imports #
# --------------- #
import os
import pandas as pd
import scipy.io as io
from scipy import sparse


def _make_peak_bed_from_multiome_features(multiome_dir, return_data=False):

    """
    Convert multiome feature.tsv to a file that only contains scATAC-seq peaks. 
    
    Parameters:
    -----------
    path
        path to feature.tsv
        type: str
    
    return_df
        defaut: False
        type: bool
        
    Returns:
    --------
    [ optional ] peaks
        pandas DataFrame of peaks
    """
    feature_df = pd.read_csv(
        os.path.join(multiome_dir, "features.tsv"),
        sep="\t",
        header=None,
        names=["ensg", "name", "feature_type", "chr", "start", "stop"],
    )

    peak_bed_filepath = os.path.join(multiome_dir, "peaks.bed")
    peak_df = feature_df.loc[feature_df.feature_type == "Peaks"][
        ["chr", "start", "stop"]
    ]
    peak_idx = peak_df.index.astype(int)
    peak_df = peak_df.reset_index(drop=True)
    peak_df.to_csv(peak_bed_filepath, sep="\t", header=None, index=False)

    mat = io.mmread(os.path.join(multiome_dir, "matrix.mtx")).toarray()[peak_idx]
    io.mmwrite(
        target=os.path.join(multiome_dir, "peak.matrix.mtx"), a=sparse.csr_matrix(mat)
    )

    if return_data:
        return feature_df, peak_df, mat