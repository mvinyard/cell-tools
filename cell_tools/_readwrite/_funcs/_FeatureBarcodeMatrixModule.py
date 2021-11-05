
# _FeatureBarcodeMatrixModule.py
__module_name__ = "_FeatureBarcodeMatrixModule.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu",])

# package imports #
# --------------- #
import anndata as a
import pandas as pd
import scipy.io
import os

def _strip_file_extension(filepath):
    return ".".join(filepath.split(".")[:-1])


def _create_feature_bc_matrix_path(
    path, files=["matrix.mtx", "barcodes.tsv", "features.tsv"]
):

    """"""

    PathDict = {}

    for file in files:
        PathDict[_strip_file_extension(file)] = os.path.join(path, file)

    return PathDict

class _FeatureBarcodeMatrix:
    def __init__(
        self,
        path,
        adata_filename="adata.h5ad",
        files=["matrix.mtx", "barcodes.tsv", "features.tsv"],
    ):

        """"""

        self.PathDict = _create_feature_bc_matrix_path(path, files)
        self.h5ad_path = os.path.join(path, adata_filename)

    def load(
        self,
        obs_names=["cell_bc"],
        var_names=["ENSG", "gene_id", "omic", "Chromosome", "Start", "End"],
    ):

        self.X = scipy.io.mmread(self.PathDict["matrix"])
        self.obs = pd.read_csv(
            self.PathDict["barcodes"], header=None, sep="\t", names=obs_names
        )
        self.var = pd.read_csv(
            self.PathDict["features"], header=None, sep="\t", names=var_names
        )

    def to_adata(self, write=True, silent=False, return_adata=False):

        self.adata = a.AnnData(self.X.toarray()).T

        self.adata.obs = self.obs
        self.adata.var = self.var

        if not silent:
            print(self.adata)

        if write:
            print("\nWriting AnnData to:\n\t {}".format(self.h5ad_path))
            self.adata.write_h5ad(self.h5ad_path)
            
        if return_adata:
            return self.adata