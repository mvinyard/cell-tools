
# _RankFeaturesModule.py
__module_name__ = "_RankFeaturesModule.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu",])

# package imports #
# --------------- #
import scanpy as sc, episcanpy as epi


from ._supporting_functions._get_score_df import _get_score_df
from ._supporting_functions._run_episcanpy_group_rank_plots import _run_episcanpy_group_rank_plots
from ._supporting_functions._write_differential_features_to_excel import _write_differential_features_to_excel


class _RankFeatures:
    def __init__(
        self, adata, groupby, n_features=100, out_prefix=False, n_plot_features=5
    ):

        """Runs rank episcanpy.tl.rank_features and sc.tl.dendrogram"""

        self.out_prefix = out_prefix
        self._adata = adata
        self._adata.raw = adata
        self.groupby = groupby
        self.n_features = n_features
        self.n_plot_features = n_plot_features
        self.group_key = "by_{}".format(groupby)

        epi.tl.rank_features(
            self._adata,
            groupby=groupby,
            omic="ATAC",
            key_added=self.group_key,
            n_features=n_features,
        )
        sc.tl.dendrogram(self._adata, groupby=groupby)

    def rank_plots(self):

        _run_episcanpy_group_rank_plots(
            self._adata,
            self.groupby,
            self.group_key,
            self.out_prefix,
            self.n_plot_features,
        )

    def get_score_df(self):

        self.score_df = _get_score_df(self._adata, self.group_key)
        
        
    def to_excel(self):

        """Write top differential features (lfc) to Excel Workbook."""
    
        _write_differential_features_to_excel(self._adata, self.score_df, self.out_prefix, self.n_features, self.groupby)
        
        