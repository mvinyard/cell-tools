
# _HarmonizeModule.py

__module_name__ = "_HarmonizeModule.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu",])


# package imports #
# --------------- #
# from harmony import harmonize
from sklearn.preprocessing import StandardScaler
import umap

# local imports #
# ------------- #
from ._supporting_funcs._plot_umap import _plot_umap


class _Harmonize:
    def __init__(self, title=False):

        """Integrate data using Harmony algorithm. Wrapped from `harmony.harmonize`."""

        self.title = title
        self.scaler = StandardScaler()

    def run(self, adata, batch_key, use_key="X_pca"):

        """"""

        self.adata = adata
        self.batch_key = batch_key
        self.use_key = use_key
        self.harmony_key = "{}_harmony".format(self.use_key)
        self.adata.obsm[self.harmony_key] = harmonize(
            self.adata.obsm[use_key], self.adata.obs, batch_key=batch_key
        )

    def umap(self, n_components=2):

        """"""

        self.umap_key = "X_umap_harmony_{}".format(self.use_key.strip("X_"))

        self.reducer = umap.UMAP(n_components=n_components)
        harmonized_input = self.adata.obsm[self.harmony_key]
        scaled_input = self.scaler.fit_transform(harmonized_input)

        self.adata.obsm[self.umap_key] = self.reducer.fit_transform(scaled_input)

    def plot(self, plot_by, colors_dict=False):

        self.plot_by = plot_by
        self.fig = _plot_umap(self.adata, self.umap_key, plot_by, colors_dict)
        self.fig.savefig(
            "./figures/{}.harmonize_on_{}.{}.pdf".format(
                self.title, self.batch_key, plot_by
            )
        )

    def save(self):
        """"""

        self.h5ad_path = "{}.harmonized_on_{}.h5ad".format(self.title, self.batch_key)
        self.adata.write_h5ad(self.h5ad_path)