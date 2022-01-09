import numpy as np
import sklearn.preprocessing
import sklearn.decomposition
import umap


def _split_test_train(adata, column="Time point", test_criteria=4):

    train = np.full(adata.shape[0], True)
    train[adata.obs.loc[adata.obs[column] == test_criteria].index.astype(int)] = False
    X_train, X_test = adata.X[train], adata.X[~train]

    return X_train.toarray(), X_test.toarray(), train


class _DimensionalityReduction:
    def __init__(
        self,
        adata,
        split=False,
        column="Time point",
        test_criteria=4,
        n_pcs=50,
        n_neighbors=30,
        n_umap_components=2,
        umap_metric="euclidean",
    ):

        self._scaler = sklearn.preprocessing.StandardScaler()
        self._pca_reducer = sklearn.decomposition.PCA(n_components=n_pcs)
        self._umap_reducer = umap.UMAP(
            n_components=n_umap_components, metric=umap_metric, n_neighbors=n_neighbors
        )
        self.split = split
        self.adata = adata
        self.adata.obsm["X_scaled"] = np.zeros(adata.shape)
        self.adata.obsm["X_pca"] = np.zeros([adata.shape[0], n_pcs])
        self.adata.obsm["X_umap"] = np.zeros([adata.shape[0], n_umap_components])

        if not self.split:
            X = adata.X
        else:
            self.X_train, self.X_test, self.train = _split_test_train(adata, column, test_criteria)

    def scale(self):

        if not self.split:
            self.adata.obsm["X_scaled"] = self.X_scaled = self._scaler.fit_transform(X)
        else:
            self.adata.obsm["X_scaled"][self.train] = self.X_scaled_train = self._scaler.fit_transform(self.X_train)
            self.adata.obsm["X_scaled"][~self.train] = self.X_scaled_test = self._scaler.transform(self.X_test)

    def PCA(self, X=False):

        if not X:
            X = self.adata.obsm["X_scaled"]

        if not self.split:
            self.adata.obsm["X_pca"] = self._pca_reducer.fit_transform(X)
        else:
            self.adata.obsm["X_pca"][self.train] = self.X_pca_train = self._pca_reducer.fit_transform(self.X_scaled_train)
            self.adata.obsm["X_pca"][~self.train] = self.X_pca_test = self._pca_reducer.transform(self.X_scaled_test)

    def UMAP(
        self, X=False, n_neighbors=False, n_umap_components=False, umap_metrics=False
    ):

        """"""

        if not X:
            X = self.adata.obsm["X_pca"]

        if not self.split:
            self.adata.obsm["X_umap"] = self._umap_reducer.fit_transform(X)
        else:
            self.adata.obsm["X_umap"][self.train] = self.X_umap_train = self._umap_reducer.fit_transform(self.X_pca_train)
            self.adata.obsm["X_umap"][~self.train] = self.X_umap_test = self._umap_reducer.transform(self.X_pca_test)


def _reduce_adata(
    adata,
    split,
    column,
    test_criteria,
    n_pcs=50,
    n_neighbors=30,
    n_umap_components=2,
    umap_metric="euclidean",
):

    """"""

    transformer = _DimensionalityReduction(
        adata,
        split,
        column,
        test_criteria,
        n_pcs,
        n_neighbors,
        n_umap_components,
        umap_metric,
    )

    print("Scaling adata.X using StandardScaler\n")
    transformer.scale()

    print("PCA: n_components={}\n".format(n_pcs))
    transformer.PCA()

    print(
        "UMAP: n_components={}, n_neighbors={}, umap_metrics={}\n".format(
            n_umap_components, n_neighbors, umap_metric
        )
    )
    transformer.UMAP()

    return transformer