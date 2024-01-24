
# -- import packages: ---------------------------------------------------------
import ABCParse
import anndata
import numpy as np
import scipy

# -- set typing: --------------------------------------------------------------
from typing import List


class RunningQuantile(ABCParse.ABCParse):
    def __init__(self, n_bins: int = 50, *args, **kwargs):

        """
        n_bins: int, default = 50
        """

        self.__parse__(locals(), public=[None])

    @property
    def x(self):
        if not hasattr(self, "_x_cp"):
            self._x_cp = self._x.copy()[self.idx]
        return self._x_cp

    @property
    def y(self):
        if not hasattr(self, "_y_cp"):
            self._y_cp = self._y.copy()[self.idx]
        return self._y_cp

    @property
    def idx(self):
        return np.argsort(self._x)

    @property
    def xi(self):
        return self.x[0]

    @property
    def xf(self):
        return self.x[-1]

    @property
    def dx(self):
        return (self.xf - self.xi) / self._n_bins

    @property
    def pad(self):
        return self.dx / 2

    @property
    def x_out(self):
        if not hasattr(self, "_x_out"):
            self._x_out = np.linspace(
                self.xi + self.pad, self.xf - self.pad, self._n_bins
            )
        return self._x_out

    @property
    def y_out(self):
        if not hasattr(self, "_y_out"):
            self._y_out = np.zeros_like(self.x_out)
        return self._y_out

    def __call__(self, x, y, fit_percentile: float = 0.1):

        """
        Calculate the quantile of y in bins of x

        Parameters
        ----------
        x

        y

        fit_percentile: float, default = 0.1

        Returns
        -------
        x_out, y_out
        """

        self.__update__(locals(), public=[None])

        for i in range(len(self.x_out)):
            ix = np.nonzero(
                (x >= self.x_out[i] - self.pad) & (x < self.x_out[i] + self.pad)
            )[0]
            if len(ix) > 0:
                self.y_out[i] = np.percentile(y[ix], self._fit_percentile)
            elif i > 0:
                self.y_out[i] = y[i - 1]
            else:
                self.y_out[i] = np.nan

        return self.x_out, self.y_out

class VScores(ABCParse.ABCParse):

    """
    Calculate v-score (above-Poisson noise statistic) for genes in the input sparse counts matrix
    Return v-scores and other stats
    """

    def __init__(
        self,
        min_mean: float = 0,
        n_bins: int = 50,
        fit_percentile: float = 0.1,
        error_wt: int = 1,
        *args,
        **kwargs
    ):

        self.__parse__(locals(), public=[None])

    @property
    def running_quantile(self):
        if not hasattr(self, "_running_quantile"):
            self._running_quantile = RunningQuantile(self._n_bins)
        return self._running_quantile

    @property
    def n_cells(self):
        self._X.shape[0]

    @property
    def mu_gene(self):
        if not hasattr(self, "_mu_gene"):
            _mu_gene = self._X.mean(axis=0).A.squeeze()
            self._gene_ix = np.nonzero(_mu_gene > self._min_mean)[0]
            self._mu_gene = _mu_gene[self._gene_ix]
        return self._mu_gene

    @property
    def FF_gene(self):
        if not hasattr(self, "_FF_gene"):
            tmp = self._X[:, self._gene_ix]
            tmp.data **= 2
            var_gene = tmp.mean(axis=0).A.squeeze() - self.mu_gene**2
            del tmp
            self._FF_gene = var_gene / self.mu_gene

        return self._FF_gene

    @property
    def data_x(self):
        return np.log(self.mu_gene)

    @property
    def data_y(self):
        return np.log(self.FF_gene / self.mu_gene)

    def _filter_xy(self, data_x, data_y):
        x, y = self.running_quantile(data_x, data_y, self._fit_percentile)
        MASK = ~np.isnan(y)
        return x[MASK], y[MASK]

    def compute_abc(self, x, y):

        g_log = lambda input: np.log(input[1] * np.exp(-input[0]) + input[2])

        h, b = np.histogram(np.log(self.FF_gene[self.mu_gene > 0]), bins=200)
        b = b[:-1] + np.diff(b) / 2
        max_ix = np.argmax(h)
        c = np.max((np.exp(b[max_ix]), 1))
        error_func = lambda b2: np.sum(abs(g_log([x, c, b2]) - y) ** self._error_wt)
        b0 = 0.1
        b = scipy.optimize.fmin(func=error_func, x0=[b0], disp=False)
        a = c / (1 + b) - 1

        return a, b, c

    @property
    def v_scores(self):
        if not hasattr(self, "_v_scores"):
            self._v_scores = self.FF_gene / (
                (1 + self.a) * (1 + self.b) + self.b * self.mu_gene
            )
        return self._v_scores

    @property
    def CV_eff(self):
        if not hasattr(self, "_CV_eff"):
            self._CV_eff = np.sqrt((1 + self.a) * (1 + self.b) - 1)
        return self._CV_eff

    @property
    def CV_input(self):
        return np.sqrt(self.b)

    def __call__(self, X: scipy.sparse.csr_matrix, *args, **kwargs):

        """
        X:
            sparse expression matrix.
        """

        self.__update__(locals(), public=[None])

        x, y = self._filter_xy(self.data_x, self.data_y)
        self.a, self.b, self.c = self.compute_abc(x, y)

        return {
            "v_scores": self.v_scores,
            "CV_eff": self.CV_eff,
            "CV_input": self.CV_input,
            "gene_ix": self._gene_ix,
            "mu_gene": self.mu_gene,
            "FF_gene": self.FF_gene,
        }
    
class HighlyVariableGenes(ABCParse.ABCParse):
    def __init__(
        self,
        min_mean: float = 0,
        n_bins: int = 50,
        fit_percentile: float = 0.1,
        error_wt: int = 1,
        *args,
        **kwargs
    ):

        self.__parse__(locals(), public=[None])

    @property
    def base_idx(self):
        if len(self._base_idx) == 0:
            self._base_idx = np.arange(self._X.shape[0])
        return self._base_idx

    def _filtering(self):

        filter_idx = self.v_scores > 0
        v_scores = self.v_scores[filter_idx]
        gene_ix = self.gene_ix[filter_idx]
        mu_gene = self.mu_gene[filter_idx]
        FF_gene = self.FF_gene[filter_idx]

        self.min_vscore = np.percentile(v_scores, self._min_vscore_percentile)
        self._idx = (
            (self._X[:, gene_ix] >= self._min_counts).sum(0).A.squeeze()
            >= self._min_cells
        ) & (v_scores >= self.min_vscore)

    @property
    def hv_idx(self):
        return self.gene_ix[self._idx]

    def _process_v_score_stats(self):
        """"""
        self._VSCORES = VScores(
            min_mean=self._min_mean,
            n_bins=self._n_bins,
            fit_percentile=self._fit_percentile,
            error_wt=self._error_wt,
        )

        v_scores_output = self._VSCORES(self._X)
        for key, val in v_scores_output.items():
            setattr(self, key, val)

    def _process_outputs(self):

        self._adata.var[self._key_added] = self._adata.var.index.astype(int).isin(
            self.hv_idx
        )

        if self._return_idx:
            return self.hv_idx

    def __call__(
        self,
        adata: anndata.AnnData,
        key_added: str = "hv_gene",
        base_idx: List = [],
        min_vscore_percentile: float = 85,
        min_counts: float = 3,
        min_cells: float = 3,
        show_vscore_plot: bool = False,
        sample_name: str = "",
        return_idx: bool = False,
        *args,
        **kwargs
    ):

        """
        Filter genes by expression level and variability. Returns
        list of filtered gene indices. Annotates `adata.var` with
        "hv_genes" column.
        """

        self.__update__(locals(), public=[None], ignore=["adata"])

        self._adata = adata
        self._X = self._adata.X

        self._process_v_score_stats()
        self._filtering()
        return self._process_outputs()


# -- API-facing function: -----------------------------------------------------
def annotate_hv_genes(
    adata: anndata.AnnData,
    key_added: str = "hv_gene",
    base_idx: List = [],
    min_vscore_percentile: float = 75,
    min_counts: float = 3,
    min_cells: float = 3,
    show_vscore_plot: bool = False,
    sample_name: str = "",
    return_idx: bool = False,
    min_mean: float = 0,
    n_bins: int = 50,
    fit_percentile: float = 0.1,
    error_wt: int = 1,
    *args,
    **kwargs
):
    """
    adata: anndata.AnnData
        Annotated single-cell data object.

    key_added: str, default = "hv_gene"
        Key added to adata.var
        
    base_idx: List = [],
    min_vscore_percentile: float = 75,
    min_counts: float = 3,
    min_cells: float = 3,
    show_vscore_plot: bool = False,
    sample_name: str = "",
    return_idx: bool = False,
    min_mean: float = 0,
    n_bins: int = 50,
    fit_percentile: float = 0.1,
    error_wt: int = 1,
    """

    hv_gene_annot = HighlyVariableGenes(
        min_mean=min_mean,
        n_bins=n_bins,
        fit_percentile=fit_percentile,
        error_wt=error_wt,
    )

    return hv_gene_annot(
        adata=adata,
        key_added=key_added,
        base_idx=base_idx,
        min_vscore_percentile=min_vscore_percentile,
        min_counts=min_counts,
        min_cells=min_cells,
        show_vscore_plot=show_vscore_plot,
        sample_name=sample_name,
        return_idx=return_idx,
    )
