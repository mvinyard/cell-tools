
import scanpy as sc
import numpy as np
import licorice
import matplotlib.pyplot as plt
import vinplots


def _build_plot(geneset_score_key):
    fig = vinplots.Plot()
    fig.construct(
        nplots=5,
        ncols=5,
        width_ratios=[2, 2, 2, 2, 0.05],
        wspace=0.15,
        figsize_width=1,
    )
    fig.modify_spines(
        ax="all",
        spines_to_delete=["top", "right"],
        spines_to_move=["bottom", "left"],
        spines_positioning_amount=15,
    )

    axes = []
    for i in range(len(fig.AxesDict[0])):
        axes.append(fig.AxesDict[0][i])
    titles = ["Unfiltered", "Filtered"]

    axes[1].set_title("Ranked {}".format(geneset_score_key))
    axes[1].set_xlabel("Cells")
    axes[1].set_ylabel("Score")

    for n, ax in enumerate([fig.AxesDict[0][2], fig.AxesDict[0][3]]):
        ax.grid(True, alpha=0.3, zorder=0)
        ax.set_xlabel("UMAP-1")
        ax.set_ylabel("UMAP-2")
        ax.set_title(titles[n], fontsize=16, y=1.03)
        
    fig.AxesDict[0][0].grid(True, alpha=0.3, zorder=0)
    fig.AxesDict[0][1].grid(True, alpha=0.3, zorder=0)
    
    return fig, axes

def _plot_score_distribution(ax, scores, n_bins):

    """"""

    lower = scores.mean() - scores.std() / 2
    upper = scores.mean() + scores.std() / 2
    bh = ax.hist(scores, range=(lower, upper), bins=n_bins, color="darkred")
    
    return lower, upper

def _plot_ranked_scores(adata, score_name, ax):

    """"""

    ranked_scores = adata.obs.sort_values(score_name, ascending=False)[
        score_name
    ].values
    ax.scatter(np.arange(len(ranked_scores)), ranked_scores, s=100, c="black")
    ax.scatter(np.arange(len(ranked_scores)), ranked_scores, s=15, c=ranked_scores)


def _plot_pre_post_filtering(adata, geneset_score_key, score_threshold, n_bins=50):

    """"""
    fig, axes = _build_plot(geneset_score_key)

    vmin, vmax = _plot_score_distribution(axes[0], adata.obs[geneset_score_key], n_bins)
    _plot_ranked_scores(adata, geneset_score_key, ax=axes[1])

    umap = adata.obsm["X_umap"]
    adata_f = adata[adata.obs[geneset_score_key] < score_threshold]
    umap_f = adata_f.obsm["X_umap"]

    im0 = axes[2].scatter(
        umap[:, 0], 
        umap[:, 1], 
        c=adata.obs[geneset_score_key], 
        zorder=10, 
        s=2,
        vmin=vmin,
        vmax=vmax,
    )
    im1 = axes[3].scatter(
        umap_f[:, 0], 
        umap_f[:, 1], 
        c=adata_f.obs[geneset_score_key], 
        zorder=10, 
        s=2,
        vmin=vmin,
        vmax=vmax,
    )
    plt.colorbar(im0, cax=axes[4], orientation="vertical")

def _print_summary(adata, adata_f):

    """"""

    n_cells_init = adata.shape[0]
    n_cells_filt = adata_f.shape[0]
    n_cells_removed = adata.shape[0] - adata_f.shape[0]

    print(
        "Cells before filtering:\t{}".format(
            licorice.font_format(str(n_cells_init), ["BOLD"])
        )
    )
    print(
        "Cells after filtering:\t{} ({} cells removed)".format(
            licorice.font_format(str(n_cells_filt), ["BOLD"]),
            licorice.font_format(str(n_cells_removed), ["BOLD"]),
        )
    )

class _GeneSetFiltering:
    def __init__(self, adata=False, gene_set=False, score_name="gene_set_score"):

        """"""

        if gene_set:
            self._gene_set = gene_set

        if adata:
            self._adata = adata

        self._score_name = score_name

    def plot_gene_set(self, adata=False, gene_set=False):

        if gene_set:
            print("gene set is....")
            self._gene_set = gene_set

        sc.pl.umap(self._adata, color=self._gene_set)

    def score_gene_set(self, adata=False, gene_set=False, score_name=False):

        if score_name:
            self._score_name = score_name

        if gene_set:
            self._gene_set = gene_set
            print("ok gene set is...")

        sc.tl.score_genes(self._adata, self._gene_set, score_name=self._score_name)
        sc.pl.umap(self._adata, color=self._score_name)

    def compare_on_score(self, geneset_score_key=False, score_threshold=0.1):

        if geneset_score_key:
            self._score_name = score_name

        self._score_threshold = score_threshold
        #         self._adata_f = self._adata[self._adata.obs[self._score_name] < self._score_threshold]
        _plot_pre_post_filtering(self._adata, self._score_name, self._score_threshold)

    def filter_on_score(self, score_threshold=0.1):

        if score_threshold != self._score_threshold:
            print(
                "Updating score threshold from {} to {}...".format(
                    self._score_threshold, score_threshold
                )
            )
        self._score_threshold = score_threshold
        adata_f =  self._adata[self._adata.obs[self._score_name] < self._score_threshold]
        _print_summary(self._adata, adata_f)

        return adata_f