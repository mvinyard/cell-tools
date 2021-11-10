import numpy as np
import vintools as v
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec


def _setup_plot():

    """"""

    plot = plt.Figure(figsize=np.array(v.pl.mpl_dimensions()) * [2, 1.2])
    gs = GridSpec(1, 2)
    ax = plot.add_subplot(gs[0, 0])
    spines = v.pl.ax_spines(ax)
    spines.delete(["top", "right"])
    spines.set_color("grey")
    spines.set_position(position_type="outward", amount=20)

    return plot, ax


def _plot_umap(adata, umap_key, plot_by, colors_dict=False):

    """"""
    
    try:
        adata.obs = adata.obs.reset_index()
    except:
        pass
    
    umap = adata.obsm[umap_key]
    if not colors_dict:
        c = v.pl.share_seq()["colors"]
    plot, ax = _setup_plot()

    for n, i in enumerate(adata.obs[plot_by].unique()):

        if colors_dict:
            c_ = colors_dict[i]
        else:
            c_ = c[n]

        idx = adata.obs.loc[adata.obs[plot_by] == i].index.astype(int)
        ax.scatter(umap[:, 0][idx], umap[:, 1][idx], c=c_, label=i, s=5, alpha=0.8)

    ax.set_title("Harmonized Data")
    ax.legend(bbox_to_anchor=(1.05, 1.05), edgecolor="white", markerscale=2)
    plt.tight_layout()
    
    return plot