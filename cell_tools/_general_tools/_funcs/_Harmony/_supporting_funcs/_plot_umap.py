
# _plot_umap.py

__module_name__ = "_plot_umap.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu",])


# package imports #
# --------------- #
import matplotlib.pyplot as plt
import numpy as np
import vinplots


def _setup_plot():

    """"""
    
    plot = vinplots.Plot()
    plot.construct(nplots=2, ncols=2, figsize_width=2, figsize_height=1.2)
    plot.style(spines_to_delete=["top", "right"], 
               color="grey", 
               spines_to_color=['bottom', 'left'], 
               spines_positioning_amount=5)
    
    ax = plot.AxesDict[0][0]

    return plot, ax


def _plot_umap(adata, umap_key, plot_by, colors_dict=False):

    """"""
    
    try:
        adata.obs = adata.obs.reset_index()
    except:
        pass
    
    umap = adata.obsm[umap_key]
    if not colors_dict:
        c = vinplots.color_palettes.SHAREseq
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