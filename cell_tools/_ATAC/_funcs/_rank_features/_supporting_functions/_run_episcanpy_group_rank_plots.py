
import episcanpy as epi

def _run_episcanpy_group_rank_plots(
    adata, group_by, group_key, out_prefix, n_plot_features
):

    epi.pl.rank_feat_groups_tracksplot(
        adata,
        key=group_key,
        save=out_prefix + ".{}.trackplot.pdf".format(group_key),
        n_features=n_plot_features,
    )

    epi.pl.rank_feat_groups_dotplot(
        adata,
        key=group_key,
        save=out_prefix + ".{}.matrix.dotplot.pdf".format(group_key),
        n_features=n_plot_features,
    )

    epi.pl.rank_feat_groups_matrixplot(
        adata,
        key=group_key,
        save=out_prefix + ".{}.matrix.heatmap.pdf".format(group_key),
        n_features=n_plot_features,
    )