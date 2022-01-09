import pandas as pd
import numpy as np
import statsmodels.sandbox.stats.multicomp
import scipy.stats
    

def _calculate_p_values(X_norm, cluster_cells, compare_cells, gene_idx):

    """"""

    X1 = X_norm[cluster_cells, :][:, gene_idx].toarray()
    X2 = X_norm[compare_cells, :][:, gene_idx].toarray()
    p_values = np.zeros(gene_idx.sum())

    for i, gene in enumerate(np.nonzero(gene_idx)[0]):
        p_values[i] = scipy.stats.ranksums(X1[:, i], X2[:, i])[1]
    p_values = statsmodels.sandbox.stats.multicomp.multipletests(
        p_values,
        alpha=0.05,
        method="fdr_bh",
    )[1]
    del X1, X2

    return p_values


def _calculate_highest_means(
    X_norm, cluster, clusters, cluster_cells, gene_idx, pseudocount
):

    """Calculate highest vs. 2nd-highest means"""

    m1 = X_norm[cluster_cells, :][:, gene_idx].mean(0).A.squeeze()
    m2 = np.zeros(gene_idx.sum())

    for c in np.unique(clusters):
        if c != cluster:
            i_test = clusters == c
            m_test = norm_counts[i_test, :][:, gene_idx].mean(0).A.squeeze()
            gene_idx_change = m_test - m2 > 0
            m2[gene_idx_change] = m_test[gene_idx_change]

    log_fc = np.log2((m1 + pseudocount) / (m2 + pseudocount))

    return log_fc

def _identify_significant_genes(
    results,
    cluster,
    gene_list,
    gene_idx,
    p_values,
    max_pval,
    log_fc,
    min_log_fold_change,
):

    """"""

    qualifying_genes = (p_values < max_pval) & (log_fc > min_log_fold_change)

    if qualifying_genes.sum() > 0:
        results["group"].extend([cluster for i in range(qualifying_genes.sum())])
        results["gene"].extend(list(gene_list[gene_idx][qualifying_genes]))
        results["p-value"].extend(list(p_values[qualifying_genes]))
        results["log2_fold_change"].extend(list(log_fc[qualifying_genes]))

    if verbose:
        print(
            "cluster={}, n_cells={}, n_diff={}".format(
                cluster, cluster_cells.sum(), qualifying_genes.sum()
            )
        )

    results_df = pd.DataFrame(results)[["group", "gene", "p-value", "log2_fold_change"]]

    return results_df

def _find_marker_genes(X_norm, genes, label_column):

    results = {"group": [], "gene": [], "p-value": [], "log2_fold_change": []}
    min_log_fold_change = np.log2(min_fold_change)
    groups = adata.obs[label_column]
    unique_groups = np.unique(adata.obs[label_column])

    for cluster in unique_groups:

        cluster_cells = groups == cluster
        compare_cells = ~cluster_cells

        gene_idx = (X_norm[cluster_cells, :] > 0).sum(
            0
        ).A.squeeze() / cluster_cells.sum() > min_fraction_expressed

        p_values = _calculate_p_values(X_norm, cluster_cells, compare_cells, gene_idx)
        log_fc = _calculate_highest_means(
            X_norm, cluster, clusters, cluster_cells, gene_idx
        )

        results_df = _identify_significant_genes(
            results,
            cluster,
            gene_list,
            gene_idx,
            p_values,
            max_pval,
            log_fc,
            min_log_fold_change,
        )
        return results_df