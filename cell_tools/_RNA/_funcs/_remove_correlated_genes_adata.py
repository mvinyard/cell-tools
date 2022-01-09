

# import packages #
# --------------- #
import numpy as np
import scipy.optimize


# local imports #
# ------------- #
from ..._utilities._funcs._sparse_data_functions import _Sparse


def _remove_correlated_genes(X, gene_list, gene_exclusion_set, test_gene_idx, min_corr = 0.1):
    """
    Remove signature-correlated genes from a list of test genes 
    
    Arguments:
    ----------
    X: scipy.sparse.csc_matrix, shape (n_cells, n_genes)
        - full counts matrix
    gene_list: numpy array, shape (n_genes,)
        - full gene list
    exclude_corr_genes_list: list of list(s)
        - Each sublist is used to build a signature. Test genes correlated
          with this signature will be removed
    test_gene_idx: 1-D numpy array
        - indices of genes to test for correlation with the 
          gene signatures from exclude_corr_genes_list
    min_corr: float (default=0.1)
        - Test genes with a Pearson correlation of min_corr or higher 
          with any of the gene sets from exclude_corr_genes_list will
          be excluded

    Returns:
    --------
        numpy array of gene indices (subset of test_gene_idx) that 
        are not correlated with any of the gene signatures
    """
    
    seed_idx_list = []
    exclusion_idx = []
    
    if type(gene_exclusion_set) != list:
        gene_exclusion_set = [gene_exclusion_set]
    
    for sublist in gene_exclusion_set:
        seed_idx_list.append(np.array([i for i in range(len(gene_list)) if gene_list[i] in sublist], dtype=int))
    
    for index_set in range(len(seed_idx_list)):
        seed_idx = seed_idx_list[index_set][X[:,seed_idx_list[index_set]].sum(axis=0).A.squeeze() > 0]
        if type(seed_idx) is int:
            seed_idx = np.array([seed_idx], dtype=int)
        elif type(seed_idx[0]) is not int:
            seed_idx = seed_idx[0]
        
        sparse = _Sparse(X[:, seed_idx])
        tmp = sparse.z_score()
        tmp = tmp.sum(1).A.squeeze()

        correlation = np.zeros(len(test_gene_idx))
        for gene_i in range(len(correlation)):
            correlation[gene_i],_ = scipy.stats.pearsonr(tmp, X[:,test_gene_idx[gene_i]].A.squeeze())
        del tmp
        
        exclusion_idx.extend([test_gene_idx[i] for i in range(len(test_gene_idx)) if (correlation[i]) >= min_corr])
    
    exclusion_idx = np.array(exclusion_idx)

    return np.array([g for g in test_gene_idx if g not in exclusion_idx], dtype=int)

def _remove_correlated_genes_adata(adata, signature_genes, query_genes="highly_variable", verbose=True):

    """
    Remove genes correlated with the expression signature of a providued signature gene list.  
    
    adata

    signature_genes
    
    query_genes
        default: "highly_variable"
    """

    df_var = adata.var.reset_index()
    
    if query_genes == "highly_variable":
        query_genes = adata.uns["highly_variable_genes_idx"]

    elif query_genes == False:
        query_genes = range(len(df_var))

    n_genes = len(query_genes)
    
    all_genes = df_var.loc[df_var["index"].isin(adata.var_names)].index.astype(int)

    gene_exclusion_set = df_var.loc[
        df_var["gene_id"].isin(signature_genes)
    ].index.astype(int)

    filtered_genes_idx = _remove_correlated_genes(
        adata.X, all_genes, gene_exclusion_set, query_genes
    )
    
    if verbose:
        print("\nRemoving genes correlated with the expression of those provided...\n")
    adata_ = adata[:, filtered_genes_idx]
    adata_.obs.index = adata_.obs.index.astype(str)
    adata_.var.index = adata_.var.index.astype(str)
    adata_ = adata_.copy()
    
    if verbose:
        n_removed = int(n_genes - adata_.shape[1])
        print("\nFrom {} genes, {} genes were removed leaving {} genes.".format(n_genes, n_removed, adata_.shape[1]))
        print("\n{}".format(adata_))
    
    return adata_