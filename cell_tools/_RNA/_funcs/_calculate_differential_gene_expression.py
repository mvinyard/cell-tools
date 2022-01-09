import statsmodels.sandbox.stats.multicomp
import scipy.stats
    
def _calculate_normalized_masked_expression(X, mask):
    return (X[mask, :] > 0).sum(0).A.squeeze() / len(mask)

def _get_gene_mask(X, mask1, mask2, min_fraction_expressed):
    
    """
    
    """
    
    X_1 = _calculate_normalized_masked_expression(X, mask1) > min_fraction_expressed
    X_2 = _calculate_normalized_masked_expression(X, mask2) > min_fraction_expressed
    
    gene_mask = X_1 | X_2
    
    print(X_1.shape, X_2.shape, gene_mask.shape)
    print(X_1.sum(), X_2.sum(), gene_mask.sum())

    return gene_mask

def _calculate_p_values(gene_mask, X1, X2):
    
    """"""
    
    p_values = np.zeros(gene_mask.sum())
    
    for i, gene in enumerate(np.nonzero(gene_mask)[0]):
        p_values[i] = scipy.stats.ranksums(X1[:,i], X2[:,i])[1]
        
    p_values = statsmodels.sandbox.stats.multicomp.multipletests(p_values, alpha=0.05, method='fdr_bh',)[1]
    
    return p_values

def _calculate_differential_gene_expression(adata, mask1, mask2, min_fraction_expressed=0.05, pseudocount=1):
    
    """"""
    
    gene_mask = _get_gene_mask(adata.X, mask1, mask2, min_fraction_expressed)
    
    print(gene_mask.sum())
    X1 = adata.X[mask1,:][:,gene_mask].toarray()
    X2 = adata.X[mask2,:][:,gene_mask].toarray()
    
    m1, m2 = X1.mean(0) + pseudocount, X2.mean(0) + pseudocount
    r = np.log2(m1 / m2)
    
    p_values = _calculate_p_values(gene_mask, X1, X2)

    df = pd.DataFrame({
            'gene': adata.var_names.values.astype(str)[gene_mask],
            'pv': p_values,
            'm1': m1 - pseudocount,
            'm2': m2 - pseudocount,
            'ratio': r
        })

    return df