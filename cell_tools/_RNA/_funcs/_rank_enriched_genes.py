
# _rank_enriched_genes.py

__module_name__ = "_rank_enriched_genes.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu",])


# import packages #
# --------------- #
import numpy as np


# local imports #
# ------------- #
from ..._utilities._funcs._sparse_data_functions import _Sparse

def _info_message(cell_mask, gene_idx):

    print("{} cells in group".format(len(cell_mask)))
    print("Considering {} genes".format(len(gene_idx)))
    
def _prepare_adata(adata, column, cell_list):
    
    df_obs = adata.obs.reset_index()
    
    if cell_list:
        return df_obs.loc[df_obs[column].isin(cell_list)].index.astype(int)
        
    else:
        return df_obs.index.astype(int)
            
def _rank_enriched_genes_array(X, 
                               gene_list, 
                               cell_mask, 
                               barcode_column="index", 
                               min_counts=3, 
                               min_cells=3, 
                               verbose=False):
    
    """"""
    
    gene_idx = (X[cell_mask, :] >= min_counts).sum(0).A.squeeze() >= min_cells

    if verbose:
        _info_message(cell_mask, gene_idx)
    
    gene_list = gene_list[gene_idx]
    sparse = _Sparse(X[:, gene_idx])
    
    means = sparse.X.mean(0).A.squeeze()
    stdev = np.sqrt(sparse.variance())
    
    scores = (sparse.X[cell_mask,:].mean(0).A.squeeze() - means) / stdev    
    o = np.argsort(-scores)
    
    return gene_list[o], scores[o]


def _rank_enriched_genes_adata(adata, 
                               gene_list=False, 
                               cell_list=False, 
                               title="gene_enrichment_scores", 
                               return_gene_scores=False, 
                               min_counts=3, 
                               min_cells=3, 
                               cell_barcode_column="index",
                               gene_id_column="gene_name",
                               verbose=True):
    
    """"""
    
    if not genes:
        genes = adata.var['gene_name'].values.tolist()
    
    cell_mask_idx = _prepare_adata(adata, column=cell_barcode_column, cell_list=cells)
        
    genes, scores = _rank_enriched_genes_array(adata.X, 
                                               gene_list=adata.var_names,
                                               cell_mask=cell_mask_idx,
                                               min_counts=min_counts,
                                               min_cells=min_cells,
                                               verbose=verbose,)
    
    
    adata.uns[title] = {}
    adata.uns[title]['genes'] = genes
    adata.uns[title]['scores'] = scores
    if verbose:
        print("\nadata.uns_key added: {}".format(title))
    
    if return_gene_scores:
        return genes, scores