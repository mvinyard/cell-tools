
import os
import scanpy
import scipy.sparse

def _pbmc3k_preprocessed(save_to_disk=False, verbose=True):
    
    """
    Download and return a Scanpy preprocessed PBMC 3k cell dataset from 10x Genomics. 
    
    Parameters:
    -----------
    save_to_disk
        type: bool
        default: False
        
    Returns:
    --------
    adata
        type: anndata
    """
    
    adata = scanpy.datasets.pbmc3k_processed()
    
    if not save_to_disk:
        path = "./data/pbmc3k_processed.h5ad"
        os.system("rm -r {}".format(os.path.dirname(path)))
    
    if adata.X.__class__.__module__ != 'scipy.sparse.csr':
        if verbose:
            print("Converting to scipy.sparse.csr_matrix...\n")
        adata.X = scipy.sparse.csr_matrix(adata.X)
    
    if verbose:
        print(adata)
    
    return adata