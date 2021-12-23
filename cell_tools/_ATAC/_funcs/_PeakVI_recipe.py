
__module_name__ = "_PeakVI_recipe.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu",])


import scvi
import anndata
import scipy
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

scvi.settings.seed = 1218

class _PeakVI:
    
    """"""
    
    def __init__(self, adata_path=False, silent=False):
        
        """"""
        
        self.silent=False
        
        self.min_cell_percentage = 0.05
        
        if adata_path:
            self.adata = anndata.read_h5ad(adata_path)
            if not silent:
                print(self.adata)
            
        
    def preprocess(self):
        
        self.min_cells = int(self.adata.shape[0] * self.min_cell_percentage)
        print("\nFeature filtering threshold of minimum: {} cells".format(self.min_cells))
        sc.pp.filter_genes(self.adata, min_cells=self.min_cells)
        
        if not self.silent:
            print("\n{}".format(self.adata))
        
    def train_model(self, batch_key="sample", save_model="trained.PeakVI.model", overwrite=False):
        
        self.batch_key = batch_key
        scvi.model.PEAKVI.setup_anndata(self.adata, batch_key=batch_key)
        pvi = scvi.model.PEAKVI(self.adata)
        pvi.train()
        pvi.save(save_model, overwrite=overwrite)
        self.save_model_path = "trained.PeakVI.model"
        
    def analysis(self, save_model_path=False, min_dist=0.2, key_added="cluster_pvi", resolution=0.2, use_rep="X_PeakVI"):
        
        if save_model_path:
            self.save_model_path = save_model_path
            
        self.pvi = scvi.model.PEAKVI.load(self.save_model_path, self.adata)
        self.adata.obsm[use_rep] = latent = pvi.get_latent_representation()
        sc.pp.neighbors(self.adata, use_rep=use_rep)
        sc.tl.umap(self.adata, min_dist=min_dist)
        sc.tl.leiden(self.adata, key_added=key_added, resolution=resolution)
        
    def UMAP(self, keys=['cluster_pvi']):
        
        self.color_keys = np.unique(keys.append(self.batch_key))
        sc.pl.umap(self.adata, color=self.color_keys)
        
    def write(self, path="peakVI.adata.h5ad"):
        
        self.adata_path = path
        adata.write(self.adata_path)
        
        if not self.silent:
            print("Saving AnnData to: {}\n{}".format(self.adata_path, self.adata))
        
