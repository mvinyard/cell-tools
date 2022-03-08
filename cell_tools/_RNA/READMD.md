# RNA


### DE

```python
df = cell.RNA.DE(adata, 
                 mask1, 
                 mask2, 
                 min_fraction_expressed=0.05,
                 pseudocount=1)
```


### filter_static_genes

```python
cell.RNA.filter_static_genes(adata,
                             base_idx=[],
                             min_var_score_percentile=85,
                             min_counts=3,
                             min_cells=3,
                             plot=True,
                             sample_name="Variable genes",
                             return_hv_genes=False,
                            )

```



### id_marker_genes

```python
marker_df = cell.rna.id_marker_genes(X_norm, genes, label_column)
```

### gene_set_filtering

A python `class` rather than a single function to enable interactivity. 
```python
GeneSetFilter = cell.rna.gene_set_filtering(adata=False, 
                            gene_set=False, 
                            score_name="gene_set_score")
GeneSetFiler.plot_gene_set(adata=False, gene_set=False)
GeneSetFiler.score_gene_set(adata=False, gene_set=False, score_name=False)
GeneSetFiler.compare_on_score(geneset_score_key=False, score_threshold=0.1)
adata_f = GeneSetFiler.filter_on_score(score_threshold=0.1)

```