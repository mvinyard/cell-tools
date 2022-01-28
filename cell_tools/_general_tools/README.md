# General Tools

### PCA

```python
cell.tl.pca(adata, n_components=2, return_adata=False, **kwargs)
```

### UMAP

```python
cell.tl.umap(adata, use="X_pca", n_components=2, return_adata=False, **kwargs)
```

### Harmony

```python
Harmony = cell.tl.Harmonize(title="sc-experiment")
Harmony.run(adata, batch_key, use_key="X_pca")
Harmony.umap(n_components=2):
Harmony.plot(plot_by = [...,], colors_dict=False)
Harmony.save():
```