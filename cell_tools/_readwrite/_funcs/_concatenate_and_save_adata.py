
import vintools as v
import anndata as a
import pandas as pd
import os


def _concatenate_and_save_adata(
    DataDict,
    merged_adata_outpath=False,
    individual_adata_outpath=False,
    merged_peak_annotations=False,
):

    """Concatenate into a single AnnData object.
    
    Load from merged features. Assumes saved structure.
    """

    v.ut.mkdir_flex(individual_adata_outpath)

    for name, adata in DataDict.items():
        adata.obs = adata.obs.reset_index(drop=True)
        adata.obs.barcoded_sample = adata.obs.barcoded_sample.astype(str)
        adata.obs = adata.obs.set_index("barcoded_sample")
        if individual_adata_outpath:
            outpath = os.path.join(
                individual_adata_outpath, "adata.atac.{}.h5ad".format(name)
            )
            adata.write_h5ad(outpath)

    adata = a.concat(list(DataDict.values()))
    if merged_peak_annotations:
        adata.var = pd.read_csv(merged_peak_annotations, index_col=0)

    if merged_adata_outpath:
        adata.write_h5ad(merged_adata_outpath)

    return adata