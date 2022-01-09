import licorice
import numpy as np


def _get_var_idx(adata, modality, modality_key):

    return np.where(adata.var[modality_key] == modality)[0]


def _get_obs_idx(adata, var_idx):

    return np.where(adata[:, var_idx].X.sum(1) > 0)[0]


def _subset_modality(
    adata, modality, modality_key="modality", copy=False, silent=False
):

    var_idx = _get_var_idx(adata, modality, modality_key)
    obs_idx = _get_obs_idx(adata, var_idx)
    adata_modal = adata[obs_idx, var_idx]

    if not silent:
        modality_title = licorice.font_format(modality + ":", ["BOLD", "CYAN"])
        print("{}\n{}\n".format(modality_title, adata_modal))

    if copy:
        return adata_modal.copy()
    else:
        return adata_modal
    
    
def _calculate_min_cell_threshold(
    adata,
    modality_key,
    modality_column_key,
    copy,
    silent,
    MinCellsDict,
    modality_fraction,
):

    """"""

    adata_subset = _subset_modality(
        adata, modality_key, modality_column_key, copy, silent
    )
    modality_cell_count = (adata_subset.X != 0).sum(0).flatten()
    MinCellsDict[modality_key] = int(adata_subset.shape[0] * modality_fraction)

    return MinCellsDict


def _filter_features_multiomic(
    adata,
    modality_key="modality",
    min_cell_threshold=0.01,
    silent=False,
):

    """"""
    
    n_features_unfiltered = adata.var.shape[1]

    modalities = adata.var[modality_key].unique()
    MinCells = {}
    filtered_features = np.array([])

    for modality in modalities:
        MinCells = _calculate_min_cell_threshold(adata, 
                                                 modality, 
                                                 modality_key,
                                                 copy=False, 
                                                 silent=silent, 
                                                 MinCellsDict=MinCells,
                                                 modality_fraction=min_cell_threshold
                                                )

    for key, threshold in MinCells.items():
        idx = np.where(adata[:, adata.var[modality_key] == key].X.sum(1) > threshold)[0]
        filtered_features = np.append(filtered_features, idx)    

    n_features_filtered = n_features_unfiltered - len(filtered_features)
    
    print("Total filtered features: {}".format(n_features_filtered))
    
    return adata[:, filtered_features]