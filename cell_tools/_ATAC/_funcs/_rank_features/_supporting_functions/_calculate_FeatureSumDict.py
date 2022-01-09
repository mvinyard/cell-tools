from tqdm import tqdm_notebook as tqdm
import pandas as pd
import numpy as np

def _format_EpiScanpy_RankDict(adata, uns_group_key):

    """Reformat the default output of EpiScanpy's rankfeatures function."""

    RankDict = {}

    for episcanpy_key in ["names", "scores", "pvals", "pvals_adj", "logfoldchanges"]:
        episcanpy_key_dict = {}
        for n, ranked_row in enumerate(adata.uns[uns_group_key][episcanpy_key]):
            episcanpy_key_dict[n] = ranked_row.tolist()
        df = pd.DataFrame.from_dict(episcanpy_key_dict, orient="index")
        df.columns = adata.uns[uns_group_key][episcanpy_key].dtype.names
        RankDict[episcanpy_key] = df

    RankDict["params"] = adata.uns[uns_group_key]["params"]

    return RankDict


def _fetch_feature_idx(df_var, feature, transcriptome_col="transcriptome"):
    return df_var[df_var[transcriptome_col] == feature].index.astype(int)


def _fetch_group_idx(df_obs, group_key, group):
    return df_obs.loc[df_obs[group_key] == group].index.astype(int)

def _calculate_feature_group_representation(
    X, df_obs, df_var, feature, group_key, group, transcriptome_col="transcriptome"
):

    """Loop through features for each sample in the RankDict
    
    """

    feature_idx = _fetch_feature_idx(df_var, feature, transcriptome_col)
    group_idx = _fetch_group_idx(df_obs, group_key, group)
    feature_vector = X[:, feature_idx]

    return feature_vector.sum(), feature_vector[group_idx].sum()

def _calculate_FeatureSumDict(
    adata, uns_group_key, group_key, transcriptome_col="transcriptome"
):

    """"""

    FeatureSumDict = {}

    RankDict = _format_EpiScanpy_RankDict(adata, uns_group_key)

    df_var = adata.var.reset_index()
    df_obs = adata.obs.reset_index()

    X = adata.X

    for group in tqdm(RankDict["names"].columns):
        FeatureSumDict[group] = {}
        FeatureSumDict[group]["feature_sum"] = []
        FeatureSumDict[group]["feature_group_sum"] = []

        group_col_ = RankDict["names"][group].values
        for feature in tqdm(group_col_):
            feature_sum, feature_group_sum = _calculate_feature_group_representation(
                X, df_obs, df_var, feature, group_key, group, transcriptome_col
            )
            FeatureSumDict[group]["feature_sum"].append(feature_sum)
            FeatureSumDict[group]["feature_group_sum"].append(feature_group_sum)
        group_sum = np.array(FeatureSumDict[group]["feature_group_sum"])
        total_sum = np.array(FeatureSumDict[group]["feature_sum"])
        FeatureSumDict[group] = pd.DataFrame(
            np.vstack([group_sum, total_sum, group_sum / total_sum])
        ).T

    return FeatureSumDict, RankDict