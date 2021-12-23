from tqdm import tqdm_notebook as tqdm
import vintools as v
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


def _fetch_group_idx(df_obs, group_key, group):
    return df_obs.loc[df_obs[group_key] == group].index.astype(int)


def _make_group_rank_df(RankDict, group):

    group_dict = {}

    for episcanpy_key in ["names", "scores", "pvals", "pvals_adj", "logfoldchanges"]:
        group_dict[episcanpy_key] = RankDict[episcanpy_key][group]

    return pd.DataFrame.from_dict(group_dict)


# def _write_multi_df_to_excel(
#     df_list, workbook_path="./new_workbook.xlsx", sheetnames=False, index=False, silent=False
# ):


def _calculate_FeatureSumDict(
    adata,
    uns_group_key,
    group_key,
    title="tumor.scATACseq",
    n_features=5000,
    transcriptome_col="transcriptome",
):

    """"""

    FeatureSumDict = {}

    RankDict = _format_EpiScanpy_RankDict(adata, uns_group_key)

    df_var = adata.var.reset_index()
    df_obs = adata.obs.reset_index()

    X = adata.X

    for group in tqdm(RankDict["names"].columns):
        FeatureSumDict[group] = {}
        feature_idx = df_var[
            df_var[transcriptome_col].isin(RankDict["names"][group].values)
        ].index.astype(int)
        group_idx = _fetch_group_idx(df_obs, group_key, group)
        feature_vector = X[:, feature_idx]
        total_sum, group_sum = (
            feature_vector.sum(axis=0),
            feature_vector[group_idx].sum(axis=0),
        )

        FeatureSumDict[group] = pd.DataFrame(
            np.vstack([group_sum, total_sum, (group_sum / total_sum)])
        ).T
        FeatureSumDict[group].columns = [
            "{}_sum".format(group),
            "total_sum",
            "group_proportion",
        ]

        group_rank_df = _make_group_rank_df(RankDict, group)
        expanded_names = group_rank_df["names"].str.split(";", expand=True)
        FeatureSumDict[group] = pd.concat(
            [group_rank_df, FeatureSumDict[group], expanded_names], axis=1
        )

    v.ut.df_to_excel(
        list(FeatureSumDict.values()),
        workbook_path="{}.top{}.logfoldchange.differential_features.{}.xlsx".format(
            title, n_features, uns_group_key
        ),
        sheetnames=list(FeatureSumDict.keys()),
    )

    return FeatureSumDict