import vintools as v
import pandas as pd
import numpy as np

def _prepare_diffCA_Dict_for_Excel(adata, score_df, groupby):

    diffCA_Dict = {}
    for sample in np.array(adata.obs[groupby].unique()):
        diffCA_Dict["{}_vs_rest".format(sample)] = score_df.filter(regex=sample)

    merged_df = pd.concat(list(diffCA_Dict.values()), axis=1)
    diffCA_Dict["all_{}_vs_rest".format(groupby)] = merged_df

    return diffCA_Dict


def _write_differential_features_to_excel(
    adata, score_df, out_prefix, n_features, groupby
):

    diffCA_Dict = _prepare_diffCA_Dict_for_Excel(adata, score_df, groupby)

    v.ut.df_to_excel(
        list(diffCA_Dict.values()),
        workbook_path="{}.top{}.logfoldchange.differential_features.{}.xlsx".format(
            out_prefix, n_features, groupby
        ),
        sheetnames=list(diffCA_Dict.keys()),
    )