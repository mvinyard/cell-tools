
import pandas as pd
import numpy as np
from ._calculate_master_merged_peaks import _calculate_master_merged_peaks

def _define_merged_peakset_each_dataset(DataDict, ClusteredRanges_df):

    for key, values in DataDict.items():
        original_single_sample = DataDict[key]["peaks"]
        DataDict[key]["merged_peaks"] = pd.merge(
            ClusteredRanges_df,
            original_single_sample,
            on=["Chromosome", "Start", "End"],
            how="inner",
        )
        
    return DataDict

def _unify_peak_matrices(DataDict, aggregate_peak_df, out):

    """
    
    Parameters:
    -----------
    DataDict
    
    merged_reduced_peak_df
    
    Return:
    -------
    
    Notes:
    ------
    """
    
    
    ClusteredRanges_df, merged_reduced_peak_df = _calculate_master_merged_peaks(aggregate_peak_df, out)
    DataDict = _define_merged_peakset_each_dataset(DataDict, ClusteredRanges_df)
    
    for a, key in enumerate(DataDict.keys()):

        DataDict[key]['merged_reduced_peak_mtx'] = merged_reduced_peak_df.merge(
            pd.DataFrame(
                DataDict[key]["matrix"].toarray(),
                index=DataDict[key]["merged_peaks"].Cluster,
            )
            .reset_index()
            .groupby("Cluster")
            .aggregate(np.sum),
            on="Cluster",
            how="left",
        ).replace(np.nan, 0)
        
    return DataDict