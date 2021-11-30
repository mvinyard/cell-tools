
# _calculate_master_merged_peaks.py

__module_name__ = "_calculate_master_merged_peaks.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu",])


# package imports #
# --------------- #
import episcanpy as epi
import numpy as np
import pandas as pd
import pyranges as pr


# local imports #
# ------------- #
from ._save_peaks_to_csv import _save_peaks_to_csv


def _pyranges_bedfile_cluster(chromosomes, starts, ends):

    ClusteredRanges_df = (
        pr.PyRanges(chromosomes=chromosomes, starts=starts, ends=ends)
        .cluster()
        .as_df()
        .drop_duplicates()
    )

    return ClusteredRanges_df


def _group_overlapping_peaks(peak_df):

    """"""

    ClusteredRanges_df = _pyranges_bedfile_cluster(
        peak_df.Chromosome, peak_df.Start, peak_df.End
    )
    grouped = ClusteredRanges_df.groupby("Cluster")

    return ClusteredRanges_df, grouped


def _calculate_master_merged_peaks(peak_superset_df, out):

    """"""

    ClusteredRanges_df, GroupedRanges_df = _group_overlapping_peaks(peak_superset_df)

    chrom_vals = GroupedRanges_df["Chromosome"]
    start_vals = GroupedRanges_df["Start"].aggregate(np.min)
    end_vals = GroupedRanges_df["End"].aggregate(np.max)

    master_merged_peaks_df = (
        pd.DataFrame([start_vals, end_vals])
        .T.reset_index()
        .merge(ClusteredRanges_df, on=["Cluster", "Start", "End"])
    )
    
    
    _save_peaks_to_csv(master_merged_peaks_df, out, filename="MergedPeaks.csv")
    
    return ClusteredRanges_df, master_merged_peaks_df