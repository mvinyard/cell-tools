
# _get_peak_superset.py

__module_name__ = "_get_peak_superset.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu",])


# package imports #
# --------------- #
import numpy as np
import pandas as pd


def _get_peak_superset(DataDict):

    """
    DataDict assumes a structure of DataDict/[sample]/peaks and DataDict/[sample]/matrix
    
    """

    PeakCoordDict = v.ut.EmptyDict(["Chromosome", "Start", "End"])
    for coord_n in PeakCoordDict.keys():
        PeakCoordDict[coord_n] = np.array([])

    for sample, data in DataDict.items():
        for coord_n in PeakCoordDict.keys():
            PeakCoordDict[coord_n] = np.append(
                PeakCoordDict[coord_n], data["peaks"][coord_n]
            )

    return pd.DataFrame.from_dict(PeakCoordDict)