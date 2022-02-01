
# _load_data.py

__module_name__ = "_load_data.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu",])


# package imports #
# --------------- #
import os
import pyrequisites as pyrex


def _save_peaks_to_csv(peak_df, out, filename="MergedPeaks.csv"):
    
    pyrex.mkdir_flex(out)
    filepath = os.path.join(out, filename)
    
    peak_df.to_csv(filepath)
    
    return filepath
