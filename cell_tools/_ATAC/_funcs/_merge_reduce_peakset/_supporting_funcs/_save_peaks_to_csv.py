import os
import vintools as v

def _save_peaks_to_csv(peak_df, out, filename="MergedPeaks.csv"):
    
    v.ut.mkdir_flex(out)
    filepath = os.path.join(out, filename)
    
    peak_df.to_csv(filepath)
    
    return filepath
