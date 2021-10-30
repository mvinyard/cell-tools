
# _merge_reduce_peakset_module.py
__module_name__ = "_merge_reduce_peakset_module.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu",])

# local imports #
# ------------- #

from ._supporting_funcs._load_data import _load_data
from ._supporting_funcs._get_peak_superset import _get_peak_superset
from ._supporting_funcs._unify_peak_matrices import _unify_peak_matrices
from ._supporting_funcs._write_DataDict_mtx_to_h5 import _write_DataDict_mtx_to_h5
from ._supporting_funcs._calculate_master_merged_peaks import _calculate_master_merged_peaks


class _PeakSet:
    
    def __init__(self, out="./"):
        
        """Instantiate"""
        
        
        self.out = out
    
    def load_10x(self, group_path, super_dir, minor_dir_10x="outs"):
        
        """
        Load samples from 10x outs.
                
        Parameters:
        -----------
        group_path
            Path to directory containing all samples
            type: str
    
        super_dir
            Local direoctory containing all samples. 
            type: str
        
        minor_dir_10x
            Directory within each sample containing the filtered_barcode_matrix
            default: "outs"
            type: str
        
        Returns:
        --------
        None
        
        Notes:
        ------
        """
            
        self.DataDict = _load_data(group_path=group_path, super_dir=super_dir, minor_dir_10x=minor_dir_10x)
        
        
    def aggregate_peaks(self):
        
        """
        Aggregate peaks across all samples.
                
        Parameters:
        -----------
        None
        
        Returns:
        --------
        None
        
        Notes:
        ------
        """

        self.aggregate_peak_df = _get_peak_superset(self.DataDict)
        
        
    def define_reduced_merged_peaks(self, out=False):
        
        """
        Defines a new unified peakset across all samples.
        
        Parameters:
        -----------
        out
        
        Returns:
        --------
        Saves new peakset to .csv in `out`.
        
        Notes:
        ------
        """
        if out:
            self.out = out
        
        self.DataDict = _unify_peak_matrices(self.DataDict, self.aggregate_peak_df, self.out)
        
        
    def write(self, out=False):
        
        """
        Write updated data matrices (with merged peaks) to .h5 files, for each sample.
        
        Parameters:
        -----------
        out
        
        Returns:
        --------
        
        Notes:
        ------
        """
        
        if out:
            self.out = out
        
        _write_DataDict_mtx_to_h5(self.DataDict, self.out)