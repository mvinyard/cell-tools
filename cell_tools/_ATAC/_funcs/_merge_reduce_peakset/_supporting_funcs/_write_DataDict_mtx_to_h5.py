import h5py
import os

def _write_matrix_h5(array, h5_filepath, name):

    """
    
    """

    with h5py.File(h5_filepath, "w") as hf:
        hf.create_dataset(name, data=array)


def _peak_df_to_array(df):

    """"""

    return df.drop(["Cluster", "Start", "End", "Chromosome"], axis=1).values


def _write_DataDict_mtx_to_h5(DataDict, out):
    
    """"""
    
    for sample_key in DataDict.keys():

        mtx_array = _peak_df_to_array(DataDict[sample_key]['merged_reduced_peak_mtx'])
        h5_filepath = os.path.join(out, "merged_reduced_peak_mtx.{}.h5".format(sample_key))
        _write_matrix_h5(mtx_array, h5_filepath, "merged.peaks.{}".format(sample_key))