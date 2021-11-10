
# _DataManagerModule.py
__module_name__ = "_DataManagerModule.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu",])

# package imports #
# --------------- #
import pandas as pd
import vintools as v
from tqdm.notebook import tqdm


# local imports #
# ------------- #
from ._read_h5 import _read_h5
from ._file_information import _get_filepaths
from ._concatenate_and_save_adata import _concatenate_and_save_adata

def _annotate_sample_in_adata(DataDict, bc_key="cell_barcode"):

    """"""

    for key in DataDict.keys():
        DataDict[key].obs["sample"] = key
        DataDict[key].obs["barcoded_sample"] = (
            DataDict[key].obs[bc_key] + "_" + DataDict[key].obs["sample"]
        )
        DataDict[key].obs_names = DataDict[key].obs["barcoded_sample"]

    return DataDict

def _print_filename_filepath(filename, filepath):

    """"""

    name_title = v.ut.format_pystring("Name: ", ["BOLD"])
    path_title = v.ut.format_pystring("Path: ", ["BOLD"])
    print("{} {:<10} | {} {}".format(name_title, filename, path_title, filepath))


def _print_pathdict_subdir(dictionary):

    for key, value in dictionary.items():
        for filetype, filepath in value.items():
            _print_filename_filepath(key, filepath)


def _print_pathdict(dictionary):

    for filename, filepath in dictionary.items():
        _print_filename_filepath(filename, filepath)


class _DataManager:
    def __init__(self):

        """"""

        self.DataDict = {}

    def get_h5_path(
        self,
        filepath_dir,
        extension=".h5",
        subdict="h5",
        unique_subdir=0,
        split_on=False,
        silent=True,
    ):

        self.PathDict = _get_filepaths(
            filepath_dir, extension, subdict, unique_subdir, split_on
        )
        if not silent:
            _print_pathdict_subdir(self.PathDict)

    def load_h5(self, extension=".h5"):

        for sample, filepaths_dict in tqdm(self.PathDict.items()):
            self.DataDict[sample] = _read_h5(self.PathDict[sample]["h5"])        
            
    def get_obs_path(
        self,
        filepath_dir,
        extension=".tsv",
        subdict=True,
        unique_subdir=0,
        split_on=False,
        silent=True,
    ):

        """ """

        self._obs_PathDict = _get_filepaths(
            filepath_dir,
            extension,
            subdict=False,
            unique_subdir=unique_subdir,
            split_on=split_on,
        )

        for sample, filepath in self._obs_PathDict.items():
            self.PathDict[sample]["obs"] = filepath
        if not silent:
            _print_pathdict_subdir(self.PathDict)

    def load_obs(self, sep="\t", header=None, names=["cell_barcode"], bc_key="cell_barcode", extension=".h5"):

        """"""

        for sample in self.PathDict.keys():
            path = self.PathDict[sample]["obs"]
            self.DataDict[sample].obs = pd.read_csv(path, sep=sep, header=header, names=names)
            
        self.DataDict = _annotate_sample_in_adata(self.DataDict, bc_key)
        
        
    def concat_save(self, merged_adata_outpath=False, individual_adata_outpath=False, merged_peak_annotations=False):
        
        """"""
        
        self.adata = _concatenate_and_save_adata(self.DataDict, 
                                                 merged_adata_outpath, 
                                                 individual_adata_outpath, 
                                                 merged_peak_annotations,)