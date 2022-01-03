__module_name__ = "_read_h5.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu",])


# package imports #
# --------------- #
from anndata import AnnData
import h5py
import licorice
import numpy as np
import pandas as pd
from scipy import sparse


def _check_hdf5_file_keys(key_list):
    
    """I guess there should be a main key. 10X uses 'matrix' for at least some h5 outputs. """
    
    if len(key_list) == 1:
        return key_list[0]

    elif len(key_list) == 0:
        print("No keys found... check file.")

    else:
        print("Too many keys identified.")
        
def _access_hdf5_file_keys(file, key='matrix', silent=True):

    key_list = [key for key in file.keys()]

    if key == None:
        key = _check_hdf5_file_keys(key_list)
        
    if not silent:
        print("{}: {}".format(licorice.font_format("h5 key", ['BOLD', 'GREEN']), key_list[0]))
    return file[key]

def _h5_data_object_to_pandas_DataFrame(data_object, level_1_keys, level_2_keys):
    
    """
    
    Assumes str encoding. 
    """
    
    tmp_dict = {}
    
    if level_2_keys:
        data_object_ = data_object[level_1_keys[0]]
        keys = level_2_keys
        
    else:
        data_object_ = data_object
        keys = level_1_keys
    
    for key in keys:
        tmp_dict[key] = np.array(data_object_[key].asstr())
        
    return pd.DataFrame.from_dict(tmp_dict)

def _h5_data_object_to_AnnData(matrix):
    
    """"""
        
    data = matrix['data']
    indices = matrix['indices']
    indptr = matrix['indptr']
    
    return AnnData(sparse.csr_matrix((data, indices, indptr)))

class _h5_to_AnnData:
    
    """"""
    
    def __init__(self, path, silent=True):
        
        """"""
        self.h5_file = h5py.File(h5_path, mode="r")
        self.silent = silent
    
    def load_matrix(self, key="matrix"):
        
        """"""
        
        self.matrix = _access_hdf5_file_keys(h5_file, key, self.silent)
        self.adata = _h5_data_object_to_AnnData(self.matrix)
    
    def load_var(self, level_1_keys=["features"], level_2_keys=['feature_type', 'genome', 'id', 'name']):
        
        """"""
        
        self.var_df = _h5_data_object_to_pandas_DataFrame(self.matrix, level_1_keys, level_2_keys)
        
    def load_obs(self, level_1_keys=["barcodes"], level_2_keys=False):
        
        """"""
        
        self.obs_df = _h5_data_object_to_pandas_DataFrame(self.matrix, level_1_keys, level_2_keys)
        
    def assemble(self, return_adata=True):
        
        """"""
        
        self.adata.obs = self.obs_df
        self.adata.var = self.var_df
        print(self.adata)
        
        return self.adata
    
def _read_h5(path, obs_keys = [["barcodes"], False], var_keys = [['features'], ['feature_type', 'genome', 'id', 'name']], matrix_key="matrix"):
    
    """"""
    
    h5 = _h5_to_AnnData(h5_path)
    h5.load_matrix(matrix_key)
    
    h5.load_obs(level_1_keys=obs_keys[0],
                level_2_keys=obs_keys[1]
               )
    
    h5.load_var(level_1_keys=var_keys[0],
                level_2_keys=var_keys[1]
               )
    
    return h5.assemble()
