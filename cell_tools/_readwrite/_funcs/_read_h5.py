from scipy import sparse
from anndata import AnnData
import h5py

def _check_hdf5_file_keys(key_list):
    
    """"""
    
    if len(key_list) == 1:
        return key_list[0]

    elif len(key_list) == 0:
        print("No keys found... check file.")

    else:
        print("Too many keys identified.")
        
def access_hdf5_file_keys(file, key=None):
    
    key_list = [key for key in file.keys()]
    
    if key == None:
        data_key = _check_hdf5_file_keys(key_list)
        return file[data_key]
        
    else:
        return file[key]
    
def _read_h5(path):
    
    """
    For the general use-case of reading an h5 file. 
    
    Parameters:
    -----------
    
    Returns:
    --------
    sparse matrix of type '<class 'numpy.float64'> with n stored elements in Compressed Sparse Row format
    
    Notes:
    ------
    (1) Not necessarily compatible with all h5 files. 
    
    (2) Assumes a single key."""
    
    
    h5_file = h5py.File(path, mode="r")
    data = access_hdf5_file_keys(h5_file)
    
    return AnnData(sparse.csr_matrix(data[()]).T)