
# _isolate_filename.py

__module_name__ = "_isolate_filename.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu",])


# package imports #
# --------------- #
import os
import re
import glob


def _isolate_filename(filepath, unique_subdir, extension, split_on):

    """
    Get the name (substring) of a filename, given a path. Most likely to be used for metadata purposes
    or organization of data in a dictionary, etc. 
    
    Parameters:
    -----------
    filepath
        path to the file of interest. 
    
    unique_subdir
        number of subdirectories that should be navigated upward. 
        For example, if path is the following:
            /home/mvinyard/data/files/files.txt
            
        0: files.txt
        1: files
        2: data:
        3: mvinyard
        4: home
        
        type: bool (False) or int
    
    extension
        file extension
        type: str
        
    split_on
        If not False, provide a tuple of a string on which to be split and an integer indicating which 
        item within the tuple should be isolated. 
        default: False
        type: bool (False) or list(str, int)
        
    
    Returns:
    --------
    filename
        name of the file. used for metadata purposes and organization. 
        type: str
        
    
    Notes:
    ------
    (1) 
    
    """

    if unique_subdir == 0:
        filename = os.path.basename(filepath)

    else:
        filename = os.path.dirname(filepath)
        for i in range(unique_subdir):
            filename = os.path.dirname(filename)
        filename = os.path.basename(filename)
    if split_on:
        filename = filename.split(split_on[0])[split_on[1]]

    filename = re.sub(extension, "", filename)

    return filename


def _get_filepaths(filepath_dir, extension, subdict=False, unique_subdir=0, split_on=False):

    """
    filepath_dir
        directory containing filepaths to be globbed.
        
    extension
        filepath extension
        e.g., .h5, .txt
        
    unique_subdir
        number of subdirectories that should be navigated upward. 
        For example, if path is the following:
            /home/mvinyard/data/files/files.txt
            
        0: files.txt
        1: files
        2: data:
        3: mvinyard
        4: home
        
        type: bool (False) or int
    
    split_on
        split_on is a two-component tuple consisting of 
        (1) the substring by which to split and (2) the 
        resulting component taken. e.g., [".", 2]
        default: False
        type: bool
        
    Returns:
    --------
    PathDict
        Dictionary with keys representative of filenames and values are corresponding filepaths. 
    """

    PathDict = {}

    filepaths = glob.glob(filepath_dir + "*" + extension)

    for path in filepaths:

        filename = _isolate_filename(
            filepath=path,
            unique_subdir=unique_subdir,
            extension=extension,
            split_on=split_on,
        )
        if subdict:
            PathDict[filename] = {}
            PathDict[filename][subdict] = path
        else:
            PathDict[filename] = path
            
    return PathDict