
import glob
import scipy.io as scipy_io
import pandas as pd
from tqdm.notebook import tqdm

def _get_sample_name_between(pathlist, super_dir, minor_dir):

    """
    If organized hierarchically where samples are on the same level, 
    use as: /super_dir/[sample]/minor_dir
    
    """

    SamplePaths = {}

    for path in pathlist:
        SamplePaths[path.split("/tumor")[-1].split("/outs")[0].strip("/")] = path
    return SamplePaths


def _load_data(group_path, super_dir, minor_dir_10x="outs", name=""):

    peak_file_paths = glob.glob(group_path + name +"*.bed")
    print("loading: {}".format(peak_file_paths))
    matr_file_paths = glob.glob(group_path + name + "*.mtx")
    print("loading: {}".format(matr_file_paths))
    

    SamplePaths_peaks = _get_sample_name_between(
        peak_file_paths, super_dir=super_dir, minor_dir=minor_dir_10x
    )
    SamplePaths_mtx = _get_sample_name_between(
        matr_file_paths, super_dir=super_dir, minor_dir=minor_dir_10x
    )
    
    DataDict = {}
    
    for sample, path in tqdm(SamplePaths_peaks.items(), desc='Loading peaks',):
        
        DataDict[sample] = {}
        DataDict[sample]['peaks'] = pd.read_csv(
            path, sep="\t", names=["Chromosome", "Start", "End"],
        )
    
    for sample, path in tqdm(SamplePaths_mtx.items(), desc='Loading data matrices (sparse)',):
        DataDict[sample]['matrix'] = scipy_io.mmread(path)

    
    return DataDict
