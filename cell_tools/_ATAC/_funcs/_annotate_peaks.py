
import episcanpy as epi

def _make_episcanpy_compatible_peaks_labels(
    adata,
    Chromosome="Chromosome",
    Start="Start",
    End="End",
    peak_label_added="peak_label",
    delimiter="_",
    unassigned=["intergenic", "unassigned"], 
    key_added='transcriptome', 
    return_adata=True,
    silent=False,
):

    """
    Take individually-formatted ['Chromosome', 'Start', 'End'] adata.var columns and 
    format into: chr_start_end under a single dataframe label, 'peak_label'.
    
    Parameters:
    -----------
    adata
        AnnData object. 
        type: anndata._core.anndata.AnnData
    
    Chromosome
        type: str
    
    Start
        type: str
    
    End
        type: str
    
    peak_label_added
        type: str
    
    delimiter
        type: str
    
    return_adata
        type: bool
    
    Returns:
    --------
    [ optional ] adata
    
    
    Notes:
    ------
    (1) Eventually, it would be nice to build out this function to handle other reformatting of peak labels 
        in an adata.var pandas DataFrame. 
    """

    adata.var = adata.var.reset_index(drop=True)

    adata.var[peak_label_added] = (
        adata.var[Chromosome].astype(str)
        + "_"
        + adata.var[Start].astype(str)
        + "_"
        + adata.var[End].astype(str)
    )
    
    adata.var_names = adata.var[peak_label_added]
    
    if not silent:
        print(adata)

    if return_adata:
        return adata
    
def _make_peak_labels_unique(adata, 
                             unassigned=["intergenic", "unassigned"], 
                             key_added='transcriptome', 
                             peak_label_added='peak_label', 
                             return_adata=True):
    
    """
    EpiScanpy assigns general labels to peaks without specific annotations (i.e., genes or features).
    This function disambiguates those labels by adding an integer to each. 
    
    Parameters:
    -----------
    adata
        AnnData object passed through `epi.tl.find_genes()`
        type: anndata._core.anndata.AnnData
    
    unassigned
        default: ["intergenic", "unassigned"]
        type: list(str)
    
    key_added
        default: 'transcriptome'
        type: str
        
    return_adata
        default: False
        type: bool
        
    Returns:
    --------
    [ optional ] adata
    
    Notes:
    ------
    (1) Should be integrated into a more well thought-out function. 
    """
    
    new_labels = []
    for n, label in enumerate(adata.var.transcript_annotation):
        if label in unassigned:
            new_labels.append(adata.var[peak_label_added][n])
        else:
            new_labels.append(label)
    adata.var[key_added] = new_labels
    adata.var_names = adata.var[key_added]
    
    if return_adata:
        return adata
    
def _annotate_peaks(adata,
                    Chromosome="Chromosome",
                    Start="Start",
                    End="End",
                    peak_label_added="peak_label",
                    delimiter="_",
                    unassigned=["intergenic", "unassigned"], 
                    compatibility_key_added='transcriptome', 
                    episcanpy_key_added="transcript_annotation", 
                    upstream=2000, 
                    feature_type='transcript', 
                    annotation='HAVANA', 
                    raw=False, 
                    gtf_filepath="/home/mvinyard/ref/hg38/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/genes/genes.gtf",
                    return_adata=False,
                    silent=False,
                   ):
    
    """
    Wrapper of `episcanpy.tl.find_genes` with some added handler and reformatting functions.
    
    
    Parameters:
    -----------
    adata
        AnnData object passed through `epi.tl.find_genes()`
        type: anndata._core.anndata.AnnData
    """
    
    adata = _make_episcanpy_compatible_peaks_labels(
        adata,
        Chromosome=Chromosome,
        Start=Start,
        End=End,
        peak_label_added=peak_label_added,
        delimiter=delimiter,
        unassigned=unassigned, 
        key_added=compatibility_key_added, 
        return_adata=True,
        silent=True,
    )
    
    print("Identifying peaks near gene TSS annotations...")
    epi.tl.find_genes(
        adata,
        gtf_file=gtf_filepath,
        key_added=episcanpy_key_added,
        upstream=upstream,
        feature_type=feature_type,
        annotation=annotation,
        raw=raw,
    )
    
    
    adata = _make_peak_labels_unique(adata, 
                                 unassigned=unassigned, 
                                 key_added=compatibility_key_added, 
                                 peak_label_added=peak_label_added,
                                 return_adata=True)
    
    if not silent:
        print(adata)
    if return_adata:
        return adata