
# _get_score_df.py

__module_name__ = "_get_score_df.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu",])


# import packages #
# --------------- #
import pandas as pd
<<<<<<< HEAD
import nuts_and_bolts
=======


def _make_empty_dict(subdicts=[]):
    
    """
    
    """
    
    Dict = {}
    
    for s in subdicts:
        Dict[s] = {}
        
    return Dict

>>>>>>> 80846f317faa84ef61ddd9b8f4f7568417716a71

def _get_score_dict(adata, group_key):

    names = adata.uns[group_key]["names"]
    scores = adata.uns[group_key]["logfoldchanges"]

<<<<<<< HEAD
    ScoreDict = nuts_and_bolts.create_empty_dict(scores.dtype.names)
=======
    ScoreDict = _make_empty_dict(scores.dtype.names)
>>>>>>> 80846f317faa84ef61ddd9b8f4f7568417716a71
    keys = list(ScoreDict.keys())

    for i in range(len(keys)):
        ScoreDict[keys[i]] = []

        for score in scores:
            ScoreDict[keys[i]].append(score[i])

        ScoreDict[keys[i] + "_feature"] = []
        for j in range(len(names)):
            ScoreDict[keys[i] + "_feature"].append(names[j][i])

    return ScoreDict


def _get_score_df(adata, group_key):

    var_df = adata.var.reset_index(drop=True)
    return pd.DataFrame.from_dict(_get_score_dict(adata, group_key))