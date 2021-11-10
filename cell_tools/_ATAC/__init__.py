
# scATAC-seq __init__.py

__module_name__ = "__init__.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu",])


from ._funcs._merge_reduce_peakset._merge_reduce_peakset_module import _PeakSet as PeakSet
from ._funcs._annotate_peaks import _annotate_peaks as annotate_peaks

from ._funcs._rank_features._RankFeaturesModule import _RankFeatures as RankFeatures