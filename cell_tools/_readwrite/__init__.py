
# readwrite __init__.py

__module_name__ = "__init__.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu",])


from ._funcs._make_peak_bed_from_multiome_features import _make_peak_bed_from_multiome_features as make_peak_bed
from ._funcs._read_10x_multiome import _read_10x_multiome as read_10x_multiome
from ._funcs._read_h5 import _read_h5 as read_h5