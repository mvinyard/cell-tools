
# utilities __init__.py

__module_name__ = "__init__.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu",])

from ._funcs._strip_10x_barcode_sequence import _strip_10x_barcode_sequence as strip_10x_barcode
from ._funcs._running_quantile import _running_quantile as running_quantile
from ._funcs._sparse_data_functions import _Sparse as sparse