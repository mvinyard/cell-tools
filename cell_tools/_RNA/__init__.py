
# scRNA-seq __init__.py

__module_name__ = "__init__.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu",])

from ._funcs._calculate_differential_gene_expression import _calculate_differential_gene_expression as DE
from ._funcs._filter_static_genes import _filter_static_genes as filter_static_genes
from ._funcs._find_marker_genes import _find_marker_genes as id_marker_genes

from ._funcs._remove_correlated_genes_adata import _remove_correlated_genes_adata as remove_correlated_genes

from ._funcs._filter_cells_by_gene_set import _GeneSetFiltering as gene_set_filtering