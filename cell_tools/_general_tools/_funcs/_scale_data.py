
__module_name__ = "_scale_data.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu"])


# import packages #
# --------------- #
from sklearn.preprocessing import StandardScaler


# local imports #
# ------------- #
from ..._utilities._funcs._check_if_numpy_array import _check_if_numpy_array


def _scale(X, **kwargs):

    """
    Scale input array. Essential step before dimension reduction.

    Parameters:
    -----------
    X
        Unscaled, dense array
        type: numpy.ndarray

    **kwargs
        Optional arguments passed to `sklearn.preprocessing.StandardScaler`

    Returns:
    --------
    scaler
        Data-fit scaler transformer.
        type: sklearn.preprocessing._data.StandardScaler

    X_scaled
        Scaled array
        type: numpy.ndarray

    Notes:
    ------
    (1) Input array must not be a sparse array.
    """

    scaler = StandardScaler(**kwargs)
    X_ = _check_if_numpy_array(X)
    X_scaled = scaler.fit_transform(X_)

    return scaler, X_scaled


def _scale_anndata(adata, return_adata=False, **kwargs):

    """
    Scale input array. Essential step before dimension reduction.

    Parameters:
    -----------
    adata
        The corresponding adata.X is an unscaled, dense array of type=numpy.ndarray
        type: anndata._core.anndata.AnnData

    **kwargs
        Optional arguments passed to `sklearn.preprocessing.StandardScaler`

    return_adata
        Optionally return the adata object. Otherwise, assumed to modified in-place.
        default: False
        type: bool

    Returns:
    --------
    [ optionally ]: adata or None where adata is modified in-place with the
        following updates:
            - adata.layers['X_scaled']
            - adata.uns['scaler']

    Notes:
    ------
    (1) adata.X input array must NOT be sparse.
    """

    adata.uns["scaler"], adata.layers["X_scaled"] = _scale(adata.X, **kwargs)

    if return_adata:
        return adata