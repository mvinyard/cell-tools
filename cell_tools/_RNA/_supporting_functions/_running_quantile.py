__module_name__ = "__init__.py"
__doc__ = """Main __init__ module."""
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu"])


# -- import packages: --------------------------------------------------------------------
import numpy as np


# -- supporting functions: ---------------------------------------------------------------
def _argsort(x, y):
    idx = np.argsort(x)
    return x[idx], y[idx]


def _dx(x, n_bins):
    return (x[-1] - x[0]) / n_bins


def _xy_binned(x, dx, n_bins):
    x_start = x[0] + (dx / 2)
    x_stop = x[-1] - (dx / 2)

    x_out = np.linspace(x_start, x_stop, n_bins)
    y_out = np.zeros(x_out.shape)
    return x_out, y_out


def _idx(x, x_out, i, dx):
    return np.nonzero((x >= x_out[i] - dx / 2) & (x < x_out[i] + dx / 2))[0]


# -- main module function: ---------------------------------------------------------------
def _running_quantile(x, y, p, n_bins):
    """
    Calculate the quantile of y in bins of x
    """

    x, y = argsort(x, y)
    dx = _dx(x, n_bins)

    x_out, y_out = _xy_binned(x, dx, n_bins)

    for i in range(len(x_out)):
        y_out[i] = np.nan
        idx = _idx(x, x_out, i, dx)
        if len(idx) > 0:
            y_out[i] = np.percentile(y[idx], p)
        elif i > 0:
            y_out[i] = y_out[i - 1]

    return x_out, y_out