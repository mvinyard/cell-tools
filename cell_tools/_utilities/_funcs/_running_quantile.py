
# _running_quantile.py

__module_name__ = "_running_quantile.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu",])


# package imports #
# --------------- #
import numpy as np


def _running_quantile(x, y, fit_percentile, n_bins):

    """
    Calculate the quantile of y in bins of x.
    
    Parameters:
    -----------
    x
    
    y
    
    fit_percentile
    
    n_bins
    
    Returns:
    --------
    x, y
    
    Notes:
    ------
    
    """

    idx = np.argsort(x)
    x = x[idx]
    y = y[idx]

    dx = (x[-1] - x[0]) / n_bins
    x_out = np.linspace(x[0] + dx / 2, x[-1] - dx / 2, n_bins)

    y_out = np.zeros(x_out.shape)

    for i in range(len(x_out)):
        ind = np.nonzero((x >= x_out[i] - dx / 2) & (x < x_out[i] + dx / 2))[0]
        if len(ind) > 0:
            y_out[i] = np.percentile(y[idx], fit_percentile)
        else:
            if i > 0:
                y_out[i] = y_out[i - 1]
            else:
                y_out[i] = np.nan

    return x_out[~np.isnan(y_out)], y_out[~np.isnan(y_out)]