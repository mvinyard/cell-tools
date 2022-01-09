
import numpy as np

def _calculate_running_quantile(x, y, p, n_bins):
    """
    Calculate the quantile of y in bins of x
    """

    idx = np.argsort(x)
    x, y = x[idx], y[idx]

    dx = (x[-1] - x[0]) / n_bins
    x_out = np.linspace(x[0]+dx/2, x[-1]-dx/2, n_bins)

    y_out = np.zeros(x_out.shape)

    for i in range(len(x_out)):
        idx = np.nonzero((x >= x_out[i]-dx/2) & (x < x_out[i]+dx/2))[0]
        if len(idx) > 0:
            y_out[i] = np.percentile(y[idx], p)
        else:
            if i > 0:
                y_out[i] = y_out[i-1]
            else:
                y_out[i] = np.nan

    return x_out, y_out