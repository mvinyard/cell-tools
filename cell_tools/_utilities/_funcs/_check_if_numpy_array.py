def _check_if_numpy_array(X):
    
    """
    If X is not a numpy array, a numpy array is returned. 
    
    Parameters:
    -----------
    X
        array
        
    Returns:
    --------
    X, potentially transformed to a numpy.ndarray if not already.
        type: numpy.ndarray
    
    Notes:
    ------
    (1) Given the implementation, this function is likely quite inflexible. But works with
        most scenarios where one encounters scipy.sparse matrices that can be converted to
        numpy arrays in this fashion. 
    
    """
    
    if type(X).__module__ != 'numpy':
        return X.toarray()
    else:
        return X