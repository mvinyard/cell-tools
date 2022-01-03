def _check_if_numpy_array(X):
    
    """"""
    
    if type(X).__module__ != 'numpy':
        return X.toarray()
    else:
        return X