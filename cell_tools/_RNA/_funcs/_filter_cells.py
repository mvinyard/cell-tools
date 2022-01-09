def _filter_cells(dictionary, filter_mask):
    
    """Filter 1- and 2-D entries in a dictionary"""
    
    for key, value in dictionary.items():
        
        if key != "meta":
            if len(value.shape) == 1:
                dictionary[key] = value[filter_mask]
            else:
                dictionary[key] = value[:, filter_mask]
                
    return dictionary