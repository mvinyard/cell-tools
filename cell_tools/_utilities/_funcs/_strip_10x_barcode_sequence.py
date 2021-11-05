def _strip_10x_barcode_sequence(df, column=False):

    """
    
    Parameters:
    -----------
    df
        pandas dataframe
        
    column
        key in the pandas dataframe for the 10x barcode sequence
        
    Returns:
    --------
    df
        modified dataframe with new column, df['stripped_bc']
    """

    barcode_sequences = []

    if column:
        bc_col = df[column].str.split("-", expand=True)
    else:
        bc_col = df.index.str.split("-", expand=True)

    for row in bc_col:
        barcode_sequences.append(row[0])

    df["stripped_bc"] = barcode_sequences

    return df