import pandas as pd
import numpy as np

# Step 1: Read the CSV file
df = pd.read_csv('helpers\CCA_2024_07_11_raw.CSV', sep=';', index_col=0)

# Step 2: Define quantile normalization function
def quantile_normalize(df_input):
    """
    Perform quantile normalization on a DataFrame.
    """
    # Rank data by sorting each column
    sorted_df = np.sort(df_input.values, axis=0)
    
    # Take the mean of the sorted values across each row (average ranks)
    mean_sorted = np.mean(sorted_df, axis=1)
    
    # Stack the DataFrame and rank the original values
    rank_mean = df_input.stack().groupby(df_input.rank(method='first').stack().astype(int)).mean()
    
    # Get ranks of the original values using the 'min' method
    ranks = df_input.rank(method='min', axis=0).stack().astype(int)
    
    # Map ranks to the averaged sorted values (mean of ranks)
    df_normalized = ranks.map(rank_mean).unstack()
    
    return df_normalized

# Step 3: Perform quantile normalization
df_normalized = quantile_normalize(df)

# Step 4: Save the normalized data to CSV
df_normalized.to_csv('quantile_normalized.csv', sep=';')

print("Quantile normalization complete and saved to 'quantile_normalized.csv'.")
