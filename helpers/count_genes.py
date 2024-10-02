import pandas as pd
import numpy as np

# Load your CSV file
csv_file = 'public\examples\metadata_CCA.csv'
df = pd.read_csv(csv_file, index_col=0)

# Transpose the DataFrame to make rows into columns
df_transposed = df.T

# Check if 'Survival time' exists as a column now
if 'Survival time' in df_transposed.columns:
    # Apply the conversion to days (assuming the time is in months)
    df_transposed['Survival time'] = df_transposed['Survival time'].apply(lambda x: np.round(float(x.split()[0]) * 30.4375) if isinstance(x, str) else x)

# Transpose back if needed or save the modified DataFrame
df_modified = df_transposed.T
df_modified.to_csv('modified_file.csv')
