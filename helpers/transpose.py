import pandas as pd

# Load the CSV file
file_path = 'metadata.csv'  # Replace with the actual file path
df = pd.read_csv(file_path)

# Transpose the DataFrame
df_transposed = df.T

# Optionally, reset the index if needed
df_transposed.columns = df_transposed.iloc[0]  # Set the first row as header
df_transposed = df_transposed.drop(df_transposed.index[0])  # Drop the old header row

# Save the transposed DataFrame to a new CSV file
df_transposed.to_csv('transposed_file.csv', index=False)

# Output the transposed DataFrame
print(df_transposed)
