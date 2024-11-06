import pandas as pd

# Read the CSV file
file_path1 = '..\miRNA_old\CCA\RAW_2024_07_11_CLEAN\CCA_2024_07_11_raw.CSV'  # Replace with your actual file path
file_path2 = '..\miRNA_old\CCA\RAW_2024_07_11_CLEAN\CCA_2024_07_11_casistica.CSV'  # Replace with your actual file path
df1 = pd.read_csv(file_path1, sep=';')
df2 = pd.read_csv(file_path2, sep=';')

# Remove the 9th column (index 8)
df1.drop(df1.columns[10], axis=1, inplace=True)
# Remove the 9th column (index 8)
df1.drop(df1.columns[8], axis=1, inplace=True)
# Remove the 9th column (index 8)
df1.drop(df1.columns[6], axis=1, inplace=True)

# Save the modified dataframe to a new CSV file
output_file_path1 = 'raw.csv'  # Replace with the desired output file name
df1.to_csv(output_file_path1, sep=';', index=False)

# Remove the 9th column (index 8)
df2.drop(df2.columns[10], axis=1, inplace=True)
# Remove the 9th column (index 8)
df2.drop(df2.columns[8], axis=1, inplace=True)
# Remove the 9th column (index 8)
df2.drop(df2.columns[6], axis=1, inplace=True)

# Save the modified dataframe to a new CSV file
output_file_path2 = 'metadata.csv'  # Replace with the desired output file name
df2.to_csv(output_file_path2, sep=';', index=False)

print("The 9th column has been removed and the file has been saved.")
