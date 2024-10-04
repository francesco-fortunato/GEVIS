import pandas as pd

# Load the CSV file into a pandas DataFrame
input_file = 'public/examples/airway_scaledcounts.csv'  # Use forward slashes
output_file = 'filtered_output.csv'  # File to save the filtered data

# Read the CSV file
df = pd.read_csv(input_file)

# Filter out rows where all the columns SRR1039508 to SRR1039521 have zero values
filtered_df = df[~((df.iloc[:, 1:] == 0.0).all(axis=1))]

# Save the filtered data to a new CSV file
filtered_df.to_csv(output_file, index=False)

print(f"Filtered data saved to {output_file}")
