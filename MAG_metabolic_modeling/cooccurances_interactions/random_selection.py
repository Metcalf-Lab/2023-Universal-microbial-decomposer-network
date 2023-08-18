import pandas as pd
import random
import sys

if len(sys.argv) != 3:
    print("Usage: python random_selection.py input_file.tsv output_file.tsv")
    sys.exit(1)

input_file_path = sys.argv[1]
output_file_path = sys.argv[2]

# Read the original TSV file
df = pd.read_csv(input_file_path, sep='\t', header=None)

# Define the number of loops and samples
num_loops = 100
num_samples = 20

# Initialize an empty list to store results
results = []

# Perform the random selection process for each loop
for loop in range(1, num_loops + 1):
    selected_indices = random.sample(range(len(df)), num_samples)
    selected_objects = df.iloc[selected_indices][0].tolist()
    
    for obj in selected_objects:
        results.append({'Loop': f'Loop{loop}', 'SelectedObjects': obj})

# Create a DataFrame from the results list
results_df = pd.DataFrame(results)

# Write results to the output TSV file
results_df.to_csv(output_file_path, sep='\t', index=False, header=False)

