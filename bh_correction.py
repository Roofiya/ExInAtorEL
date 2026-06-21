import os
import pandas as pd
import numpy as np
from statsmodels.stats.multitest import multipletests
import argparse

# Set up argument parsing
parser = argparse.ArgumentParser(description="BH Correction of count files.")
parser.add_argument("-o", "--bh_output_dir", required=True, help="Path to the output directory.")
parser.add_argument("-i", "--count_dir", required=True, help="Path to the count files directory.")
parser.add_argument("-g", "--geneset", required=False, help="Gene set name, if applicable.")
args = parser.parse_args()

# Set variables based on the arguments
BH_OUTPUT_DIR = args.bh_output_dir
COUNT_DIR = args.count_dir
GENESET = args.geneset

# Print directory structure to confirm paths (for debugging)
print("Output Directory:", BH_OUTPUT_DIR)
print("Count Directory:", COUNT_DIR)
print("Gene Set:", GENESET)


# Create the output directory if it doesn't exist
if not os.path.exists(BH_OUTPUT_DIR):
    os.makedirs(BH_OUTPUT_DIR)

# Traverse the input directory and its subdirectories using os.walk root: Current directory, dir: Subdirectories, files: individual cvs
for root, dirs, files in os.walk(COUNT_DIR):
    # Process only CSV files with the prefix 'output_'
    csv_files = [f for f in files if f.startswith('na_') and f.endswith('.csv')]
    
    if not csv_files:
        continue

    # Create corresponding subdirectory in the output directory using os.path.relpath
    relative_path = os.path.relpath(root, COUNT_DIR)
    output_subdirectory = os.path.join(BH_OUTPUT_DIR, relative_path)

    if not os.path.exists(output_subdirectory):
        os.makedirs(output_subdirectory)

    # Combine p-values from all files in the current directory
    all_pvalues = []
    file_pvalues_indices = []  # To keep track of the indices for each file's p-values
    file_first_columns = []  # To keep track of the first column for each p-value

    for file in csv_files:
        file_path = os.path.join(root, file)
        df = pd.read_csv(file_path)
        
        
        # Filter rows where both element mutations (M) and non-element mutations are >= 1
        filtered_df = df[(df['Element Mutations(M)'] >= 1) & (df['Non Element Mutations:(m)'] >= 1)]
        pvalues = filtered_df['p-value FET'].values
        all_pvalues.extend(pvalues)
        file_pvalues_indices.append(len(pvalues))
        first_column_values = df.iloc[:, 0].values  # Get the values of the first column
        file_first_columns.extend(first_column_values)  # Track the first column values for each p-value

    # Apply the BH correction
    rejected, pvalues_corrected, _, _ = multipletests(all_pvalues, alpha=0.1, method='fdr_bh')

    # Assign the corrected p-values back to their original files/tests
    corrected_pvalues_index = 0
    for file, num_pvalues in zip(csv_files, file_pvalues_indices):
        file_path = os.path.join(root, file)
        df = pd.read_csv(file_path)
        
        # Filter rows where both element mutations (M) and non-element mutations are >= 1
        filtered_df = df[(df['Element Mutations(M)'] >= 1) & (df['Non Element Mutations:(m)'] >= 1)]
        corrected_pvalues = pvalues_corrected[corrected_pvalues_index:corrected_pvalues_index + num_pvalues]
        
        # Initialize new columns with NaN
        df['Corrected P-Value'] = np.nan
        df['-log10(fdradjp)'] = np.nan

        # Only fill in rows that passed the filter
        mask = (df['Element Mutations(M)'] >= 1) & (df['Non Element Mutations:(m)'] >= 1)
        df.loc[mask, 'Corrected P-Value'] = corrected_pvalues
        df.loc[mask, '-log10(fdradjp)'] = -np.log10(corrected_pvalues)
        corrected_pvalues_index += num_pvalues
        
        # Save the updated DataFrame to the corresponding subdirectory in the output directory
        output_file = os.path.join(output_subdirectory, file)
        df.to_csv(output_file, index=False)

    # Print the results in a tabular format, including the first column from the CSV files
    print(f"{'Folder':<20}{'Original P-Value':<20}{'Significant':<10}{'Corrected P-Value':<20}")
    print('-' * 80)
    for first_col, original_p, is_rejected, corrected_p in zip(file_first_columns, all_pvalues, rejected, pvalues_corrected):
        print(f"{first_col:<20}{original_p:<20}{is_rejected:<10}{corrected_p:<20}")