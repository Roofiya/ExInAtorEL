import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import regex as re
import sys
import os
import time
import csv
import multiprocessing as mp
import random
from math import *
import os.path
import subprocess
import scipy.stats as stats
import argparse

pd.set_option('display.max_columns', None)
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import matplotlib.colors as mcolors


# Set up argument parsing
parser = argparse.ArgumentParser(description="Process output count files.")
parser.add_argument("-i", "--base_output_dir", required=True, help="Path to the output directory.")
parser.add_argument("-o", "--count_dir", required=True, help="Path to the count files directory.")
parser.add_argument("-g", "--geneset", required=False, help="Gene set name, if applicable.")
args = parser.parse_args()

# Set variables based on the arguments
BASE_OUTPUT_DIR = args.base_output_dir
COUNT_DIR = args.count_dir
GENESET = args.geneset

# Print directory structure to confirm paths (for debugging)
print("Output Directory:", BASE_OUTPUT_DIR)
print("Count Directory:", COUNT_DIR)
print("Gene Set:", GENESET)

output_store = os.path.join(COUNT_DIR,f"output_{os.path.basename(BASE_OUTPUT_DIR)}.csv")
output_store_na = os.path.join(COUNT_DIR,f"na_output_{os.path.basename(BASE_OUTPUT_DIR)}.csv")

dirs = [d for d in os.listdir(BASE_OUTPUT_DIR) if os.path.isdir(os.path.join(BASE_OUTPUT_DIR, d))]
print(dirs)

# Rename the folders

path = BASE_OUTPUT_DIR
for filename in os.listdir(path):
    if filename.startswith('19_'):
        new_filename = filename[3:]
        os.rename(os.path.join(path, filename), os.path.join(path, new_filename))

for filename in os.listdir(path):
    if filename.endswith('.bed'):
        new_filename = filename[:len(filename)-4]
        os.rename(os.path.join(path, filename), os.path.join(path, new_filename))
        print(new_filename)

# Get the current working directory
current_directory = os.getcwd()

# Print the current directory
print(current_directory)

# Define the directory where the TSV files are stored
directory = BASE_OUTPUT_DIR
print(directory)


##########################################   kmer_counts_data.tsv     ##############################################

# Create an empty DataFrame to store the combined data
combined_data = pd.DataFrame()

# Loop through each folder in the directory
for folder in os.listdir(directory):
    folder_path = os.path.join(directory, folder)
    #print(folder_path)
    if os.path.isdir(folder_path):
        # Loop through each TSV file in the folder
        for file in os.listdir(folder_path):
            file_path = os.path.join(folder_path, file)
            if file.endswith('kmer_counts_output.tsv'):
                # Read the TSV file into a DataFrame and append it to the combined data
                data = pd.read_csv(file_path, sep='\t')
                combined_data = pd.concat([combined_data, data], axis=0, ignore_index=True)

# Write the combined data to a new TSV file
output_path = os.path.join(COUNT_DIR,'kmer_counts_data.tsv')
combined_data.to_csv(output_path, sep='\t', index=False)

#######################################  counts_data.tsv   ######################################################
combined_data1 = pd.DataFrame()
# Loop through each folder in the directory
for folder in os.listdir(directory):
    folder_path = os.path.join(directory, folder)
    if os.path.isdir(folder_path):
        # Loop through each TSV file in the folder
        for file in os.listdir(folder_path):
            file_path = os.path.join(folder_path, file)
            if file == 'counts_output.tsv':
                # Read the TSV file into a DataFrame and append it to the combined data.
                # Main.py writes an 'Element' column as the first column; only add
                # it from the folder name if an older file without it is encountered.
                data = pd.read_csv(file_path, sep='\t')
                if 'Element' not in data.columns:
                    data.insert(0, 'Element', folder)
                combined_data1 = pd.concat([combined_data1, data], axis=0, ignore_index=True)

# Write the combined data to a new TSV file
output_path = os.path.join(COUNT_DIR,'counts_data.tsv')
combined_data1.to_csv(output_path, sep='\t', index=False)

####################################    ExInAtor_Gene_List     ###############################################################
combined_data1 = pd.DataFrame()
# Loop through each folder in the directory
for folder in os.listdir(directory):
    folder_path = os.path.join(directory, folder)
    if os.path.isdir(folder_path):
        # Loop through each TSV file in the folder
        for file in os.listdir(folder_path):
            file_path = os.path.join(folder_path, file)
            if file == 'ExInAtor_Gene_List.txt':
                # Read the TSV file into a DataFrame and append it to the combined data
                data = pd.read_csv(file_path, sep='\t')
                data['Element'] = folder
                columns = data.columns.tolist()
                columns = ['Element'] + columns[:-1]
                data = data[columns]
                combined_data1 = pd.concat([combined_data1, data], axis=0, ignore_index=True)

# Write the combined data to a new TSV file
output_path = os.path.join(COUNT_DIR,'ExInAtor_Gene_List.tsv')
combined_data1.to_csv(output_path, sep='\t', index=False)

# read in the count file as a pandas dataframe
count_df = pd.read_csv(os.path.join(COUNT_DIR,"counts_data.tsv"), sep='\t')
# read in the count file as a pandas dataframe
kmer_df = pd.read_csv(os.path.join(COUNT_DIR,"kmer_counts_data.tsv"), sep='\t')


count_df=count_df.sort_values('Element', ascending=True).reset_index(drop=True)
count_df
count_df.to_csv(output_store,index=False)

f = count_df
f['log2OR'] = np.log2(f['OR-FET'])
f['p-value FET'] = f['p-value FET'].replace(0.000000e+00, 0.000000000000000001)
f['10logp'] = -np.log10(f['p-value FET'])
f.to_csv(output_store_na,index=False)