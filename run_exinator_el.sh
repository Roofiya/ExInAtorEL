#!/bin/bash
# Script to run ExInAtorEL

# Usage message for when the script is run incorrectly
if [ "$#" -lt 5 ]; then
    echo "Usage: $0 <element_file> <mutations_file> <gtf_file> <GENESET> <strand>"
    exit 1
fi

# Base directory that holds the shared input files (blacklist, genome, etc.)
# >>> INSERT YOUR PATH <<<
BASE_DIR="/path/to/base_dir"

# Root directory where all pipeline outputs will be written
# >>> INSERT YOUR PATH <<<
PIPELINE_DIR="/path/to/pipeline_outputs"

# Set input BED file and Mutations file names from the passed arguments
BED_FILE_NAME=$1
MUTATIONS_FILE_NAME=$2
GTF_NAME=$3
GENE_SET=$4
STRAND=$5

# Convert strand to directory name
if [[ "$STRAND" == "true" ]]; then
    STRAND_DIR="strand_specific"
elif [[ "$STRAND" == "false" ]]; then
    STRAND_DIR="unstranded"
else
    echo "Error: Strand must be either 'true' or 'false'."
    exit 1
fi

# Paths to Python version and script locations
PYTHON="python2.7"
# Python 3 interpreter from your post-processing conda env  >>> INSERT YOUR PATH <<<
PYTHON3="/path/to/conda/envs/2022_base/bin/python3"
MAIN_SCRIPT="Main.py"
PROCESS_SCRIPT="process_outputs.py"
BH_CORRECTION_SCRIPT="bh_correction.py"
HEATMAP_SCRIPT="heatmaps.py"

# Input files (BED, mutations and GTF come in as arguments; the rest are shared)
BED_FILE="${BED_FILE_NAME}"
MUTATION_FILE="${MUTATIONS_FILE_NAME}"
GTF_FILE="${GTF_NAME}"
BLACKLIST_FILE="${BASE_DIR}/path/to/blacklist.bed"
FASTA_FILE="${BASE_DIR}/path/to/genome.fa"
WHOLE_GTF="${BASE_DIR}/path/to/whole_genome_lncRNA.gtf"
CHROM_FILE="${BASE_DIR}/path/to/chromosomes.bed"
KMERS_FILE="${BASE_DIR}/path/to/3mers.txt"

# Use the BED file base name as the element type and output folder name
ELEMENT_TYPE=$(basename "$BED_FILE_NAME" .bed)

# Set output directories for each step
OUTPUT_DIR="${PIPELINE_DIR}/Outputs/${STRAND_DIR}/${GENE_SET}/${ELEMENT_TYPE}"
COUNT_DIR="${PIPELINE_DIR}/Count_Files/${STRAND_DIR}/${GENE_SET}"
BASE_OUTPUT_DIR=$(dirname "$OUTPUT_DIR")
BH_CORRECTED_DIR="${PIPELINE_DIR}/BH_corrected/${STRAND_DIR}/${GENE_SET}/"
BH_DIR="${PIPELINE_DIR}/BH_corrected/${STRAND_DIR}"
HEATMAP_DIR="${PIPELINE_DIR}/Heatmaps/${STRAND_DIR}"

# Parameters
CORES=6
N=119
BOOTSTRAPS=10000

# Create the output directories if they don't exist
mkdir -p "$OUTPUT_DIR"
mkdir -p "$COUNT_DIR"
mkdir -p "$BH_CORRECTED_DIR"
mkdir -p "$HEATMAP_DIR"

# Check that the BED and mutation files exist
if [ ! -f "$BED_FILE" ]; then
    echo "Error: The specified BED file '$BED_FILE' does not exist."
    exit 1
fi

if [ ! -f "$MUTATION_FILE" ]; then
    echo "Error: The specified mutations file '$MUTATION_FILE' does not exist."
    exit 1
fi

# Step 1: Run Main.py with the specified BED file to generate enrichment outputs
echo "Running Main.py with BED file $BED_FILE and Mutations file $MUTATION_FILE..."
$PYTHON "$MAIN_SCRIPT" -i "$MUTATION_FILE" -o "$OUTPUT_DIR" -f "$FASTA_FILE" -g "$GTF_FILE" -e "$BED_FILE" -r "$BLACKLIST_FILE" -s "$CHROM_FILE" -k "$KMERS_FILE" -w "$WHOLE_GTF" -c $CORES -n $N -b $BOOTSTRAPS -ss $STRAND
if [ $? -ne 0 ]; then
    echo "Error: Main.py execution failed."
    exit 1
fi
echo "Main.py completed successfully."

# Step 1a: Switch from the py27 environment to the Python 3 post-processing env
echo "Deactivating py27 environment and activating 2022_base environment..."
conda deactivate py27
conda activate 2022_base

echo "Using Python interpreter: $($PYTHON3 --version)"
which $PYTHON3

# Step 2: Process output folders to generate na_output_*.csv files
echo "Step 2: Processing output folders to generate na_output_*.csv files..."
$PYTHON3 "$PROCESS_SCRIPT" -i "$BASE_OUTPUT_DIR" -o "$COUNT_DIR" -g "$GENE_SET"
if [ $? -ne 0 ]; then
    echo "Error: process_outputs.py execution failed."
    exit 1
fi
echo "Output processing completed successfully."

# Step 3: Run bh_correction.py on the Count Files
echo "Running BH correction on files in $COUNT_DIR..."
$PYTHON3 "$BH_CORRECTION_SCRIPT" -i "$COUNT_DIR" -o "$BH_CORRECTED_DIR" -g "$GENE_SET"
if [ $? -ne 0 ]; then
    echo "Error: bh_correction.py execution failed."
    exit 1
fi
echo "BH correction completed successfully. Results saved in $BH_CORRECTED_DIR."

# Step 4: Run heatmaps.py on the BH corrected Files
echo "Running heatmaps.py on files in $BH_DIR..."
# Plot only the cohort from this run. Omit --cohorts (or run heatmaps.py by hand) to build a multi-cohort comparison from cohorts in $BH_DIR.
$PYTHON3 "$HEATMAP_SCRIPT" -i "$BH_DIR" -o "$HEATMAP_DIR" --cohorts "$GENE_SET"
if [ $? -ne 0 ]; then
    echo "Error: heatmaps.py execution failed."
    exit 1
fi
echo "Heatmaps generated successfully. Results saved in $HEATMAP_DIR."

# Completion message
echo "Pipeline execution finished. Outputs in $OUTPUT_DIR, count files in $BASE_OUTPUT_DIR, BH-corrected files in $BH_CORRECTED_DIR and heatmaps in $HEATMAP_DIR."