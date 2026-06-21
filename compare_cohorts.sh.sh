#!/bin/bash
# Make a multi-cohort comparison heatmap without retyping paths.
#
# Usage:
#   ./compare_cohorts.sh cohort1 cohort2 cohort3 ...
#
# Example:
#   ./compare_cohorts.sh AC BLCA BRCA
#
# Edit the two paths below once to match your setup (same values as
# BH_DIR and HEATMAP_DIR in run_exinator_el.sh).

# >>> INSERT YOUR PATHS <<<
BH_DIR="/path/to/pipeline_outputs/BH_corrected/strand_specific"
HEATMAP_DIR="/path/to/pipeline_outputs/Heatmaps/strand_specific"

PYTHON3="python3"          # or the full path to your 2022_base python3
HEATMAP_SCRIPT="heatmaps.py"

if [ "$#" -lt 1 ]; then
    echo "Usage: $0 <cohort1> [cohort2] [cohort3] ..."
    echo "  (no cohorts given -> plots ALL cohorts found in $BH_DIR)"
    exit 1
fi

# Join the cohort arguments into a comma-separated list
COHORTS=$(IFS=,; echo "$*")

echo "Comparing cohorts: $COHORTS"
$PYTHON3 "$HEATMAP_SCRIPT" -i "$BH_DIR" -o "$HEATMAP_DIR" --cohorts "$COHORTS"