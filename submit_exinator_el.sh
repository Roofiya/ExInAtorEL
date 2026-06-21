#!/bin/bash -l
#BSUB -q short
#BSUB -P re_gecip_cancer_pan
#BSUB -n 1
#BSUB -W 2:40
#BSUB -J GEL[1-252] #Change based on elements analysed
#BSUB -o /path/to/logs/job.%J.%I.out
#BSUB -e /path/to/logs/job.%J.%I.err
#BSUB -R "rusage[mem=10000] span[hosts=1]"
#BSUB -M 10000

BED_FILES=($(ls /path/to/element_bed_files/your_element_set/*.bed))

cd /path/to/your_working_dir/

module load R/4.3.3
conda activate py27

BED_FILE=${BED_FILES[$((LSB_JOBINDEX - 1))]}

BASE_DIR="/path/to/base_dir"
MUTATION_FILE="${BASE_DIR}/path/to/your_mutations_sorted_uniq.bed"
STRAND=true
GTF="${BASE_DIR}/path/to/your_annotation.gtf"
GENESET="your_cohort_name"

if [ -z "$GENESET" ] || [ -z "$MUTATION_FILE" ] || [ -z "$GTF" ]; then
    echo "Error: Missing arguments. Usage: $0 GENESET MUTATION_FILE GTF"
    exit 1
fi

if [ ! -f "$BED_FILE" ]; then
    echo "Error: BED file not found at $BED_FILE"
    exit 1
fi

if [ ! -f "$GTF" ]; then
    echo "Error: GTF file not found at $GTF"
    exit 1
fi

if [ ! -f "$MUTATION_FILE" ]; then
    echo "Error: Mutation file not found for $BED_FILE"
    exit 1
fi

./run_exinator_el.sh "$BED_FILE" "$MUTATION_FILE" "$GTF" "$GENESET" "$STRAND"
