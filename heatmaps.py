import os
import argparse
from pathlib import Path
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)

# Set up argument parsing
parser = argparse.ArgumentParser(description="Cross-cohort enrichment heatmap.")
parser.add_argument("-i", "--bh_corrected_dir", required=True,
                    help="Directory containing one sub-folder per cohort, each with na_output_<cohort>.csv")
parser.add_argument("-o", "--heatmap_dir", required=True,
                    help="Output directory for the heatmap PDF.")
parser.add_argument("--cohorts", default=None,
                    help="Comma-separated cohort names to plot, in order. "
                         "If omitted, cohorts are auto-detected from the input directory.")
parser.add_argument("--labels", default=None,
                    help="Comma-separated display labels matching --cohorts. "
                         "If omitted, the cohort folder names are used.")
args = parser.parse_args()

BH_CORRECTED_DIR = args.bh_corrected_dir
HEATMAP_DIR = args.heatmap_dir

print("Input Directory:", BH_CORRECTED_DIR)
print("Output Directory:", HEATMAP_DIR)


def cohort_csv(cohort):
    """Path to a cohort's BH-corrected table."""
    return os.path.join(BH_CORRECTED_DIR, cohort, f'na_output_{cohort}.csv')


# --- Decide which cohorts to plot --------------------------------------------
# Explicit list via --cohorts, otherwise auto-detect every sub-folder that
# actually contains an na_output_<cohort>.csv file.
if args.cohorts:
    cancer_types = [c.strip() for c in args.cohorts.split(',') if c.strip()]
else:
    cancer_types = sorted(
        d.name for d in Path(BH_CORRECTED_DIR).iterdir()
        if d.is_dir() and os.path.exists(cohort_csv(d.name))
    )

if not cancer_types:
    raise SystemExit(
        f"No cohorts found in {BH_CORRECTED_DIR}. "
        f"Expected sub-folders each containing na_output_<cohort>.csv, "
        f"or pass --cohorts explicitly."
    )

print("Cohorts to plot:", cancer_types)

# --- Read each cohort, skipping any that are missing -------------------------
combined_log2OR_df = pd.DataFrame()
combined_pvalues_df = pd.DataFrame()
loaded = []

for cancer in cancer_types:
    csv_file_path = cohort_csv(cancer)
    if not os.path.exists(csv_file_path):
        print(f"WARNING: skipping '{cancer}' - file not found: {csv_file_path}")
        continue
    try:
        df = pd.read_csv(csv_file_path, usecols=['Element', 'log2OR', 'Corrected P-Value'])
    except Exception as e:
        print(f"WARNING: skipping '{cancer}' - could not read {csv_file_path}: {e}")
        continue
    if df.empty:
        print(f"WARNING: skipping '{cancer}' - no rows in {csv_file_path}")
        continue
    df['CancerType'] = cancer
    combined_log2OR_df = pd.concat([combined_log2OR_df, df], ignore_index=True)
    combined_pvalues_df = pd.concat([combined_pvalues_df, df], ignore_index=True)
    loaded.append(cancer)

if not loaded:
    raise SystemExit("No cohort tables could be loaded - nothing to plot.")

# --- Resolve display labels (default to the cohort names) --------------------
if args.labels:
    labels = [l.strip() for l in args.labels.split(',')]
    if len(labels) != len(loaded):
        print(f"WARNING: {len(labels)} labels for {len(loaded)} cohorts - "
              f"falling back to cohort names.")
        labels = loaded
else:
    labels = loaded

# --- Pivot: elements as rows, cohorts as columns -----------------------------
combined_log2OR_df_pivot = combined_log2OR_df.pivot(index='Element', columns='CancerType', values='log2OR')
combined_pvalues_df_pivot = combined_pvalues_df.pivot(index='Element', columns='CancerType', values='Corrected P-Value')

# Dictionary for renaming element identifiers to readable names
rename_dict = {
    'crss': 'Conserved RNA secondary structures',
    'rbp_combined': 'RNA-binding protein interaction sites',
    'ridls_strand': 'Repeat Insertion Domains',
    'rnadna': 'RNA:DNA interaction sites',
    'repeats': 'Repeat regions',
    'encori': 'MicroRNA response elements',
    'mircode.highcons': 'Highly conserved microRNA binding sites',
    'mircode.mediumcons': 'Medium conserved microRNA binding sites',
    '100way': 'Conservation 100 species',
    '30way': 'Conservation 30 species',
    '20way': 'Conservation 20 species',
    '7way': 'Conservation 7 species',
    'phastcons': 'Evolutionary conserved elements',
    'rnamod': 'RNA modification sites'
}
combined_log2OR_df_pivot = combined_log2OR_df_pivot.rename(index=rename_dict)
combined_pvalues_df_pivot = combined_pvalues_df_pivot.rename(index=rename_dict)

# Keep only cohorts that actually made it into the pivot, preserving order
loaded = [c for c in loaded if c in combined_log2OR_df_pivot.columns]
labels = labels[:len(loaded)] if args.labels else loaded
combined_log2OR_df_pivot = combined_log2OR_df_pivot[loaded]
combined_pvalues_df_pivot = combined_pvalues_df_pivot[loaded]

# Mask depleted (log2OR <= 0) or non-significant (q >= 0.1) cells
mask = (combined_log2OR_df_pivot <= 0) | (combined_pvalues_df_pivot >= 0.1)

# --- Plot --------------------------------------------------------------------
n_cols = combined_log2OR_df_pivot.shape[1]
n_rows = combined_log2OR_df_pivot.shape[0]
figSize = (max(1.1 * n_cols + 2, 4), max(0.3 * n_rows + 1, 3))
fig, ax = plt.subplots(figsize=figSize, gridspec_kw={'wspace': 0}, dpi=100)

sns.heatmap(combined_log2OR_df_pivot, mask=mask, annot=combined_pvalues_df_pivot,
            cbar_kws={'label': 'log2(odds)\nenrichment score'},
            square=False, cmap='Blues',
            vmax=combined_log2OR_df_pivot.max().max(), vmin=0,
            cbar=True, linewidths=.5, xticklabels=labels, linecolor='grey')
ax.set_title('Cohort-specific enrichment in RNA elements\n')
plt.xlabel('Cohort')
plt.xticks(rotation=45)
plt.ylabel('Elements')
for side in ('top', 'bottom', 'left', 'right'):
    ax.spines[side].set_visible(True)

os.makedirs(HEATMAP_DIR, exist_ok=True)
plt.savefig(f'{HEATMAP_DIR}/Motif_cancer_elements.pdf', format='pdf', bbox_inches='tight')
print(f"Saved: {HEATMAP_DIR}/Motif_cancer_elements.pdf")