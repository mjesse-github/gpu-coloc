import os
import pandas as pd
from tqdm import tqdm
import argparse

def process_gwas_file(file_path, directory_name, mat_dir, summary_file_path):
    os.makedirs(mat_dir, exist_ok=True)

    try:
        df = pd.read_csv(file_path, sep='\t', low_memory=False, on_bad_lines='skip')
    except Exception as e:
        print(f"Error reading file {file_path}: {e}")
        return

    trait_data = df.groupby(df['molecular_trait_id'])
    
    for trait_id, group in tqdm(trait_data, desc=f"Processing molecular traits for {file_path}"):
        for i in range(1, 11):
            lbf_column = f'lbf_variable{i}'
            if lbf_column in group.columns:
                process_signal(group, directory_name, trait_id, mat_dir, lbf_column, i, summary_file_path)

def process_signal(group, directory_name, trait_id, mat_dir, lbf_column, lbf_index, summary_file_path):
    df_filtered = group[['molecular_trait_id', 'region', 'variant', 'chromosome', 'position', lbf_column]].copy()
    df_filtered.rename(columns={lbf_column: 'lbf'}, inplace=True)
    
    if (df_filtered['lbf'] == 0).all():
        return
    
    signal_strength = df_filtered['lbf'].abs().max()

    if signal_strength < 5:
        return
        
    chromosome = df_filtered['chromosome'].astype(str).iloc[0]

    if chromosome != "1":
        return

    location_min = df_filtered['position'].min()
    location_max = df_filtered['position'].max()

    signal = f"{directory_name}_{trait_id}_L{lbf_index}"
    output_file_name = f"{signal}.pickle"
    output_file_path = os.path.join(mat_dir, output_file_name)

    df_filtered['lbf'] = pd.to_numeric(df_filtered['lbf'], errors='coerce')

    mat1_df = pd.DataFrame(df_filtered.set_index('variant')['lbf']).T
    variant_id = mat1_df.T["lbf"].idxmax()

    mat1_df.to_pickle(output_file_path)

    summary_data = pd.DataFrame([{
        'signal': signal,
        'chromosome': chromosome,
        'location_min': location_min,
        'location_max': location_max,
        'signal_strength': signal_strength,
        'lead_variant': variant_id
    }])
    
    header_needed = not os.path.exists(summary_file_path)  
    summary_data.to_csv(summary_file_path, sep='\t', mode='a', header=header_needed, index=False)

parser = argparse.ArgumentParser(description="Create files and summary for signals")

parser.add_argument("--input_table", type=str, required=True, help="Path to TSV file containing 'id' and 'path' columns, e.g., 'eqtl_paths.tsv'.")
parser.add_argument("--output", type=str, required=True, help="Directory to store signal files, e.g., 'eqtl_signals'.")
parser.add_argument("--summary", type=str, required=True, help="File to write the summary, e.g., 'eqtl_summary.tsv'.")

args = parser.parse_args()

metadata = pd.read_csv(args.input_table, sep='\t')

# Process each file
for row in tqdm(metadata.itertuples(index=False), desc="Processing files"):
    dataset_id = row.id  # Dataset ID, e.g., QTD1
    file_path = row.path  # File path

    process_gwas_file(file_path, dataset_id, args.output, args.summary)