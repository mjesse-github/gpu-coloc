import os
import pandas as pd
from tqdm import tqdm

def process_gwas_file(file_path, directory_name, mat_dir, summary_file_path):
    signals_dir = os.path.join(mat_dir, "mat")
    os.makedirs(signals_dir, exist_ok=True)
    
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

    if signal_strength < 1:
        return
    
    chromosome = df_filtered['chromosome'].iloc[0]
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

metadata_path = 'tabix_ftp_paths.tsv'  # Adjust this if the file path changes
metadata = pd.read_csv(metadata_path, sep='\t')  # Assuming TSV format; adjust if different

metadata = metadata[metadata["quant_method"].isin(["ge", "microarray"])]

for _, row in tqdm(metadata.iterrows(), desc="processing QTD"):
    r7 = ['QTS000032', 'QTS000033', 'QTS000034', 'QTS000036', 'QTS000037', 'QTS000038', 'QTS000039', 'QTS000040', 'QTS000041', 'QTS000042']

    if row["study_id"] in r7:
        QTS_directory = "eQTL_Catalogue_r7"
    else: 
        QTS_directory = "eQTL_Catalogue_r6"

    QTS = row["study_id"]
    QTD = row["dataset_id"]
    QTS_PATH = f"/gpfs/space/projects/eQTLCatalogue/qtlmap/{QTS_directory}/susie/{QTS}/{QTD}/{QTD}.lbf_variable.txt.gz"

    mat_dir = "mat"
    summary_file_path = "metadata.tsv"

    process_gwas_file(QTS_PATH, QTD, mat_dir, summary_file_path)
