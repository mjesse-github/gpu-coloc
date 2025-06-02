import os
import pandas as pd
from tqdm import tqdm
import argparse

def process_gwas_file(file_path, signals_directory, summary_file_path):

    try:
        # df = pd.read_csv(file_path, sep='\t', low_memory=False, on_bad_lines='skip')
        df = pd.read_csv(file_path, compression='gzip', sep='\t') 

    except Exception as e:
        print(f"Error reading file {file_path}: {e}")
        return
    
    df = df[df["maf"]>= 0.01]
    df["variant"] = df["chromosome"]+"_"+df["position"].astype(str)+"_"+df["allele1"]+"_"+df["allele2"]
    trait_data = df.groupby(df['trait'])
    
    # for trait_id, group in tqdm(trait_data, desc=f"Processing molecular traits for {file_path}", leave=False):
    for trait_id, group in trait_data:
        for region in group["region"].unique():

            region_df = group[group["region"]==region]

            lbfs = [column for column in group.columns if "lbf_variable" in column]

            for lbf_column in lbfs:
                lbf_index = int(lbf_column[12:])
                df_filtered = region_df[['trait', 'region', 'variant', 'chromosome', 'position', lbf_column]].copy()
                df_filtered.sort_values(by="position", ascending=True, inplace=True)

                df_filtered.rename(columns={lbf_column: 'lbf'}, inplace=True)

                if (df_filtered['lbf'] == 0).all():
                    continue

                signal_strength = df_filtered['lbf'].max()
                if signal_strength < 5:
                    continue
                
                chromosome = df_filtered['chromosome'].iloc[0]
                location_min = df_filtered['position'].min()
                location_max = df_filtered['position'].max()

                signal = f"{trait_id}_{region}_L{lbf_index}"

                output_file_name = f"{signal}.pickle"
                output_file_path = os.path.join(signals_directory, output_file_name)

                df_filtered.loc[:, 'lbf'] = pd.to_numeric(df_filtered['lbf'], errors='coerce')

                mat_df = pd.DataFrame(df_filtered.set_index('variant')['lbf']).T

                variant_id = mat_df.T["lbf"].idxmax()

                mat_df.to_pickle(output_file_path)

                summary_data = pd.DataFrame([{
                    'signal': signal,
                    'chromosome': chromosome[3:],
                    'location_min': location_min,
                    'location_max': location_max,
                    'signal_strength': signal_strength,
                    'lead_variant': variant_id
                }])

                header_needed = not os.path.exists(summary_file_path)  
                summary_data.to_csv(summary_file_path, sep='\t', mode='a', header=header_needed, index=False)

parser = argparse.ArgumentParser(description="Seperate and map GWAS signals")

parser.add_argument("--summary", type=str, required=True, help="Path where to write summary, e.g., 'gwas_summary.tsv'.")
parser.add_argument("--output", type=str, required=True, help="Dir where to write the signals, e.g., 'gwas_signals'.")
parser.add_argument("--input_table", type=str, required=True, help="table of inputs, e.g., 'gwas_table.tsv'.")

args = parser.parse_args()

signals_directory = args.output
summary_file_path = args.summary

os.makedirs(signals_directory, exist_ok=True)

input_table = pd.read_csv(args.input_table, sep="\t")

for gwas in tqdm(input_table.itertuples(index=False), desc=f"Processing GWAS's"):
    process_gwas_file(gwas.path, signals_directory, summary_file_path)

