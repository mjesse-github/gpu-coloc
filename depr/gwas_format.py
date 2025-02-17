import argparse
import os
import pandas as pd
from tqdm import tqdm

parquet_records = []

def create_parquet(meta_sub, index, chrom, chrom_dir):
    max_gap = 0
    i1 = 0
    j1 = 0

    positions = meta_sub['location_min'].tolist()
    positions_sorted = sorted(positions)

    for i in range(len(positions_sorted)-1):
        gap = positions_sorted[i+1] - positions_sorted[i]
        if gap > max_gap:
            max_gap = gap
            i1 = i
            j1 = i+1

    if max_gap > 500_000:
        split_point_1 = positions_sorted[i1]
        split_point_2 = positions_sorted[j1]

        df_part1 = meta_sub[meta_sub["location_min"] <= split_point_1].copy()
        df_part2 = meta_sub[meta_sub["location_min"] >= split_point_2].copy()

        index = create_parquet(df_part1, index, chrom, chrom_dir)      
        index = create_parquet(df_part2, index, chrom, chrom_dir)  

        return index

    signals = meta_sub.index.tolist()
    mat_files = [f"{signals_dir}/{sig}.pickle" for sig in signals]

    snp_set = set()
    for mat_file in mat_files:
        df_tmp = pd.read_pickle(mat_file)
        snp_set.update(df_tmp.columns.tolist())
        del df_tmp

    columns = list(meta_sub.columns) + sorted(snp_set)
    combined_df = pd.DataFrame(index=meta_sub.index, columns=columns)

    for col in meta_sub.columns:
        combined_df[col] = meta_sub[col]

    combined_df.iloc[:, len(meta_sub.columns):] = -1e+6

    combined_array = combined_df.to_numpy()
    snp_columns = {
        snp: idx 
        for idx, snp in enumerate(combined_df.columns[len(meta_sub.columns):], start=len(meta_sub.columns))
    }

    for mat_file in mat_files:
        signal_name = os.path.splitext(os.path.basename(mat_file))[0]  
        df_mat = pd.read_pickle(mat_file)
        
        row_idx = combined_df.index.get_loc(signal_name)
        for snp_col, value in zip(df_mat.columns, df_mat.iloc[0].values):
            if snp_col in snp_columns:
                combined_array[row_idx, snp_columns[snp_col]] = value

        del df_mat

    combined_df = pd.DataFrame(combined_array, index=combined_df.index, columns=combined_df.columns)
    combined_df.reset_index(inplace=True)

    min_loc = combined_df["location_min"].min()
    max_loc = combined_df["location_max"].max()

    parquet_filename = f"chr{chrom}_met_group_{index}_region_{min_loc}-{max_loc}.parquet"
    parquet_path = os.path.join(chrom_dir, parquet_filename)

    combined_df.to_parquet(parquet_path, engine="pyarrow")

    parquet_records.append({
        "chromosome": chrom,
        "group": index,
        "n_signals": combined_df.shape[0],
        "min_position": min_loc,
        "max_position": max_loc,
        "parquet_file": parquet_path
    })

    return index + 1

parser = argparse.ArgumentParser(description="Format GWAS signals")

parser.add_argument("--input", type=str, required=True, help="Directory where are signal files, e.g., 'eqtl_signals'.")
parser.add_argument("--output", type=str, required=True, help="Directory where to save formatted matrixes, e.g., 'formatted_eqtls'.")
parser.add_argument("--input_summary", type=str, required=True, help="Path to the summary eQTL summary file, e.g., 'eqtl_summary.tsv'.")
parser.add_argument("--output_summary", type=str, help="Path, e.g., 'parquet_summary.tsv'.")

args = parser.parse_args()

summary_file_path = args.input_summary 
signals_dir = args.input
output_dir = args.output

metadata = pd.read_csv(summary_file_path, sep="\t")
os.makedirs(output_dir, exist_ok=True)

metadata["chromosome"] = metadata["chromosome"].astype(str)
chromosomes = metadata["chromosome"].unique()

for chrom in tqdm(chromosomes, desc="Processing chromosomes"):
    chrom_dir = f"{output_dir}/{chrom}"
    os.makedirs(chrom_dir, exist_ok=True)

    meta_sub = metadata[
        (metadata["chromosome"] == chrom) #& 
        # (metadata["signal_strength"] > 7)
    ].copy()

    meta_sub.set_index("signal", inplace=True)
    meta_sub.sort_values(by="location_min", ascending=True, inplace=True)

    create_parquet(meta_sub, 0, chrom, chrom_dir)

if args.output_summary:
    df_parquet_meta = pd.DataFrame(parquet_records)
    df_parquet_meta.to_csv(args.output_summary, sep="\t", index=False)

print("Done.")
