import argparse
import os
import math
import pandas as pd
from tqdm import tqdm

parquet_records = []
signals_dir = None 

def process_group(meta_group, index, chrom, chrom_dir, group_id=None):
    signals = meta_group.index.tolist()
    mat_files = [os.path.join(signals_dir, f"{sig}.{input_type}") for sig in signals]
    
    min_loc = meta_group["location_min"].min()
    max_loc = meta_group["location_max"].max() if "location_max" in meta_group.columns else meta_group["location_min"].max()

    snp_set = set()
    for mat_file in mat_files:
        if input_type == "feather":
            df_tmp = pd.read_feather(mat_file)
        else:  
            df_tmp = pd.read_pickle(mat_file)
        snp_set.update(df_tmp.columns.tolist())
        del df_tmp

    columns = list(meta_group.columns) + sorted(snp_set)
    combined_df = pd.DataFrame(index=meta_group.index, columns=columns)
    for col in meta_group.columns:
        combined_df[col] = meta_group[col]
    combined_df.iloc[:, len(meta_group.columns):] = -1e6

    combined_array = combined_df.to_numpy()
    snp_columns = {snp: idx for idx, snp in enumerate(combined_df.columns[len(meta_group.columns):], start=len(meta_group.columns))}

    for mat_file in mat_files:
        signal_name = os.path.splitext(os.path.basename(mat_file))[0]
        if input_type == "feather":
            df_mat = pd.read_feather(mat_file)
        else:
            df_mat = pd.read_pickle(mat_file)
        row_idx = combined_df.index.get_loc(signal_name)
        for snp_col, value in zip(df_mat.columns, df_mat.iloc[0].values):
            if snp_col in snp_columns:
                combined_array[row_idx, snp_columns[snp_col]] = value
        del df_mat

    combined_df = pd.DataFrame(combined_array, index=combined_df.index, columns=combined_df.columns)
    combined_df.reset_index(inplace=True)
    
    if group_id is not None:
        parquet_filename = f"chr{chrom}_group_{group_id}.parquet"
    else:
        parquet_filename = f"chr{chrom}_met_group_{index}_region_{min_loc}-{max_loc}.parquet"
    parquet_path = os.path.join(chrom_dir, parquet_filename)
    
    combined_df.to_parquet(parquet_path)
    parquet_records.append({
        "chromosome": chrom,
        "group": group_id if group_id is not None else index,
        "n_signals": combined_df.shape[0],
        "min_position": min_loc,
        "max_position": max_loc,
        "parquet_file": parquet_path
    })
    
    return index + 1

def create_parquet(meta_sub, index, chrom, chrom_dir):
    meta_sub.sort_values(by="location_min", inplace=True)
    positions = meta_sub["location_min"].tolist()

    if len(positions) >= 2:
        positions_sorted = sorted(positions)
        for i in range(len(positions_sorted)):
            gap = positions_sorted[i] - positions_sorted[i-1]
            if gap > 1_000_000:
                split_point = positions_sorted[i - 1] 
                df_part1 = meta_sub[meta_sub["location_min"] <= split_point].copy()
                df_part2 = meta_sub[meta_sub["location_min"] > split_point].copy()
                index = create_parquet(df_part1, index, chrom, chrom_dir)
                index = create_parquet(df_part2, index, chrom, chrom_dir)
                return index

    if len(meta_sub) > 1000:
        signals = meta_sub.index.tolist()
        total = len(signals)
        chunk_size = 1000
        n_groups = math.ceil(total / chunk_size)
        for group in range(n_groups):
            start = group * chunk_size
            end = min(start + chunk_size, total)
            meta_group = meta_sub.loc[signals[start:end]].copy()
            index = process_group(meta_group, index, chrom, chrom_dir, group_id=index)
        return index

    index = process_group(meta_sub, index, chrom, chrom_dir, group_id=index)
    return index


def main():
    parser = argparse.ArgumentParser(
        description="Process signals with recursive gap splitting (>500k) and chunking (max 1000 signals)"
    )
    parser.add_argument("--input", type=str, required=True, help="Directory containing signal pickle files")
    parser.add_argument("--output", type=str, required=True, help="Directory to save parquet files")
    parser.add_argument("--input_summary", type=str, required=True, help="Path to summary TSV file")
    parser.add_argument("--output_summary", type=str, help="Path to write parquet summary TSV")
    parser.add_argument("--input_type", type=str, help="Type of input files ('pickle' or 'feather')", default="pickle", choices=["pickle", "feather"])

    args = parser.parse_args()

    global signals_dir
    signals_dir = args.input

    global input_type
    input_type = args.input_type

    os.makedirs(args.output, exist_ok=True)
    metadata = pd.read_csv(args.input_summary, sep="\t")
    metadata["chromosome"] = metadata["chromosome"].astype(str)
    chromosomes = metadata["chromosome"].unique()

    for chrom in tqdm(chromosomes, desc="Processing chromosomes"):
        group_index = 0
        chrom_dir = os.path.join(args.output, chrom)
        os.makedirs(chrom_dir, exist_ok=True)
        meta_sub = metadata[metadata["chromosome"] == chrom].copy()
        meta_sub.set_index("signal", inplace=True)
        meta_sub.sort_values(by="location_min", inplace=True)
        group_index = create_parquet(meta_sub, group_index, chrom, chrom_dir)

    if args.output_summary:
        pd.DataFrame(parquet_records).to_csv(args.output_summary, sep="\t", index=False)
    print("Done.")

if __name__ == "__main__":
    main()
