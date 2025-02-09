import os
import math
import pandas as pd
from tqdm import tqdm

metadata = pd.read_csv("metadata.tsv", sep="\t")

parquet_dir = "p"
mat_dir = "mat"

os.makedirs(parquet_dir, exist_ok=True)

parquet_records = []

metadata["chromosome"] = metadata["chromosome"].astype(str)

chromosomes = metadata["chromosome"].unique()

for chrom in tqdm(chromosomes, desc="Processing chromosomes"):
    chrom_dir = f"{parquet_dir}/{chrom}"
    os.makedirs(chrom_dir, exist_ok=True)

    meta_sub = metadata[
        (metadata["chromosome"] == chrom) & 
        (metadata["signal_strength"] > 7)
    ].copy()

    meta_sub.set_index("signal", inplace=True)

    meta_sub.sort_values(by="location_min", ascending=False, inplace=True)

    signals = meta_sub.index.tolist()
    total_signals = len(signals)
    print(f"Chrom {chrom}: {total_signals} signals pass filter.")

    chunk_size = 1000
    n_groups = math.ceil(total_signals / chunk_size)

    for group_i in range(n_groups):
        start_idx = group_i * chunk_size
        end_idx = min(start_idx + chunk_size, total_signals)
        group_signals = signals[start_idx:end_idx]

        print(f"  Group {group_i+1}/{n_groups}, signals {start_idx}..{end_idx-1}")

        meta_group = meta_sub.loc[group_signals].copy()
        mat_files = [f"{mat_dir}/{sig}.pickle" for sig in group_signals]

        min_loc = meta_group["location_min"].min()
        max_loc = meta_group["location_min"].max()  

        snp_set = set()
        for mat_file in tqdm(mat_files, desc=f"Collecting SNPs in group {group_i+1}", leave=False):
            df_tmp = pd.read_pickle(mat_file)
            snp_set.update(df_tmp.columns.tolist())
            del df_tmp

        print(f"    => {len(snp_set)} unique SNPs in group {group_i+1}")

        columns = list(meta_group.columns) + sorted(snp_set)
        combined_df = pd.DataFrame(index=meta_group.index, columns=columns)

        for col in meta_group.columns:
            combined_df[col] = meta_group[col]

        combined_df.iloc[:, len(meta_group.columns):] = -1e+6

        combined_array = combined_df.to_numpy()
        snp_columns = {
            snp: idx 
            for idx, snp in enumerate(combined_df.columns[len(meta_group.columns):], start=len(meta_group.columns))
        }

        for mat_file in tqdm(mat_files, desc=f"Reading SNP values in group {group_i+1}", leave=False):
            signal_name = os.path.splitext(os.path.basename(mat_file))[0]  
            df_mat = pd.read_pickle(mat_file)
            
            row_idx = combined_df.index.get_loc(signal_name)
            for snp_col, value in zip(df_mat.columns, df_mat.iloc[0].values):
                if snp_col in snp_columns:
                    combined_array[row_idx, snp_columns[snp_col]] = value

            del df_mat

        combined_df = pd.DataFrame(combined_array, index=combined_df.index, columns=combined_df.columns)
        combined_df.reset_index(inplace=True)

        parquet_filename = f"chr{chrom}_group_{group_i+1}.parquet"
        parquet_path = os.path.join(chrom_dir, parquet_filename)

        combined_df.to_parquet(parquet_path, engine="pyarrow")

        parquet_records.append({
            "chromosome": chrom,
            "group": group_i + 1,
            "n_signals": len(group_signals),
            "min_position": min_loc,
            "max_position": max_loc,
            "parquet_file": parquet_path
        })

df_parquet_meta = pd.DataFrame(parquet_records)
df_parquet_meta.to_csv("parquet_metadata.tsv", sep="\t", index=False)

print("Done.")
