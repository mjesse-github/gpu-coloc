import argparse
import math
import os
import torch
import pandas as pd
import numpy as np
from tqdm import tqdm
import pyarrow.parquet as pq

def logdiff_torch(a, b):
    mx = torch.maximum(a, b)
    val = torch.exp(a - mx) - torch.exp(b - mx)
    mask = (val <= 0)
    out = mx + torch.log(torch.where(mask, torch.tensor(float('nan'), device=a.device), val))
    out[mask] = float('nan')
    return out

def coloc_bf_bf_torch(
    bf1_cpu, bf2_cpu,
    p1=1e-4, p2=1e-4, p12=5e-6,
    device="mps"
):
    if isinstance(bf1_cpu, pd.Series):
        bf1_cpu = bf1_cpu.to_frame().T
    if isinstance(bf2_cpu, pd.Series):
        bf2_cpu = bf2_cpu.to_frame().T

    isnps = list(set(bf1_cpu.columns).intersection(bf2_cpu.columns) - {"null"})
    if not isnps:
        return {
            "summary": pd.DataFrame({"nsnps": [np.nan]}),
            "pp_3d": None,
            "pp_H4_matrix": None,
            "priors": {"p1": p1, "p2": p2, "p12": p12}
        }

    bf1_arr = torch.tensor(bf1_cpu[isnps].values, dtype=torch.float32, device=device)
    bf2_arr = torch.tensor(bf2_cpu[isnps].values, dtype=torch.float32, device=device)

    N, M = bf1_arr.shape
    K, _ = bf2_arr.shape

    bf1_3d = bf1_arr.unsqueeze(1) 
    bf2_3d = bf2_arr.unsqueeze(0) 
    sum_3d = bf1_3d + bf2_3d 

    sum_3d_logexp = torch.logsumexp(sum_3d, dim=2)  

    l1_sum = torch.logsumexp(bf1_arr, dim=1)  
    l2_sum = torch.logsumexp(bf2_arr, dim=1)  

    l1_sum_2d = l1_sum.unsqueeze(1).expand(N, K)
    l2_sum_2d = l2_sum.unsqueeze(0).expand(N, K)

    p1_t  = torch.tensor(p1,  dtype=torch.float32, device=device)
    p2_t  = torch.tensor(p2,  dtype=torch.float32, device=device)
    p12_t = torch.tensor(p12, dtype=torch.float32, device=device)

    lH0_2d = torch.zeros((N, K), device=device)
    lH1_2d = torch.log(p1_t) + l1_sum_2d
    lH2_2d = torch.log(p2_t) + l2_sum_2d
    lH4_2d = torch.log(p12_t) + sum_3d_logexp
    lH3_2d = torch.log(p1_t) + torch.log(p2_t) + logdiff_torch(l1_sum_2d + l2_sum_2d, sum_3d_logexp)

    all_abf_3d = torch.stack([lH0_2d, lH1_2d, lH2_2d, lH3_2d, lH4_2d], dim=0) 
    denom_2d   = torch.logsumexp(all_abf_3d, dim=0)                           
    pp_abf_3d  = torch.exp(all_abf_3d - denom_2d.unsqueeze(0))             

    pp_H3_2d = pp_abf_3d[3]  
    pp_H4_2d = pp_abf_3d[4]  

    i_coords = torch.arange(N, device=device).unsqueeze(1).expand(N, K).flatten() 
    j_coords = torch.arange(K, device=device).unsqueeze(0).expand(N, K).flatten() 

    pp_H3_flat        = pp_H3_2d.flatten()
    pp_H4_flat        = pp_H4_2d.flatten()

    i_coords_cpu    = i_coords.cpu().numpy()
    j_coords_cpu    = j_coords.cpu().numpy()
    pp_H3_cpu       = pp_H3_flat.cpu().numpy()
    pp_H4_cpu       = pp_H4_flat.cpu().numpy()

    summary_df = pd.DataFrame({
        "idx1": i_coords_cpu,
        "idx2": j_coords_cpu,
        "PP.H3": pp_H3_cpu,
        "PP.H4": pp_H4_cpu,
    })

    return {
        "summary": summary_df,
        "pp_3d": pp_abf_3d.cpu().numpy(),
        "pp_H4_matrix": pp_H4_2d.cpu().numpy(),
        "priors": {"p1": p1, "p2": p2, "p12": p12}
    }


def logsum(arr):
    max_val = np.max(arr)
    return max_val + np.log(np.sum(np.exp(arr - max_val)))

def logbf_to_pp(df, pi, last_is_null):
    n = df.shape[1] - 1 if last_is_null else df.shape[1]
    
    if isinstance(pi, (int, float)):
        if pi > 1 / n:
            pi = 1 / n
        pi = np.append(np.repeat(pi, n), 1 - n * pi) if last_is_null else np.repeat(pi, n)
    
    if any(pi == 0):
        pi[pi == 0] = 1e-16
        pi = pi / np.sum(pi)
    
    if last_is_null:
        df = df.subtract(df.iloc[:, -1], axis=0)
    
    priors = np.tile(np.log(pi), (df.shape[0], 1))
    
    denom = np.apply_along_axis(logsum, 1, df.values + priors)
    
    denom_df = pd.DataFrame(np.tile(denom, (df.shape[1], 1)).T, index=df.index, columns=df.columns)
    
    result = np.exp(df.values + priors - denom_df.values)
    
    return pd.DataFrame(result, index=df.index, columns=df.columns)

def trim(bf1, bf2, p1=1e-4, p2=1e-4, overlap_min=0.5, silent=True):
    if isinstance(bf1, pd.Series):
        bf1 = bf1.to_frame().T
    if isinstance(bf2, pd.Series):
        bf2 = bf2.to_frame().T

    isnps = list(set(bf1.columns).intersection(set(bf2.columns)).difference(['null']))

    if not isnps:
        if not silent:
            print("No common SNPs found.")
        return pd.DataFrame({'nsnps': [np.nan]})
    
    pp1 = logbf_to_pp(bf1, p1, last_is_null=True)
    pp2 = logbf_to_pp(bf2, p2, last_is_null=True)

    bf1 = bf1[isnps]
    bf2 = bf2[isnps]

    prop1 = pp1[isnps].sum(axis=1) / pp1.loc[:, pp1.columns != "null"].sum(axis=1)

    prop2 = pp2[isnps].sum(axis=1) / pp2.loc[:, pp2.columns != "null"].sum(axis=1)

    todo = pd.DataFrame([(i, j) for i in range(bf1.shape[0]) for j in range(bf2.shape[0])], columns=['i', 'j'])

    drop = [prop1[todo['i'][k]] < overlap_min or prop2[todo['j'][k]] < overlap_min for k in range(len(todo))]

    if all(drop):
        if not silent:
            print("Warning: SNP overlap too small between datasets: too few SNPs with high posterior in one trait represented in other")

        return pd.DataFrame({'nsnps': [np.nan]})

    return todo[~pd.Series(drop)].reset_index(drop=True)

def coloc_loop(
    mat1: pd.DataFrame,
    mat2: pd.DataFrame,
    metadata1: pd.DataFrame,
    metadata2: pd.DataFrame,
    num_chunks1=0,
    num_chunks2=0,
    device="cuda",
    p1=1e-4, p2=1e-4, p12=1e-6, H4_threshold=0.8
):

    try:
        overlapping_pairs = trim(mat1, mat2)
        valid_pairs = set(overlapping_pairs[["i", "j"]].itertuples(index=False, name=None))
    except:
        print("Possible error in trim function")
        return pd.DataFrame()

    if overlapping_pairs.empty:
        return pd.DataFrame()
    
    chunk_size = 100

    mat1_chunks = []
    meta1_chunks = []
    start1_idx = 0

    N1 = mat1.shape[0]

    for i in range(num_chunks1):
        end1_idx = start1_idx + chunk_size if i < (num_chunks1 - 1) else N1
        mat1_chunk = mat1.iloc[start1_idx:end1_idx, :].copy()
        meta1_chunk = metadata1.iloc[start1_idx:end1_idx, :].copy()

        mat1_chunks.append(mat1_chunk)
        meta1_chunks.append(meta1_chunk)
        start1_idx = end1_idx 

    mat2_chunks = []
    meta2_chunks = []
    start2_idx = 0

    N2 = mat2.shape[0]

    for i in range(num_chunks2):
        end2_idx = start2_idx + chunk_size if i < (num_chunks2 - 1) else N2
        mat2_chunk = mat2.iloc[start2_idx:end2_idx, :].copy()
        meta2_chunk = metadata2.iloc[start2_idx:end2_idx, :].copy()

        mat2_chunks.append(mat2_chunk)
        meta2_chunks.append(meta2_chunk)
        start2_idx = end2_idx

    all_results = []

    total_pairs = []

    for i in range(num_chunks1):
        for j in range(num_chunks2):
            total_pairs.append((i, j))

    for pair in tqdm(total_pairs, desc="All chunk pairs", leave=False):    
        out = coloc_bf_bf_torch(
            bf1_cpu=mat1_chunks[pair[0]],
            bf2_cpu=mat2_chunks[pair[1]],
            p1=p1, p2=p2, p12=p12,
            device=device
        )
        if out is None or out["summary"] is None:
            continue

        summary_df = out["summary"]

        summary_df.loc[:, "idx1"] = summary_df["idx1"] + pair[0] * 100
        summary_df.loc[:, "idx2"] = summary_df["idx2"] + pair[1] * 100

        summary_df = summary_df[summary_df.apply(lambda row: (row["idx1"], row["idx2"]) in valid_pairs, axis=1)]

        summary_df = summary_df[summary_df["PP.H4"] >= H4_threshold].reset_index(drop=True)

        if summary_df.empty:
            continue

        summary_df["signal1"] = metadata1["signal"].iloc[
            summary_df["idx1"]
        ].values

        summary_df["lead1"] = metadata1["lead_variant"].iloc[
            summary_df["idx1"]
        ].values


        summary_df["signal2"] = metadata2["signal"].iloc[
            summary_df["idx2"]
        ].values

        summary_df["lead2"] = metadata2["lead_variant"].iloc[
            summary_df["idx2"]
        ].values

        summary_df = summary_df[summary_df["signal1"] != summary_df["signal2"]].reset_index(drop=True)

        summary_df.drop(columns=["idx1", "idx2"], inplace=True)

        all_results.append(summary_df)

    if all_results:
        final_df = pd.concat(all_results, ignore_index=True)
    else:
        final_df = pd.DataFrame()

    return final_df

def main():
    parser = argparse.ArgumentParser(description="Run coloc")

    parser.add_argument("--dir1", type=str, required=True, help="First directory of directories of parquet files, e.g., 'formatted_eqtls'.")
    parser.add_argument("--dir2", type=str, required=True, help="Second directory of directories of parquet files, e.g., 'formatted_metabolites'.")
    parser.add_argument("--results", type=str, required=True, help="File to write the colocalization results, e.g., 'results.tsv'.")
    parser.add_argument("--p12", type=float, required=True, help="p12 prior, e.g. 1e-6")
    parser.add_argument("--H4", type=float, required=False, help="Threshold for H4, e.g. 0.8", default=0.8)

    args = parser.parse_args()

    p12 = args.p12
    H4_threshold = args.H4

    if torch.cuda.is_available():
        device = torch.device("cuda")
    elif torch.backends.mps.is_available():
        device = torch.device("mps")
    else:
        device = torch.device("cpu")

    for root, dirs, files in os.walk(args.dir1):
        for directory in tqdm(dirs, desc="chromosomes"):
            try:
                dir_path = os.path.join(root, directory)
                met_files = os.listdir(dir_path)

                ge_dir_path = os.path.join(args.dir2, directory)
                files = os.listdir(ge_dir_path)

                met_cache = {}
                ge_cache = {}

                for i in range(len(met_files)):
                    # met_cache[i] = pd.read_parquet(f"{dir_path}/{met_files[i]}")
                    pf = pq.ParquetFile(
                        f"{dir_path}/{met_files[i]}",
                        thrift_string_size_limit=2**31-1,
                        thrift_container_size_limit=2**31-1,
                    )

                    table = pf.read()
                    met_cache[i] = table.to_pandas()

                for i in range(len(files)):
                    # ge_cache[i] = pd.read_parquet(f"{ge_dir_path}/{files[i]}")
                    pf = pq.ParquetFile(
                        f"{ge_dir_path}/{files[i]}",
                        thrift_string_size_limit=2**31-1,
                        thrift_container_size_limit=2**31-1,
                    )

                    table = pf.read()
                    ge_cache[i] = table.to_pandas()

                for i in tqdm(range(len(met_files)), desc="processing met", leave=False):
                    input1 = met_cache[i]
                    metadata1 = input1.iloc[:, :6].copy()  
                    mat1 = input1.iloc[:, 6:].copy()

                    min_pos_1 = metadata1['location_min'].min()
                    max_pos_1 = metadata1['location_max'].max()

                    for j in tqdm(range(len(files)), desc="running files", leave=False):
                        input2 = ge_cache[j]
                        metadata2 = input2.iloc[:, :6].copy()  
                        mat2 = input2.iloc[:, 6:].copy()

                        min_pos_2 = metadata2['location_min'].min()
                        max_pos_2 = metadata2['location_max'].max()

                        if max_pos_1 < min_pos_2 or max_pos_2 < min_pos_1:
                            continue

                        final_results = coloc_loop(
                            mat1=mat1,
                            mat2=mat2,
                            metadata1=metadata1,
                            metadata2=metadata2,
                            num_chunks1=math.ceil(mat1.shape[0]/100),
                            num_chunks2=math.ceil(mat2.shape[0]/100),
                            device=device,
                            p1=1e-4,
                            p2=1e-4,
                            p12=p12,
                            H4_threshold=H4_threshold,
                        )
                        
                        output_file=args.results

                        if final_results is None or final_results.empty:
                            continue

                        if not os.path.exists(output_file):
                            final_results.to_csv(output_file, sep="\t", index=False, mode='w', header=True)
                        else:
                            final_results.to_csv(output_file, sep="\t", index=False, mode='a', header=False)
            except Exception as e:
                print(f"Error reading files from {ge_dir_path}: {e}, possibly missing")
                continue

if __name__ == "__main__":
    main()