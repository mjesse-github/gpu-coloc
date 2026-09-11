import argparse
import math
import os
import torch
import pandas as pd
import numpy as np
from tqdm import tqdm
import pyarrow.parquet as pq
import time

SENTINEL = -1e6


def bf_to_pip(bf_df, device):
    """Convert a log-BF matrix to PIPs once, at load time.

    The "null" column is dropped before the softmax, matching what the old
    post-intersection code did. Each row sums to 1 over this file's SNPs;
    after intersecting with a second dataset rows sum to <= 1.
    """
    cols = [c for c in bf_df.columns if c != "null"]
    if not cols:
        return pd.DataFrame(index=bf_df.index)

    arr = torch.tensor(bf_df[cols].to_numpy(dtype=np.float32), device=device)
    arr = arr.masked_fill(arr <= SENTINEL / 2, float("-inf"))
    pip = torch.nan_to_num(torch.softmax(arr, dim=-1), nan=0.0)

    return pd.DataFrame(pip.cpu().numpy(), columns=cols, index=bf_df.index)


def CLPA(
    bf1_cpu, bf2_cpu,
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
            "clpa_matrix": None,
        }

    pip1 = torch.tensor(bf1_cpu[isnps].values, dtype=torch.float32, device=device)
    pip2 = torch.tensor(bf2_cpu[isnps].values, dtype=torch.float32, device=device)

    N, _ = pip1.shape
    K, _ = pip2.shape

    # CLPA[i, j] = sum_s min(pip1[i, s], pip2[j, s])
    clpa_2d = torch.empty(N, K, dtype=torch.float32, device=device)
    for i in range(N):
        clpa_2d[i] = torch.minimum(pip1[i], pip2).sum(dim=-1)

    # CLPP[i, j] = sum_s pip1[i, s] * pip2[j, s]
    clpp_2d = pip1 @ pip2.T

    i_coords = torch.arange(N, device=device).unsqueeze(1).expand(N, K).flatten()
    j_coords = torch.arange(K, device=device).unsqueeze(0).expand(N, K).flatten()

    summary_df = pd.DataFrame({
        "idx1": i_coords.cpu().numpy(),
        "idx2": j_coords.cpu().numpy(),
        "CLPA": clpa_2d.flatten().cpu().numpy(),
        "CLPP": clpp_2d.flatten().cpu().numpy(),
    })

    return {
        "summary": summary_df,
    }


def CLPA_loop(
    mat1: pd.DataFrame,
    mat2: pd.DataFrame,
    metadata1: pd.DataFrame,
    metadata2: pd.DataFrame,
    n_tests,
    chunk_size=100,
    num_chunks1=0,
    num_chunks2=0,
    device="cuda",
    coloc_time=0
):

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

        start = time.time()

        out = CLPA(
            bf1_cpu=mat1_chunks[pair[0]],
            bf2_cpu=mat2_chunks[pair[1]],
            device=device
        )

        end = time.time()
        coloc_time += end - start

        if out is None or out["summary"] is None:
            continue

        summary_df = out["summary"]

        if {"idx1", "idx2"} - set(summary_df.columns):
            continue

        summary_df.loc[:, "idx1"] = summary_df["idx1"] + pair[0] * chunk_size
        summary_df.loc[:, "idx2"] = summary_df["idx2"] + pair[1] * chunk_size

        n_tests += summary_df.shape[0]

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

    return final_df, n_tests, coloc_time


def main():
    run_start = time.time()

    n_tests = 0

    parser = argparse.ArgumentParser(description="Run coloc")

    parser.add_argument("--dir1", type=str, required=True, help="First directory of directories of parquet files, e.g., 'formatted_eqtls'.")
    parser.add_argument("--dir2", type=str, required=True, help="Second directory of directories of parquet files, e.g., 'formatted_metabolites'.")
    parser.add_argument("--results", type=str, required=True, help="File to write the colocalization results, e.g., 'results.tsv'.")
    parser.add_argument("--verbose", action="store_true", help="Print timing and test info")
    parser.add_argument("--CPU", action="store_true", help="Force Torch calculations on CPU")
    parser.add_argument("--chunk_size", type=int, required=False, help="number of signals in a chunk", default=1000)

    args = parser.parse_args()

    IO_time = 0
    coloc_time = 0

    if args.CPU:
        device = torch.device("cpu")
    else:
        if torch.cuda.is_available():
            device = torch.device("cuda")
        elif torch.backends.mps.is_available():
            device = torch.device("mps")
        else:
            device = torch.device("cpu")

    chunk_size = args.chunk_size

    if args.verbose:
        print(f"using {device} backend")
        print(f"threads: {torch.get_num_threads()}")
        print(f"signals in a chunk: {chunk_size}")

    for root, dirs, _ in os.walk(args.dir1):
        for directory in tqdm(dirs, desc="chromosomes"):
            try:
                dir1_path = os.path.join(root, directory)
                dir1_files = os.listdir(dir1_path)

                dir2_path = os.path.join(args.dir2, directory)
                dir2_files = os.listdir(dir2_path)

                dir1_cache = {}

                start = time.time()

                for i in range(len(dir1_files)):
                    pf = pq.ParquetFile(
                        os.path.join(dir1_path, dir1_files[i]),
                        thrift_string_size_limit=2**31-1,
                        thrift_container_size_limit=2**31-1,
                    )

                    table = pf.read().to_pandas()

                    # PIP conversion happens once per file, not once per pair
                    dir1_cache[i] = (
                        table.iloc[:, :6].copy(),
                        bf_to_pip(table.iloc[:, 6:], device),
                    )

                    del table

                end = time.time()
                IO_time += end - start

                for j in tqdm(range(len(dir2_files)), desc="processing outer files", leave=False):
                    start = time.time()

                    pf = pq.ParquetFile(
                        os.path.join(dir2_path, dir2_files[j]),
                        thrift_string_size_limit=2**31-1,
                        thrift_container_size_limit=2**31-1,
                    )

                    table = pf.read().to_pandas()

                    end = time.time()
                    IO_time += end - start

                    metadata2 = table.iloc[:, :6].copy()
                    pip2 = bf_to_pip(table.iloc[:, 6:], device)

                    del table

                    min_pos_2 = metadata2['location_min'].min()
                    max_pos_2 = metadata2['location_max'].max()

                    for i in tqdm(range(len(dir1_files)), desc="processing inner files", leave=False):
                        metadata1, pip1 = dir1_cache[i]

                        min_pos_1 = metadata1['location_min'].min()
                        max_pos_1 = metadata1['location_max'].max()

                        if max_pos_1 < min_pos_2 or max_pos_2 < min_pos_1:
                            continue

                        final_results, n_tests, coloc_time = CLPA_loop(
                            mat1=pip1,
                            mat2=pip2,
                            metadata1=metadata1,
                            metadata2=metadata2,
                            n_tests=n_tests,
                            chunk_size=chunk_size,
                            num_chunks1=math.ceil(pip1.shape[0]/chunk_size),
                            num_chunks2=math.ceil(pip2.shape[0]/chunk_size),
                            device=device,
                            coloc_time=coloc_time
                        )

                        output_file = args.results

                        if final_results is None or final_results.empty:
                            continue

                        if not os.path.exists(output_file):
                            final_results.to_csv(output_file, sep="\t", index=False, mode='w', header=True)
                        else:
                            final_results.to_csv(output_file, sep="\t", index=False, mode='a', header=False)
            except Exception as e:
                print(f"Error while using files from {dir2_path}: {e}")
                continue

        dirs.clear()

    if args.verbose:
        print(f"{n_tests} pairs tested for colocalisation")
        print(f"IO time: {IO_time}")
        print(f"coloc time: {coloc_time}")
        print(f"Total elapsed time: {time.time() - run_start}")


if __name__ == "__main__":
    main()