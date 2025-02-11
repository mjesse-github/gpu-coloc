import os
import numpy as np
import pandas as pd
from tqdm import tqdm

def approx_bf_estimates(variant, z, V, type, suffix=None, sdY=1, effect_priors={'quant': 0.15, 'cc': 0.2}):
    sd_prior = effect_priors['quant'] * sdY if type == "quant" else effect_priors['cc']
    
    r = sd_prior**2 / (sd_prior**2 + V)
    lABF = 0.5 * (np.log(1 - r) + (r * z**2))
    
    ret = pd.DataFrame({"variant": variant, 'lbf': lABF})
    
    if suffix is not None:
        ret.columns = [f"{col}.{suffix}" for col in ret.columns]

    ret.loc[:, 'lbf'] = pd.to_numeric(ret['lbf'], errors='coerce')

    ret_mat = pd.DataFrame(ret.set_index('variant')['lbf']).T

    return ret_mat

summary_file_path = "met_metadata.tsv"

signals_dir = "met_mat"
os.makedirs(signals_dir, exist_ok=True)

lead_variants = pd.read_csv("meta_EUR_all_lead_variants.tsv", sep="\t")

for metabolite, group in lead_variants.groupby("metabolite"):
    regions = group

    metabolite_df = pd.read_parquet(f"/gpfs/space/projects/alasoo/metabolite_gwas/EstBB_UKB_meta-analysis/metaanalysis/Est_vs_EUR_meta/{metabolite}_EstBB_UKBB_EUR_metaanalysis.parquet")

    for _, region in tqdm(regions.iterrows(), desc="processing regions"):
        chrom = region["CHR"] 
        region_start = region["POS"] - 1_000_000
        region_end = region["POS"] + 1_000_000
        region_met = metabolite_df[metabolite_df["CHROM"]==chrom]

        if chrom == 23:
            chrom = "X"

        region_signal = region_met[(region_met["GENPOS"]>=region_start) & (region_met["GENPOS"]<=region_end)]

        region_signal["MAF"] = np.minimum(region_signal["AC"] / (2 * region_signal["N"]), 1 - region_signal["AC"] / (2 * region_signal["N"]))

        region_signal = region_signal[region_signal["MAF"]>=0.01]

        if region_signal.empty: 
            continue

        region_signal["variant"] = (
        "chr" + region_signal["CHROM"].astype(int).astype(str) + "_" +
        region_signal["GENPOS"].astype(int).astype(str) + "_" +
        region_signal["ALLELE0"].astype(str) + "_" +
        region_signal["ALLELE1"].astype(str)
        )

        region_signal["z"]=region_signal["Effect"] / region_signal["StdErr"]
        region_signal["V"]= region_signal["StdErr"]**2

        mat = approx_bf_estimates(region_signal["variant"],region_signal["z"], region_signal["V"], "quant") 
        
        signal_strength = mat.T["lbf"].max()
        variant_id = mat.T["lbf"].idxmax()

        if signal_strength < 5:
            continue

        signal = f"{metabolite}_chr{chrom}:{region_start}-{region_end}"

        summary_data = pd.DataFrame([{
            'signal': signal,
            'chromosome': chrom,
            'location_min': region_start,
            'location_max': region_end,
            'signal_strength': signal_strength,
            'lead_variant': variant_id
        }])

        header_needed = not os.path.exists(summary_file_path)  
        summary_data.to_csv(summary_file_path, sep='\t', mode='a', header=header_needed, index=False)

        output_file_name = f"{signal}.pickle"
        output_file_path = os.path.join(signals_dir, output_file_name)

        mat.to_pickle(output_file_path)

