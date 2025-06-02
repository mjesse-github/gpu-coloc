# gpu-coloc

**gpu-coloc** is a GPU-accelerated implementation of the Bayesian colocalization algorithm (COLOC), providing identical results to R's coloc.bf\_bf at approximately 1000x greater speed.

## Citation

If you use **gpu-coloc**, please cite: *(citation placeholder)*

## Installation

Clone the repository:

```bash
git clone https://github.com/mjesse-github/gpu-coloc
```

### Dependencies

Install required Python libraries:

```bash
pip install -r requirements.txt
```

For Linux x64 servers, we recommend using our Singularity container:
*(Singularity link placeholder)*

## Testing Installation

Run:

```bash
bash test.sh
```

## Workflow

Note: The following example assumes gpu-coloc is downloaded into your working directory. Adjust paths accordingly if downloaded elsewhere.

Variants must follow a uniform naming convention, as the COLOC algorithm requires consistent naming. Use the format: chr[chromosome]_[position]_[ref]_[alt]. Perform any renaming prior to Step 1 below. We use chromosome X, not 23.

1. **Prepare Signals and Summary Files**

   * **Signals files**: Each signal should be saved in `[signal].pickle` format, containing variants and their respective log Bayes Factors (lbf).

Example format:

```
variant	chrX_153412224_C_A	chrX_153412528_C_T	...
lbf	-0.060991	-1.508802	...
```

* **Summary file**: Tab-separated file with the structure below:

```
signal	chromosome	location_min	location_max	signal_strength	lead_variant
QTD000141_ENSG00000013563_L1	X	153412224	155341332	12.1069377174147	chrX_154403855_T_G
...
```

Example naming convention:

* `gwas_summary.tsv`
* Signals in directory `gwas_signals/[signal].pickle`

Scripts `gwas_signals.py` and `eqtl_signals.py` are provided as examples, but may require adjustments.

2. **Formatting data for colocalization:**

```bash
python3 gpu-coloc/format.py --input [path_to_signals] --input_summary [summary_file] --output [output_folder]
```

The formatting algorithm creates a seperate folder for each chromosome into the .parquet file output folder, so it is important to keep in mind that both .parquet folders have the same chromosome naming convention (for chromosome X we use X not 23). Also variants should hold a uniform naming convention, as COLOC algorithm is dependent on this, we use the convention chr[chromosme]_[position]_[ref]_[alt], however renaming should be done before the first step.

3. **Run colocalization analysis:**

```bash
python3 gpu-coloc/coloc.py --dir1 [formatted_dataset_1] --dir2 [formatted_dataset_2] --results [results_output] --p12 1e-6 --H4 0.8
```
