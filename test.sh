python3 summary_and_signals_examples/gwas_signals.py --summary example/gwas_summary.tsv --output example/gwas_signals --input_table example_data/example_gwas_table.tsv
python3 summary_and_signals_examples/eqtl_signals.py --summary example/eqtl_summary.tsv --output example/eqtl_signals --input_table example_data/example_eQTL_table.tsv
gpu-coloc -f --input example/gwas_signals --output example/formatted_gwas --input_summary example/gwas_summary.tsv --output_summary example/gwas_files_summary.tsv
gpu-coloc -f --input example/eqtl_signals --output example/formatted_eqtls --input_summary example/eqtl_summary.tsv --output_summary example/eqtl_files_summary.tsv
gpu-coloc -r --dir1 example/formatted_eqtls --dir2 example/formatted_gwas --results example/example_results.tsv --p12 1e-6 --H4 0.8