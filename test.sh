python3 gwas_signals.py --summary example/gwas_summary.tsv --output example/gwas_signals --input_table example_data/example_gwas_table.tsv
python3 eqtl_signals.py --summary example/eqtl_summary.tsv --output example/eqtl_signals --input_table example_data/example_eQTL_table.tsv
python3 format.py --input example/gwas_signals --output example/formatted_gwas --input_summary example/gwas_summary.tsv --output_summary example/gwas_files_summary.tsv
python3 format.py --input example/eqtl_signals --output example/formatted_eqtls --input_summary example/eqtl_summary.tsv --output_summary example/eqtl_files_summary.tsv
python3 coloc.py --dir1 example/formatted_eqtls --dir2 example/formatted_gwas --results example/example_results.tsv