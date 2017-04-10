# fit model to predict number of rare variants (MAF < 0.1%) for the BRIDGE data
# this script is for the 5,034 unrelated europeans

bridge_tsv=~/scratch/popgen/data/synonymous_all_chrs.EUR.tsv
genes=./data/gencode.v19.CDS.min_10_coverage.gene_union.txt
gene_mutation_rates=./data/coding_mutation_rates.gencode_v19_CDS.txt
bridge_model_out=~/scratch/popgen/synonymous/bridge_synonymous_lm.EUR.RData
bridge_pop_size=5034

cd ~/software/selection2/R
Rscript generate_synonymous_rare_var_lm.Rscript \
--variants $bridge_tsv --genes $genes \
--gene_mutation_rates $gene_mutation_rates \
--model_out $bridge_model_out --pop_size $bridge_pop_size

