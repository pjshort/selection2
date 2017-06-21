# selection2

# Requirements and Data
R version 3.2.2 or higher is recommended, and only the GenomicRanges packages is required.

`install.packages("GenomicRanges")`

A file with tri-nucleotide mutation rates from Samocha et. al, 2014 and whole genome bisulfite sequencing data from embryonic stem cells and sperm is in the `/data` folder. 

The ES-cell data was downloaded from http://egg2.wustl.edu/roadmap/data/byDataType/dnamethylation/WGBS/FractionalMethylation_bigwig/ and the UCSC software toolkit was used to transform the bigWig file into a bedGraph and then into a bed file. The sperm methylome data was downloaded from "Genomic Hypomethylation in the Human Germline Associates with Selective Structural Mutability in the Human Genome
Li, J., et al. 2012" and lifted over from hg18 to hg19 using the UCSC liftover tool.

Any set of elements and variants in hg19 coordinates can be used. The elements are expected to have at least columns named `chr, start, stop` and the variants should have columns `chr, pos, ref, alt, allele_count, CQ` where allele count is the number of times the variant was observed in the population and CQ is the VEP consequence (most damaging of all transcripts).

# Running the scripts

First, the synonymous variants should be used to generate a model to predict the number of rare variants with MAF < 0.1% for a given sequence based on its mutability.

The example provided, `synonymous_model_BRIDGE.sh` uses synonymous sites from deep whole genomes from 5,034 individuals in the BRIDGE cohort.

`Rscript generate_synonymous_rare_var_lm.Rscript \
--variants /path/to/variants --genes /path/to/genes/or/exons \
--gene_mutation_rates /path/to/mutation/rate/per/gene \
--model_out /path/to/output --pop_size pop_size`

Next, this model can be used to run the `selection2.R` script which takes the following command line arguments:

`Rscript selection2.R \
--vars /path/to/variants \
--elements /path/to/elements \
--null_model /path/to/model/from/above \
--output /path/to/output \
--pop_size pop_size
`

This script was written to be easily parallelized, for instance in slurm by running 22 separate jobs (one for each chromosome) with a different set of input variants. The set of elements can be kept the same and the script will only generate a score for those for which there are variants included in the input. For instance, if you only pass the chromosome 1 variants, obs/exp ratio and Z-score will only be generated for elements on chromosome 1.
