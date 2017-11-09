# selection2

# Requirements and Data
R version 3.2.2 or higher is recommended. Install the following packages:

`install.packages("GenomicRanges")
install.packages("plyr")`
install.packages("optparse")`

`source("https://bioconductor.org/biocLite.R")
biocLite("BSgenome.Hsapiens.UCSC.hg19")`

A file with tri-nucleotide mutation rates from Samocha et. al, 2014 and whole genome bisulfite sequencing data from embryonic stem cells and sperm is in the `/data` folder. 

The ES-cell data was downloaded from http://egg2.wustl.edu/roadmap/data/byDataType/dnamethylation/WGBS/FractionalMethylation_bigwig/ and the UCSC software toolkit was used to transform the bigWig file into a bedGraph and then into a bed file. The sperm methylome data was downloaded from "Genomic Hypomethylation in the Human Germline Associates with Selective Structural Mutability in the Human Genome
Li, J., et al. 2012" and lifted over from hg18 to hg19 using the UCSC liftover tool.

Any set of elements and variants in hg19 coordinates can be used. The elements are expected to have at least columns named `chr, start, stop` and the variants should have columns `chr, pos, ref, alt, allele_count, CQ` where allele count is the number of times the variant was observed in the population and CQ is the VEP consequence (most damaging of all transcripts).

# Running the scripts

First, the synonymous variants should be used to generate a model to predict the number of rare variants with MAF < 0.1% for a given sequence based on its mutability.

The following slurm jobs will generate these files for BRIDGE, gnomad (requires some pre-processing to get the synonymous variants from these studies):

`scripts/fit_synonymous_models_BRIDGE_slurm
scripts/fit_synonymous_models_gnomAD_slurm`

Next, these models can be used to run `selection2.scalable.R` which uses a config file designed to make it easier to scale the number of deep WGS used to meta-analyze.

`
cd ./R
Rscript selection2.scalable.R \
--config_file /path/to/config
--elements /path/to/elements \
--chromosome <number>
--output /path/to/output \
`

An example script to run on >1 million elements can be found in `scripts/selection2_FLAGSHIP`
