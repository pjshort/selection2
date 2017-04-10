# selection2

# Requirements and Data
R version 3.2.2 or higher is recommended, and only the GenomicRanges packages is required.

`install.packages("GenomicRanges")`

A file with tri-nucleotide mutation rates from Samocha et. al, 2014 and whole genome bisulfite sequencing data from embryonic stem cells and sperm is in the `/data` folder. 

The ES-cell data was downloaded from http://egg2.wustl.edu/roadmap/data/byDataType/dnamethylation/WGBS/FractionalMethylation_bigwig/ and the UCSC software toolkit was used to transform the bigWig file into a bedGraph and then into a bed file. The sperm methylome data was downloaded from "Genomic Hypomethylation in the Human Germline Associates with Selective Structural Mutability in the Human Genome
Li, J., et al. 2012" and lifted over from hg18 to hg19 using the UCSC liftover tool.

Any set of elements and variants in hg19 coordinates can be used. The elements are expected to have at least columns named `chr, start, stop` and the variants should have columns `chr, pos, ref, alt, allele_count` where allele count is the number of times the variant was observed in the population.
