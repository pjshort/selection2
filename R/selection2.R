# pass a tsv
# and the set of elements to calculate constraint in
# and variants from unaffected parents

# uses the dddMAPS library

library(stringr)
library(plyr)
library(optparse)
source("./annotation_tools.R")
source('./mutation_null_model.R')
library(BSgenome.Hsapiens.UCSC.hg19)
#source("~/software/dddMAPS/dddMAPS/MAPS.R")


### command line options
option_list <- list(
  make_option("--vars", help = "Pass list of variants in unaffected parents"),
  make_option("--elements", help = "Pass elements in which DNMs fall"),
  make_option("--null_model", default = "../data/ddd_synonymous_lm.RData"),
  make_option("--output", help = "location to save data frame with windows stats"),
  make_option("--pop_size", default=8082, help="Size of the population (number of individuals). Used to calculate allele count cutoff for rare vars (<0.1%).")
)

args <- parse_args(OptionParser(option_list=option_list))

# load in regions file with required columns: chr, start, stop
if (summary( file(args$elements) )$class == "gzfile") {
  elements <- read.table(gzfile(args$elements), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
} else {
  elements <- read.table(args$elements, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
}

colnames(elements)[1:3] = c("chr", "start", "stop")
if (any(!grepl("^chr", elements$chr))) {
  elements$chr = paste0("chr", elements$chr)
}

vars = read.table(args$vars, header = TRUE, sep = "\t")
vars$chr = paste0("chr", vars$chr)
vars = subset(vars, (nchar(as.character(vars$alt)) == 1) & (nchar(as.character(vars$ref)) == 1))
vars = subset(vars, filter == 'PASS')

print(args)
load(args$null_model)


# subset to relevant elements for running in array by chromosome
vars_chromosomes = unique(vars$chr)

elements = subset(elements, chr %in% vars_chromosomes)


# add region id
if (!("region_id" %in% colnames(elements))) {
  elements$region_id = paste0(elements$chr, ":", elements$start, "-", elements$stop)
}

# get sequence context
if (!("seq" %in% colnames(elements))) {
  elements$seq = as.character(getSeq(Hsapiens, elements$chr, elements$start, elements$stop))
}

elements = subset(elements, !grepl("N", seq))

print(head(elements))

# for each element, need to pick out the sites that are CpGs and retrieve the prop methylated
# load in ESC WGBS data
wgbs = read.table("~/scratch/CpG/data/E008_WGBS_FractionalMethylation.bedGraph", header = FALSE, sep = "\t")
colnames(wgbs) = c("chr", "start", "stop", "prop_methylated")
wgbs$pos = wgbs$start + 1
wgbs = wgbs[,c("chr", "pos", "prop_methylated")]

print(head(wgbs))

wgbs = subset(wgbs, chr %in% paste0("chr", seq(1,22)))
wgbs$chr = factor(wgbs$chr, levels = paste0("chr", seq(1,22)), ordered = TRUE)

wgbs = filter_with_bed(wgbs, elements)
wgbs = get_region_id_multi_overlap(wgbs, elements)

print(head(wgbs))

wgbs_split = split(wgbs, wgbs$region_id)
cpg_positions = sapply(wgbs_split, function(df) c(df$pos - elements$start[elements$region_id == df$region_id[1]] + 1))
prop_methylated = sapply(wgbs_split, function(df) c(df$prop_methylated))

seqs = elements$seq[match(names(prop_methylated), elements$region_id)]
seqs = split(seqs, f = names(prop_methylated))

# get probability per element with cpg adjustment
#load(methylation_correction_model.RData)
methyl_df = read.table("~/scratch/CpG/methylation_effect_on_obs_exp.ESC.txt", header = TRUE, sep = "\t")
methyl_df$prop_methylated = seq(0.025, 0.975, 0.05)
methylation_correction_model = lm(obs_exp_ratio ~ prop_methylated, methyl_df)

print(methyl_df)

elements = elements[match(names(prop_methylated), elements$region_id),]
elements$p_snp_null = 2 * mapply(p_sequence_meth, seqs, cpg_positions, prop_methylated, MoreArgs = list("correction_model" = methylation_correction_model))  # takes sequence, sites that are cpgs, prop methylated at those sites

# add expected per element
elements$expected = predict(synonymous_lm, elements)

print(head(elements))

# add observed per element
v = filter_with_bed(vars, elements)
print('vars pre allele count filter')
print(head(vars))
v = subset(v, allele_count < args$pop_size*2*0.001)
print('vars post allele count filter')
print(head(v))

v = get_region_id_multi_overlap(v, elements)
o = ddply(v, "region_id", function(df) data.frame(observed = sum(df$allele_count)))

print(head(v))

elements = subset(elements, region_id %in% o$region_id)
elements$observed = o$observed[match(o$region_id, elements$region_id)]

print(tail(elements[order(elements$observed),]))

# add observed/expected
elements$obs_exp_ratio = elements$observed/elements$expected



# calculate Z score from observed and expected
var_Z_score = function(observed, expected){
  Xsq_vals = (observed- expected)^2/expected
  excess = ifelse(observed > expected, -1, 1)
  Z = sqrt(Xsq_vals) * excess
  
  # use trimmed z scores to get standard deviation
  Z_trimmed = Z[ (Z > -5) & (Z < 5)]
  Z_trimmed = Z_trimmed[!is.na(Z_trimmed)]
  Z_sd = sd(Z_trimmed)
  
  # divide ALL Z scores by sd from middle set
  Z_normalized = Z/Z_sd
  
  return(Z_normalized)
}


# add Z score
elements$z_score = var_Z_score(elements$observed, elements$expected)


# add phastcons100 score
library(phastCons100way.UCSC.hg19)
element_intervals = GRanges(seqnames=elements$chr, IRanges(start = elements$start, width = elements$stop - elements$start + 1))
elements$phastcons100 = scores(phastCons100way.UCSC.hg19, element_intervals)


# output annotated elements
write.table(elements, file = args$output, col.names = TRUE, sep = "\t", row.names = FALSE, quote = FALSE)

