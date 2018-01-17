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

option_list <- list(
  make_option("--elements", help = "Elements to annotate with constraint scores - expects 1-based coordinates."),
  make_option("--chromosome", help = "Chromosome to use (for running in parallel)"),
  make_option("--config_file", help = "File that describes the list of studies to use. Name of the study will be used to name columns."),
  make_option("--output", help = "where to save the elements annotated with obs/exp and z scores"),
  make_option("--coverage_correction", default=F, action='store_true', help = "use coverage correction models (in config file) to make coverage correction to observed variants"),
  make_option("--bed_file", default=F, action='store_true', help = "use this option if elements is bed file. assumes first three columns are chr, start, stop"),
  make_option("--rds_file", default=F, action='store_true', help = "use this option if elements is an RDS file - used in BRIDGE regulome project")
)

args <- parse_args(OptionParser(option_list=option_list))

source(args$config_file)

# load in regions file with required columns: chr, start, stop
if (args$bed_file) {
  elements <- read.table(args$elements, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  colnames(elements)[1:3] = c("chr", "start", "end")
} else if (args$rds_file) {
  elements <- readRDS(args$elements)
} else {
  if (summary( file(args$elements) )$class == "gzfile") {
    elements <- read.table(gzfile(args$elements), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  } else {
    elements <- read.delim(args$elements)
  }
}

# make all of the in the form chr1, chr2, etc.
if (any(!grepl("^chr", elements$chr))) {
  elements$chr = paste0("chr", elements$chr)
}

if (!(grepl("chr", args$chromosome))) {
  elements = subset(elements, chr %in% paste0("chr", args$chromosome))
} else {
  elements = subset(elements, chr %in% args$chromosome)
}


if (nrow(elements) == 0){
  stop("Did not find any elements after filtering by chromosome.")
}

# add region id
if (!("region_id" %in% colnames(elements))) {
  elements$region_id = paste0(elements$chr, ":", elements$start, "-", elements$end)
}

print("Elements passed:")
print(head(elements))

# get sequence context - some sequences will have N in them - we will have to give mutation rate of NA
if (!("seq" %in% colnames(elements))) {
  elements$seq = as.character(getSeq(Hsapiens, elements$chr, elements$start-1, elements$end+1))
}

# add phastcons100 score
library(phastCons100way.UCSC.hg19)
element_intervals = GRanges(seqnames=elements$chr, IRanges(start = elements$start, width = elements$end - elements$start + 1))
elements$phastcons100 = scores(phastCons100way.UCSC.hg19, element_intervals)

# add the base mutation rate
# for each element, need to pick out the sites that are CpGs and retrieve the prop methylated
# load in ESC WGBS data
wgbs = read.table(gzfile("~/scratch/CpG/data/E008_WGBS_FractionalMethylation.bedGraph.gz"), header = FALSE, sep = "\t")
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

# get probability per element with cpg adjustment
#load(methylation_correction_model.RData)
# old version, fit to BRIDGE data
#methyl_df = read.table("~/scratch/CpG/methylation_effect_on_obs_exp.ESC.txt", header = TRUE, sep = "\t")

#new version, fit to gnomAD data
methyl_df = read.table("/home/pjs90/scratch/CpG/methylation_effect_on_obs_exp.ESC_sperm_matched.txt", header = TRUE, sep = "\t")
methyl_df$prop_methylated = seq(0.025, 0.975, 0.05)
methylation_correction_model = lm(obs_exp_ratio ~ prop_methylated, methyl_df)

print(methyl_df)

seqs = elements$seq[match(names(prop_methylated), elements$region_id)]
seqs = split(seqs, f = names(prop_methylated))

# set the sequence mutability
elements$p_snp_null = sapply(elements$seq, p_sequence)

# save the triplet p_snp_null
elements$p_snp_null_no_methyl_correction = elements$p_snp_null

# correct based on CpGs methylation
elements$p_snp_null[match(names(prop_methylated), elements$region_id)] = mapply(p_sequence_meth, seqs, cpg_positions, prop_methylated, MoreArgs = list("correction_model" = methylation_correction_model))  # takes sequence, sites that are cpgs, prop methylated at those sites

print(head(elements))

elements$filter = "PASS"

for (study in studies) {
  print(study)
  starting_cols = colnames(elements)  

  # load in the null model to generate expected number
  load(study$null_model)  # loads as synonymous_lm

  # add expected per element
  elements$expected = predict(synonymous_lm, elements)

  # load in the variants to get actual observed
  vars = read.delim(sprintf(study$vars, args$chromosome))
  if (any(!grepl("^chr", vars$chr))) {
    vars$chr = paste0("chr", vars$chr)
  }
  
  vars = subset(vars, (nchar(as.character(vars$alt)) == 1) & (nchar(as.character(vars$ref)) == 1))
  
  print(head(vars))

  # add observed per element
  v = filter_with_bed(vars, elements)
  print('vars pre allele count filter')
  print(head(vars))
  mask = (v[,study$AC_column_name] < study$pop_size*2*0.001) & (v[,study$AC_column_name] > 0)
  v = subset(v, mask)
  print('vars post allele count filter')
  print(head(v))

  v = get_region_id_multi_overlap(v, elements)
  o = ddply(v, "region_id", function(df) data.frame(observed = sum(df$filter %in% study$pass_flags), observed_low_qual = sum(!(df$filter %in% study$pass_flags))))

  #elements = subset(elements, region_id %in% o$region_id)
  # get number of observed variants
  elements$observed = 0
  elements$observed[match(o$region_id, elements$region_id)] = o$observed
  
  # this will be used later to flag elements with many low quality variant calls
  elements$observed_low_qual = 0
  elements$observed_low_qual[match(o$region_id, elements$region_id)] = o$observed_low_qual

  # add observed/expected
  elements$obs_exp_ratio = elements$observed/elements$expected

  elements$low_qual_prop = elements$observed_low_qual/(elements$observed_low_qual + elements$observed)
  elements$filter[elements$low_qual_prop > 0.5 & !(elements$filter == "PASS")] = paste0(elements$filter, ";", "low_quality_0.5_", study$study_name)[elements$low_qual_prop > 0.5 & !(elements$filter == "PASS")]
  elements$filter[elements$low_qual_prop > 0.5 & !(elements$filter == "PASS")] = paste0("low_quality_0.5_", study$study_name)  

  # calculate Z score from observed and expected
  var_Z_score = function(observed, expected){
    Xsq_vals = (observed- expected)^2/expected
    excess = ifelse(observed > expected, 1, -1)
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

  # rename all of the columns using study name
  colnames(elements)[!(colnames(elements) %in% starting_cols)] = paste0(colnames(elements)[!(colnames(elements) %in% starting_cols)], "_", study$study_name)
}

print("Elements after adding variant counts and expected:")
print(head(elements))

# now, produce meta obs/exp
elements$meta_observed = rowSums(elements[,grepl("observed", colnames(elements)) & !grepl("low_qual", colnames(elements))])
elements$meta_expected = rowSums(elements[,grepl("expected", colnames(elements))])
elements$meta_obs_exp_ratio = elements$meta_observed/elements$meta_expected
elements$meta_z_score = var_Z_score(elements$meta_observed, elements$meta_expected)

if (args$coverage_correction) {
  load(studies$gnomad$gnomad_cov_correction)
  load(studies$BRIDGE$bridge_cov_correction)
  elements$expected_BRIDGE_cov_corrected = elements$expected_BRIDGE * predict(bridge_cov_correction_loess, elements)
  elements$obs_exp_ratio_BRIDGE_cov_corrected = elements$observed_BRIDGE/elements$expected_BRIDGE_cov_corrected

  elements$expected_gnomad_cov_corrected = elements$expected_gnomad * predict(gnomad_cov_correction_loess, elements)
  elements$obs_exp_ratio_gnomad_cov_corrected = elements$observed_gnomad/elements$expected_gnomad_cov_corrected

  elements$meta_expected_cov_corrected = (elements$expected_gnomad_cov_corrected + elements$expected_BRIDGE_cov_corrected)
  elements$meta_obs_exp_ratio_cov_corrected = (elements$observed_gnomad + elements$observed_BRIDGE)/elements$meta_expected_cov_corrected
  elements$meta_z_score_cov_corrected = sqrt(((elements$observed_gnomad + elements$observed_BRIDGE) - elements$meta_expected_cov_corrected)^2/elements$meta_expected_cov_corrected)

  for (study in studies) {
    cov = elements[,colnames(elements) %in% paste0("median_coverage_", study$study_name)]
    cov[is.na(cov)] = 0  # replace missing coverage with 0
    elements$filter[cov < 20 & !(elements$filter == "PASS")] = paste0(elements$filter, ";", "low_coverage_", study$study_name)[cov < 20 & !(elements$filter == "PASS")]
    elements$filter[cov < 20 & (elements$filter == "PASS")] = paste0("low_coverage_", study$study_name)
  }
}


# output annotated elements
write.table(elements, file = args$output, col.names = TRUE, sep = "\t", row.names = FALSE, quote = FALSE)
