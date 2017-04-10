library(GenomicRanges)

filter_with_bed <- function(de_novos, bed){
  # assumes that the bed input has at least three columns with chr, start, end as the first three (all others ignored)

  if (any(!grepl("^chr", de_novos$chr))) {
    de_novos$chr = paste0("chr", de_novos$chr)
  }
  if (any(!grepl("^chr", bed[,1]))) {
    bed[ ,1] = paste0("chr", bed[ ,1])
  }

  dn = GRanges(seqnames=Rle(de_novos$chr), ranges = IRanges(start = de_novos$pos, end = de_novos$pos))
  gen = GRanges(seqnames=Rle(bed[,1]), ranges = IRanges(start = bed[,2], end = bed[,3]))

  # find overlap between denovos and genomicus annotation
  hits = findOverlaps(dn, gen)
  dn_hits_idx = unique(queryHits(hits)) # get index of de novos with genomicus hit

  de_novos = de_novos[dn_hits_idx, ]

  return(de_novos)
}

exclude_with_bed <- function(de_novos, bed){
  # assumes that the bed input has at least three columns with chr, start, end as the first three (all others ignored)

  if (any(!grepl("^chr", de_novos$chr))) {
    de_novos$chr = paste0("chr", de_novos$chr)
  }
  if (any(!grepl("^chr", bed[,1]))) {
    bed[ ,1] = paste0("chr", bed[ ,1])
  }

  dn = GRanges(seqnames=Rle(de_novos$chr), ranges = IRanges(start = de_novos$pos, end = de_novos$pos))
  gen = GRanges(seqnames=Rle(bed[,1]), ranges = IRanges(start = bed[,2], end = bed[,3]))

  # find overlap between denovos and genomicus annotation
  hits = findOverlaps(dn, gen)
  dn_hits_idx = unique(queryHits(hits)) # get index of de novos with genomicus hit

  de_novos = de_novos[-dn_hits_idx, ]

  return(de_novos)
}



get_region_id <- function(de_novos, CNEs){
  # assumes that the CNE input has at least three columns with chr, start, end as the first three (all others ignored)

  if (any(!grepl("^chr", de_novos$chr))) {
    de_novos$chr = paste0("chr", de_novos$chr)
  }
  if (any(!grepl("^chr", CNEs$chr))) {
    CNEs$chr = paste0("chr", CNEs$chr)
  }

  if ("end" %in% colnames(de_novos)){ # region instead of de novos - use the first position of the region to get closest gene!
    dn = GRanges(seqnames=Rle(de_novos$chr), ranges = IRanges(start = de_novos$start, end = de_novos$start))
  } else {
    dn = GRanges(seqnames=Rle(de_novos$chr), ranges = IRanges(start = de_novos$pos, end = de_novos$pos))
  }

  cne = GRanges(seqnames=Rle(CNEs$chr), ranges = IRanges(start = CNEs$start, end = CNEs$stop))

  # find overlap between denovos and annotated CNEs
  hits = findOverlaps(dn, cne)
  dn_hits_idx = queryHits(hits) # get index of de novos
  CNE_hits_idx = subjectHits(hits) # get index of CNEs


  return(CNEs$region_id[CNE_hits_idx])
}

get_region_id_multi_overlap <- function(vars, elements){
  # get the region ID in cases where the elements used might be overlapping (for tiling selection test)

  if (any(!grepl("^chr", vars$chr))) {
    vars$chr = paste0("chr", vars$chr)
  }
  if (any(!grepl("^chr", elements$chr))) {
    elements$chr = paste0("chr", elements$chr)
  }

  if ("end" %in% colnames(vars)){ # region instead of de novos - use the first position of the region to get closest gene!
    v = GRanges(seqnames=Rle(vars$chr), ranges = IRanges(start = vars$start, end = vars$start))
  } else {
    v = GRanges(seqnames=Rle(vars$chr), ranges = IRanges(start = vars$pos, end = vars$pos))
  }

  e = GRanges(seqnames=Rle(elements$chr), ranges = IRanges(start = elements$start, end = elements$stop))

  # find overlap between denovos and annotated elements
  hits = findOverlaps(v, e)
  v_hits_idx = queryHits(hits) # get index of de novos
  e_hits_idx = subjectHits(hits) # get index of elements

  vars = vars[v_hits_idx,]  # will repeat rows of v if there are overlapping elements
  vars$region_id = elements$region_id[e_hits_idx]

  return(vars)
}
