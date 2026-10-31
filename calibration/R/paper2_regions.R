# Paper 2: region classes of SNVs relative to MANE Select transcripts (GENCODE v50), with the rules of script 30
# (main article, amendment A10) extended to coding SNVs. Sourced by scripts 54 onwards.
# Noncoding priority when a SNV falls in several transcripts: canonical splice > splice region intronic > 5'UTR >
# 3'UTR > near-intronic > deep intronic > upstream 5 kb > downstream 5 kb; otherwise intergenic.
# Coding SNVs (MANE CDS positions): most severe consequence over MANE transcripts, stop gained > start lost >
# stop lost > missense > synonymous; synonymous split into exon edge (within 3 nt of an internal CDS junction,
# as in script 30) and other.
suppressPackageStartupMessages({ library(data.table) })
REGION_CLASSES <- c("canonical_splice", "splice_region_intronic", "5UTR", "3UTR", "near_intronic", "deep_intronic",
                    "upstream_5kb", "downstream_5kb", "intergenic", "missense", "synonymous_exon_edge", "synonymous_other",
                    "stop_gained", "start_lost", "stop_lost")
CDIR <- "C:/cadd/data/coding"

mane_structure <- function() {
  gtf <- fread("C:/cadd/annot/gencode.v50.basic.annotation.gtf.gz", sep = "\t", header = FALSE, skip = "chr1\t",
               col.names = c("chr", "src", "type", "start", "end", "score", "strand", "frame", "attr"))
  gtf <- gtf[grepl('tag "MANE_Select"', attr, fixed = TRUE) & type %chin% c("transcript", "exon", "CDS")]
  gtf[, `:=`(chrom = sub("^chr", "", chr), tx = sub('.*transcript_id "([^"]+)".*', "\\1", attr), gene = sub('.*gene_name "([^"]+)".*', "\\1", attr))]
  gtf <- gtf[!grepl("_PAR_Y$", tx)][chrom == "M", chrom := "MT"][, c("chr", "src", "score", "frame", "attr") := NULL]
  txs <- gtf[type == "transcript", .(chrom, tx, gene, strand, tstart = start, tend = end)]
  exons <- gtf[type == "exon", .(chrom, tx, start, end)]
  cds <- gtf[type == "CDS", .(cds_start = min(start), cds_end = max(end)), by = tx]
  introns <- exons[order(tx, start)][, .(start = end[-.N] + 1L, end = start[-1] - 1L), by = .(chrom, tx)][start <= end]
  list(txs = txs, exons = exons, cds = cds, introns = introns)
}

# q: data.table with chrom, pos (positions outside MANE CDS). Returns chrom, pos, gene, class (script 30 rules).
classify_noncoding <- function(q, m) {
  q <- unique(q[, .(chrom, pos)])[, `:=`(s = pos, e = pos)]; setkey(q, chrom, s, e)
  cls <- list()
  ex <- merge(merge(m$exons, m$cds, by = "tx"), m$txs[, .(tx, strand, gene)], by = "tx")
  setnames(ex, c("start", "end"), c("s", "e")); ex <- ex[!is.na(s) & !is.na(e)]; setkey(ex, chrom, s, e)
  o <- foverlaps(q, ex, nomatch = NULL)
  cls[[1]] <- o[, .(chrom, pos, gene, class = fifelse((strand == "+" & pos < cds_start) | (strand == "-" & pos > cds_end), "5UTR",
                                         fifelse((strand == "+" & pos > cds_end) | (strand == "-" & pos < cds_start), "3UTR", "exonic_other")))]
  it <- m$introns[m$txs[, .(tx, gene)], on = "tx", nomatch = NULL]; setnames(it, c("start", "end"), c("s", "e"))
  it <- it[!is.na(s) & !is.na(e)]; setkey(it, chrom, s, e)
  o <- foverlaps(q, it, nomatch = NULL)[, d := pmin(pos - s + 1L, e - pos + 1L)]
  cls[[2]] <- o[, .(chrom, pos, gene, class = fcase(d <= 2, "canonical_splice", d <= 8, "splice_region_intronic", d <= 50, "near_intronic", default = "deep_intronic"))]
  fl <- m$txs[, .(chrom, tx, gene, strand, s = pmax(1L, tstart - 5000L), e = tend + 5000L, tstart, tend)]; setkey(fl, chrom, s, e)
  o <- foverlaps(q, fl, nomatch = NULL)[pos < tstart | pos > tend]
  cls[[3]] <- o[, .(chrom, pos, gene, class = fifelse((strand == "+" & pos < tstart) | (strand == "-" & pos > tend), "upstream_5kb", "downstream_5kb"))]
  pri <- c("canonical_splice", "splice_region_intronic", "5UTR", "3UTR", "near_intronic", "deep_intronic", "upstream_5kb", "downstream_5kb", "exonic_other")
  reg <- rbindlist(cls)[, p := match(class, pri)][order(chrom, pos, p)][, .SD[1], by = .(chrom, pos)][, p := NULL]
  merge(q[, .(chrom, pos)], reg, by = c("chrom", "pos"), all.x = TRUE)[is.na(class), `:=`(class = "intergenic", gene = NA_character_)]
}

# Exon-edge coding positions (within 3 nt of an internal CDS junction), as in script 30.
exon_edge_positions <- function() {
  pos_all <- readRDS(file.path(CDIR, "mane_cds_positions.rds")); setorder(pos_all, transcript_id, cds_pos)
  pos_all[, jump := c(abs(diff(pos)) != 1L, FALSE), by = transcript_id][, idx := seq_len(.N), by = transcript_id]
  jn <- pos_all[jump == TRUE, .(transcript_id, kk = idx)]
  fl2 <- rbindlist(lapply(c(-2:0, 1:3), function(o) jn[, .(transcript_id, idx = kk + o)]))
  pos_all[, near := FALSE][fl2, on = .(transcript_id, idx), near := TRUE]
  unique(pos_all[near == TRUE, .(chrom, pos)])
}

# x: SNVs (chrom, pos, ref, alt) at MANE CDS positions of one chromosome. Returns the coding class per SNV.
classify_coding_chr <- function(x, edge) {
  f <- file.path(CDIR, sprintf("coding_snv_%s.rds", x$chrom[1])); if (!file.exists(f)) return(NULL)
  d <- readRDS(f)[, .(chrom, pos, ref, alt, gene_name, consequence)][x[, .(chrom, pos, ref, alt)], on = c("chrom", "pos", "ref", "alt"), nomatch = NULL]
  sev <- c("stop_gained", "start_lost", "stop_lost", "missense", "synonymous", "stop_retained")
  d[, r := match(consequence, sev)]
  setorder(d, chrom, pos, ref, alt, r); d <- unique(d, by = c("chrom", "pos", "ref", "alt"))
  d <- d[consequence != "stop_retained"]
  d[, class := consequence][class == "synonymous", class := "synonymous_other"]
  d[edge, on = .(chrom, pos), class := fifelse(class == "synonymous_other", "synonymous_exon_edge", class)]
  d[, .(chrom, pos, ref, alt, gene = gene_name, class)]
}
