# Builds the coding-sequence coordinates of every GENCODE v50 MANE Select transcript:
# one row per coding nucleotide (CDS + stop codon), in transcript order, with its position in
# the CDS. Output: C:/cadd/data/coding/mane_cds_positions.rds (+ transcript summary .tsv)
# Usage: Rscript 08_build_mane_cds.R
suppressPackageStartupMessages(library(data.table))
gtf_file <- "C:/cadd/annot/gencode.v50.basic.annotation.gtf.gz"
outdir   <- "C:/cadd/data/coding"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

gtf <- fread(gtf_file, sep = "\t", header = FALSE, skip = "chr1\t",
             col.names = c("chr", "src", "type", "start", "end", "score", "strand", "frame", "attr"))
gtf <- gtf[type %chin% c("CDS", "stop_codon") & grepl('tag "MANE_Select"', attr, fixed = TRUE)]
att <- function(k) sub(paste0('.*', k, ' "([^"]+)".*'), "\\1", gtf$attr)
gtf[, `:=`(gene_name = att("gene_name"), gene_id = att("gene_id"),
           transcript_id = att("transcript_id"), gene_type = att("gene_type"))]
gtf[, attr := NULL]
# PAR copies on chrY carry the suffix _PAR_Y; keep the chrX copy only
gtf <- gtf[!grepl("_PAR_Y$", transcript_id)]
cat("MANE Select transcripts with CDS:", uniqueN(gtf$transcript_id), "| genes:", uniqueN(gtf$gene_id), "\n")

# Expand each feature to single nucleotides and order them 5'->3' along the transcript
pos <- gtf[, .(pos = seq.int(start, end)), by = .(chr, strand, gene_name, gene_id, gene_type, transcript_id, type, start)]
setorder(pos, transcript_id, pos)
pos[strand == "-", ord := rev(seq_len(.N)), by = transcript_id]
pos[strand == "+", ord := seq_len(.N), by = transcript_id]
setorder(pos, transcript_id, ord)
pos[, cds_pos := ord][, c("ord", "start") := NULL]
pos[, chrom := sub("^chr", "", chr)][chrom == "M", chrom := "MT"][, chr := NULL]
setcolorder(pos, c("chrom", "pos", "strand", "gene_name", "gene_id", "gene_type", "transcript_id", "cds_pos", "type"))

summ <- pos[, .(chrom = chrom[1], strand = strand[1], gene_name = gene_name[1], gene_type = gene_type[1],
                cds_length = .N, has_stop_feature = any(type == "stop_codon")), by = transcript_id]
cat("coding nucleotides:", nrow(pos), "| unique genomic positions:", uniqueN(pos, by = c("chrom", "pos")), "\n")
cat("CDS length multiple of 3:", summ[cds_length %% 3 == 0, .N], "of", nrow(summ), "\n")
cat("transcripts with stop_codon feature:", summ[has_stop_feature == TRUE, .N], "\n")
print(summ[, .N, by = gene_type][order(-N)])

saveRDS(pos, file.path(outdir, "mane_cds_positions.rds"))
fwrite(summ, file.path(outdir, "mane_transcripts_summary.tsv"), sep = "\t")
cat("Written:", file.path(outdir, "mane_cds_positions.rds"), "\n")
