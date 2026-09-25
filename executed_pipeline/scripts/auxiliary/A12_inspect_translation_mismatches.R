# Inspects the MANE Select transcripts whose translated reference CDS differs from the GENCODE v50
# protein (qc_translation.tsv from 10_annotate_consequences.R): gene, each differing residue
# (ours vs GENCODE) and the genomic codon, to classify the causes (e.g. selenocysteine UGA).
suppressPackageStartupMessages({ library(data.table); library(Biostrings) })
dir <- "C:/cadd/data/coding"
qc  <- fread(file.path(dir, "qc_translation.tsv"))[identical == FALSE | is.na(identical)]
pos <- readRDS(file.path(dir, "mane_cds_positions.rds"))[transcript_id %in% qc$transcript_id]
prot <- readAAStringSet("C:/cadd/annot/gencode.v50.pc_translations.fa.gz")
names(prot) <- vapply(strsplit(names(prot), "|", fixed = TRUE), `[`, "", 2)
out <- list()
for (tx in qc$transcript_id) {
  ch <- qc[transcript_id == tx, chrom]
  snv <- readRDS(file.path(dir, sprintf("coding_snv_%s.rds", ch)))[transcript_id == tx]
  cod <- unique(snv[, .(protein_pos, aa_ref, ref_codon)])[order(protein_pos)]
  g <- as.character(prot[[tx]]); gv <- strsplit(g, "")[[1]]
  ours <- cod[aa_ref != "*" | protein_pos < max(protein_pos)]
  gene <- pos[transcript_id == tx, gene_name[1]]
  if (nrow(ours) != length(gv)) {
    out[[tx]] <- data.table(transcript_id = tx, gene = gene, protein_pos = NA, ours = NA, gencode = NA,
                            codon = NA, note = sprintf("length differs: ours %d vs GENCODE %d; CDS %d nt",
                                                       nrow(ours), length(gv), qc[transcript_id == tx, cds_nt]))
    next
  }
  d <- which(ours$aa_ref != gv)
  out[[tx]] <- data.table(transcript_id = tx, gene = gene, protein_pos = ours$protein_pos[d],
                          ours = ours$aa_ref[d], gencode = gv[d], codon = ours$ref_codon[d], note = "")
}
res <- rbindlist(out)
res[, cause := fcase(ours == "*" & gencode == "U" & codon == "TGA", "selenocysteine (UGA recoded)",
                     !is.na(ours) & ours != "*", "reference codon differs from GENCODE protein",
                     default = "other")]
print(res[, .N, by = cause])
print(res, nrows = 200)
fwrite(res, file.path(dir, "qc_translation_mismatches.tsv"), sep = "\t")
