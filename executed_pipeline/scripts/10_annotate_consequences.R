# Consequence of every possible coding SNV on its MANE Select transcript.
# Reference codons are rebuilt from the reference alleles of the CADD file itself (so coordinates
# and alleles are consistent with the scores), translated with the standard genetic code (the
# vertebrate mitochondrial code for chrMT). QC: every reference CDS is translated and compared with
# the GENCODE v50 protein sequence of the same transcript.
# Input : mane_cds_positions.rds, cadd_coding_<chrom>.rds (08, 09)
# Output: coding_snv_<chrom>.rds (one row per SNV x MANE transcript), qc_translation.tsv
# Usage: Rscript 10_annotate_consequences.R
suppressPackageStartupMessages({ library(data.table); library(Biostrings) })
dir <- "C:/cadd/data/coding"
pos_all <- readRDS(file.path(dir, "mane_cds_positions.rds"))

comp <- c(A = "T", C = "G", G = "C", T = "A")
code_std <- GENETIC_CODE
code_mt  <- getGeneticCode("SGC1")          # vertebrate mitochondrial
aa3 <- c(A="Ala",R="Arg",N="Asn",D="Asp",C="Cys",Q="Gln",E="Glu",G="Gly",H="His",I="Ile",L="Leu",K="Lys",
         M="Met",F="Phe",P="Pro",S="Ser",T="Thr",W="Trp",Y="Tyr",V="Val","*"="Ter")

# Reference proteins from GENCODE (header: >ENSP|ENST|ENSG|...), keyed by transcript ID
prot <- readAAStringSet("C:/cadd/annot/gencode.v50.pc_translations.fa.gz")
names(prot) <- vapply(strsplit(names(prot), "|", fixed = TRUE), `[`, "", 2)

qc_all <- list()
for (ch in unique(pos_all$chrom)) {
  f <- file.path(dir, sprintf("cadd_coding_%s.rds", ch))
  if (!file.exists(f)) { cat("no CADD file for chromosome", ch, "\n"); next }
  cadd <- readRDS(f)
  pos  <- pos_all[chrom == ch]
  code <- if (ch == "MT") code_mt else code_std

  # Reference base on the coding strand, per transcript position
  refpos <- unique(cadd[, .(chrom, pos, ref)])
  pos <- merge(pos, refpos, by = c("chrom", "pos"), all.x = TRUE)
  pos[, cbase := ifelse(strand == "+", ref, comp[ref])]
  setorder(pos, transcript_id, cds_pos)
  pos[, `:=`(codon = (cds_pos - 1L) %/% 3L + 1L, frame = (cds_pos - 1L) %% 3L + 1L)]
  cod <- pos[, .(ref_codon = paste(cbase, collapse = "")), by = .(transcript_id, codon)]

  # QC: translate the full reference CDS (without the final stop) and compare with GENCODE
  qc <- cod[, {
    seq <- paste(ref_codon, collapse = "")
    ok_len <- nchar(seq) %% 3 == 0 && !grepl("[^ACGT]", seq)
    tr <- if (ok_len) as.character(translate(DNAString(seq), genetic.code = code, no.init.codon = FALSE,
                                             if.fuzzy.codon = "X")) else NA_character_
    tr_nostop <- sub("\\*$", "", tr)
    ref_p <- if (transcript_id %in% names(prot)) as.character(prot[[transcript_id]]) else NA_character_
    list(cds_nt = nchar(seq), translated_len = nchar(tr_nostop), gencode_len = nchar(ref_p),
         identical = !is.na(tr_nostop) && !is.na(ref_p) && tr_nostop == ref_p,
         n_diff = if (!is.na(tr_nostop) && !is.na(ref_p) && nchar(tr_nostop) == nchar(ref_p))
                    sum(strsplit(tr_nostop, "")[[1]] != strsplit(ref_p, "")[[1]]) else NA_integer_,
         internal_stops = if (is.na(tr_nostop)) NA_integer_ else lengths(regmatches(tr_nostop, gregexpr("\\*", tr_nostop))))
  }, by = transcript_id]
  qc[, chrom := ch]
  qc_all[[ch]] <- qc

  # Every SNV x transcript: alternative codon and consequence
  pos <- merge(pos, cod, by = c("transcript_id", "codon"))
  snv <- merge(pos[, .(chrom, pos, strand, gene_name, gene_id, transcript_id, cds_pos, codon, frame, ref_codon)],
               cadd, by = c("chrom", "pos"), allow.cartesian = TRUE)
  snv[, calt := ifelse(strand == "+", alt, comp[alt])]
  snv[, alt_codon := paste0(ifelse(frame == 1L, calt, substr(ref_codon, 1, 1)),
                            ifelse(frame == 2L, calt, substr(ref_codon, 2, 2)),
                            ifelse(frame == 3L, calt, substr(ref_codon, 3, 3)))]
  snv[, `:=`(aa_ref = unname(code[ref_codon]), aa_alt = unname(code[alt_codon]))]
  snv[, consequence := fcase(
    codon == 1L & aa_ref == "M" & aa_alt != "M", "start_lost",
    aa_ref == "*" & aa_alt != "*",               "stop_lost",
    aa_ref == "*" & aa_alt == "*",               "stop_retained",
    aa_ref != "*" & aa_alt == "*",               "stop_gained",
    aa_ref == aa_alt,                            "synonymous",
    default =                                    "missense")]
  snv[, hgvsp := paste0("p.", aa3[aa_ref], codon, ifelse(aa_ref == aa_alt, "=", aa3[aa_alt]))]
  snv[, c("calt", "frame") := NULL]
  setnames(snv, "codon", "protein_pos")
  setcolorder(snv, c("chrom", "pos", "ref", "alt", "RawScore", "PHRED", "gene_name", "gene_id",
                     "transcript_id", "strand", "cds_pos", "protein_pos", "ref_codon", "alt_codon",
                     "aa_ref", "aa_alt", "hgvsp", "consequence"))
  saveRDS(snv, file.path(dir, sprintf("coding_snv_%s.rds", ch)))
  cat(sprintf("chr%-2s transcripts %5d | QC identical %5d | SNV x transcript rows %s\n", ch, nrow(qc),
              qc[identical == TRUE, .N], format(nrow(snv), big.mark = ",")))
  print(snv[, .N, by = consequence][order(-N)])
}
qc <- rbindlist(qc_all)
fwrite(qc, file.path(dir, "qc_translation.tsv"), sep = "\t")
cat("\nQC: transcripts identical to GENCODE protein:", qc[identical == TRUE, .N], "of", nrow(qc), "\n")
cat("Non-identical transcripts:\n"); print(qc[identical == FALSE | is.na(identical)])
