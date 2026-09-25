# Corrects the consequence of coding SNVs at three kinds of special codons identified by the
# translation QC (auxiliary A12, qc_translation_mismatches.tsv):
#   1. Selenocysteine UGA codons (GENCODE protein has U): reference residue is Sec (U); an SNV is
#      stop_gained if the new codon is a stop, otherwise missense (Sec -> X).
#   2. Start codons (codon 1), including non-AUG initiation codons: any SNV is start_lost
#      (Sequence Ontology / Ensembl VEP definition).
#   3. Incomplete terminal codons (< 3 nt): incomplete_terminal_codon_variant.
# Files are modified in place (coding_snv_<chrom>.rds); a table of every changed row is written.
suppressPackageStartupMessages({ library(data.table); library(Biostrings) })
dir <- "C:/cadd/data/coding"
mm  <- fread(file.path(dir, "qc_translation_mismatches.tsv"))
sec <- unique(mm[cause == "selenocysteine (UGA recoded)", .(transcript_id, protein_pos)])
aa3 <- c(A="Ala",R="Arg",N="Asn",D="Asp",C="Cys",Q="Gln",E="Glu",G="Gly",H="His",I="Ile",L="Leu",K="Lys",
         M="Met",F="Phe",P="Pro",S="Ser",T="Thr",W="Trp",Y="Tyr",V="Val","*"="Ter",U="Sec")
changes <- list()
for (f in list.files(dir, "^coding_snv_.*\\.rds$", full.names = TRUE)) {
  d <- readRDS(f)
  d[, row_id := .I]
  before <- d[, .(row_id, consequence_before = consequence, hgvsp_before = hgvsp)]
  # 1. Selenocysteine
  s <- d[sec, on = .(transcript_id, protein_pos), nomatch = NULL, which = TRUE]
  if (length(s)) {
    d[s, aa_ref := "U"]
    d[s, consequence := ifelse(aa_alt == "*", "stop_gained", "missense")]
    d[s, hgvsp := paste0("p.Sec", protein_pos, aa3[aa_alt])]
  }
  # 2. Start codon: any SNV is start_lost
  d[protein_pos == 1L, consequence := "start_lost"]
  # 3. Incomplete terminal codon
  d[nchar(ref_codon) < 3, `:=`(consequence = "incomplete_terminal_codon_variant", aa_ref = NA, aa_alt = NA, hgvsp = NA)]
  ch <- merge(before, d[, .(row_id, chrom, pos, ref, alt, gene_name, transcript_id, protein_pos, ref_codon,
                            alt_codon, consequence, hgvsp)], by = "row_id")[consequence != consequence_before |
                                                                                 hgvsp != hgvsp_before |
                                                                                 xor(is.na(hgvsp), is.na(hgvsp_before))]
  changes[[basename(f)]] <- ch
  d[, row_id := NULL]
  n_na <- d[is.na(aa_ref) & consequence != "incomplete_terminal_codon_variant", .N]
  if (n_na) stop("Unexplained NA residues in ", basename(f))
  saveRDS(d, f)
}
ch <- rbindlist(changes)
fwrite(ch, file.path(dir, "consequence_corrections.tsv"), sep = "\t")
cat("Rows corrected:", nrow(ch), "\n")
print(ch[, .N, by = .(consequence_before, consequence)][order(-N)])
