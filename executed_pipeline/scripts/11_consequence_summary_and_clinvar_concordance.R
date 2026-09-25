# (1) Genome-wide counts of coding SNV consequences (MANE Select), (2) number of non-ATG start
# codons, and (3) concordance of our consequences with the ClinVar molecular consequence (MC)
# field for ClinVar SNVs located in MANE Select coding sequence.
suppressPackageStartupMessages(library(data.table))
dir <- "C:/cadd/data/coding"
snv <- rbindlist(lapply(list.files(dir, "^coding_snv_.*\\.rds$", full.names = TRUE), function(f)
  readRDS(f)[, .(chrom, pos, ref, alt, transcript_id, gene_name, protein_pos, ref_codon, consequence)]))
cat("SNV x transcript rows:", nrow(snv), "| unique SNVs:", uniqueN(snv, by = c("chrom", "pos", "alt")), "\n")
print(snv[, .N, by = consequence][order(-N)][, pct := round(100 * N / sum(N), 2)][])
st <- unique(snv[protein_pos == 1L, .(transcript_id, gene_name, ref_codon)])
cat("\nStart codons:", nrow(st), "| non-ATG:", st[ref_codon != "ATG", .N], "\n")
print(st[ref_codon != "ATG", .N, by = ref_codon][order(-N)])

# Variant-level consequence (a SNV in two overlapping genes can have two)
vc <- snv[, .(ours = paste(sort(unique(consequence)), collapse = "&")), by = .(chrom, pos, ref, alt)]

cv <- fread("C:/cadd/annot/clinvar.vcf.gz", sep = "\t", skip = "#CHROM", select = c(1, 2, 4, 5, 8),
            col.names = c("chrom", "pos", "ref", "alt", "info"), colClasses = list(character = 1))
cv <- cv[nchar(ref) == 1 & nchar(alt) == 1]
cv[, mc := sub(".*(^|;)MC=([^;]*).*", "\\2", info)][!grepl("(^|;)MC=", info), mc := NA]
cv[, info := NULL]
m <- merge(vc, cv, by = c("chrom", "pos", "ref", "alt"))
cat("\nClinVar SNVs in MANE coding sequence:", nrow(m), "| with MC:", m[!is.na(mc), .N], "\n")
so <- c(missense = "missense_variant", synonymous = "synonymous_variant", stop_gained = "nonsense",
        stop_lost = "stop_lost", start_lost = "initiator_codon_variant|start_lost",
        stop_retained = "stop_retained|synonymous_variant")
m <- m[!is.na(mc)]
m[, main := sub("&.*", "", ours)]
m <- m[main %in% names(so)]
m[, agree := mapply(function(o, x) grepl(so[[o]], x), main, mc)]
res <- m[, .(n = .N, agree = sum(agree), pct_agree = round(100 * mean(agree), 2)), by = .(ours = main)][order(-n)]
print(res)
cat("Overall concordance:", round(100 * mean(m$agree), 2), "%\n")
cat("\nMost frequent ClinVar MC among discordant missense calls:\n")
print(m[main == "missense" & !agree, .N, by = mc][order(-N)][1:10])
fwrite(res, file.path(dir, "clinvar_mc_concordance.tsv"), sep = "\t")
