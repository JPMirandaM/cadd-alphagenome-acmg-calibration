suppressPackageStartupMessages(library(data.table))
cv <- fread("C:/cadd/annot/clinvar.vcf.gz", sep = "\t", skip = "#CHROM", select = c(1, 2, 4, 5, 8),
            col.names = c("chr", "pos", "ref", "alt", "info"), colClasses = list(character = 1))
cv <- cv[nchar(ref) == 1 & nchar(alt) == 1 & ref %chin% c("A","C","G","T") & alt %chin% c("A","C","G","T")]
get <- function(k) { m <- regmatches(cv$info, regexpr(paste0("(^|;)", k, "=[^;]*"), cv$info)); out <- rep(NA_character_, nrow(cv)); out[grepl(paste0("(^|;)", k, "="), cv$info)] <- sub(paste0("^;?", k, "="), "", m); out }
cv[, sig := get("CLNSIG")][, rev := get("CLNREVSTAT")][, mc := get("MC")]
cv[, clase := fcase(sig %chin% c("Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic"), "P/LP",
                    sig %chin% c("Benign", "Likely_benign", "Benign/Likely_benign"), "B/LB", default = "otra")]
cv[, estrellas := fcase(grepl("practice_guideline", rev), "4", grepl("reviewed_by_expert_panel", rev), "3",
                        grepl("multiple_submitters,_no_conflicts", rev), "2",
                        grepl("single_submitter|conflicting_classifications", rev), "1", default = "0")]
cv[, consecuencia := fcase(grepl("missense_variant", mc), "missense",
                           grepl("nonsense|stop_gained", mc), "stop_gained",
                           grepl("splice_donor|splice_acceptor", mc), "splice_canonico",
                           grepl("synonymous", mc), "sinonima",
                           grepl("intron", mc), "intronica",
                           grepl("UTR", mc), "UTR",
                           default = "otra/sin_MC")]
cat("SNV en ClinVar:", nrow(cv), "\n\nP/LP y B/LB con >= 1 estrella, por consecuencia:\n")
print(dcast(cv[clase != "otra" & estrellas >= "1"], consecuencia ~ clase, fun.aggregate = length))
cat("\nMissense con >= 2 estrellas:\n")
print(cv[clase != "otra" & estrellas >= "2" & consecuencia == "missense", .N, by = clase])
