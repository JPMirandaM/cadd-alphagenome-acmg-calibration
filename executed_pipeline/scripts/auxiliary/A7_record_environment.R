# Records the computational environment for publication:
#   package_inventory.tsv  every installed R package with version, repository and library path
# Usage: Rscript A7_record_environment.R <output_dir>
out <- commandArgs(trailingOnly = TRUE)[1]
ip <- installed.packages(fields = c("Repository", "biocViews"))
inv <- data.frame(
  package    = ip[, "Package"],
  version    = ip[, "Version"],
  repository = ifelse(!is.na(ip[, "Repository"]), ip[, "Repository"],
                      ifelse(!is.na(ip[, "biocViews"]), "Bioconductor", ifelse(ip[, "Priority"] %in% c("base", "recommended"), "R (bundled)", NA))),
  built      = ip[, "Built"],
  library    = ip[, "LibPath"],
  row.names = NULL
)
inv <- inv[order(tolower(inv$package)), ]
write.table(inv, file.path(out, "package_inventory.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("Packages recorded:", nrow(inv), "\n")
