# Multiple sequence alignment of ATP7B coding sequences with the msa package (MUSCLE)
#
# Run from the assignment-1-atp7b folder, e.g.:  Rscript scripts/msa_muscle.R
# Packages: msa, Biostrings, ape (msa and Biostrings come from Bioconductor)

library(msa)
library(Biostrings)
library(ape)

fasta_file <- "data/orthologs_and_isoforms.fasta"
out_dir <- "results/msa"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# STEP 1: load the FASTA file
cat("Loading sequences from:", fasta_file, "\n")
sequences <- readDNAStringSet(fasta_file)

cat("\nSEQUENCE INFORMATION\n")
cat("Number of sequences:", length(sequences), "\n")
cat("Sequence names:", names(sequences), "\n")
cat("Sequence lengths:", width(sequences), "\n\n")

# STEP 2: multiple sequence alignment
# alignment <- msa(sequences, method = "ClustalW")
alignment <- msa(sequences, method = "Muscle")
cat("Alignment completed!\n\n")

# STEP 3: view the results
cat("=== MULTIPLE SEQUENCE ALIGNMENT ===\n")
print(alignment)

alignment_matrix <- as.matrix(alignment)
cat("Alignment length:", ncol(alignment_matrix), "\n")
cat("Number of sequences:", nrow(alignment_matrix), "\n")

# STEP 4: save the results
write.phylip(alignment, file.path(out_dir, "alignment_results.phy"))

alignment_fasta <- as(alignment, "DNAStringSet")
writeXStringSet(alignment_fasta, file.path(out_dir, "alignment_results.fasta"))

capture.output(print(alignment, show = "complete"),
               file = file.path(out_dir, "alignment_results.txt"))

cat("\n=== RESULTS SAVED in", out_dir, "===\n")
print(list.files(out_dir))
