# Part A: questions about the short-read FASTQ file
# Data: SRR36584168 (short-read WGS of a Klebsiella pneumoniae phage)
# Download with the SRA Toolkit:  fastq-dump SRR36584168
# Packages: ShortRead and Biostrings (Bioconductor), ggplot2
#
# Note: the original script was lost, so this version re-creates the analysis from the saved
# outputs in results/part_a. Re-run it to confirm that it gives the same results.

library(ShortRead)
library(Biostrings)
library(ggplot2)

fastq_file <- "SRR36584168.fastq"
reads <- readFastq(fastq_file)

# I. How many reads are in the file?
num_reads <- length(reads)
print(paste("Number of reads:", num_reads))

# II. Identifier, sequence and quality of the first read
first_read <- reads[1]
print("Identifier:")
print(id(first_read))
print("Sequence:")
print(sread(first_read))
print("Quality:")
print(quality(first_read))

# III. How many times does GCGATCAACTCGCTGG appear?
seqs <- sread(reads)
subseq <- DNAString("GCGATCAACTCGCTGG")
count <- vcountPattern(subseq, seqs, fixed = TRUE)
total_count <- sum(count)
print(paste("Subsequence appears", total_count, "times"))

# IV. Extract the first 1000 reads
first_1000 <- reads[1:1000]
writeFastq(first_1000, "first_1000_reads.fastq", compress = FALSE)

# V. Box plot of the base qualities by position (first 1000 reads)
qual <- as(quality(first_1000), "matrix")
qual_df <- data.frame(
  position = rep(seq_len(ncol(qual)), each = nrow(qual)),
  quality = as.vector(qual)
)
qual_df <- qual_df[!is.na(qual_df$quality), ]

ggplot(qual_df, aes(x = position, y = quality, group = position)) +
  geom_boxplot(color = "navy", outlier.size = 0.5) +
  labs(title = "Per Base Sequence Quality Box Plot",
       x = "Position in Read (bp)", y = "Phred Quality Score") +
  theme_minimal()

# VI. Distribution of read lengths (density plot)
read_lengths <- width(reads)

ggplot(data.frame(length = read_lengths), aes(x = length)) +
  geom_density(fill = "steelblue", alpha = 0.5) +
  geom_vline(xintercept = mean(read_lengths), color = "red", linetype = "dashed", linewidth = 1) +
  labs(title = "Read Length Distribution", x = "Read Length (bp)", y = "Density") +
  theme_minimal()
