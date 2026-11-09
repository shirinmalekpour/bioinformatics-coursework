# Bioinformatics Coursework

Assignments from *Introduction to Bioinformatics* (University of Tehran). The work uses public databases and web tools (ClinVar, NCBI, Clustal Omega, BLAST, AlphaFold DB, SWISS-MODEL) together with Python, R and command-line tools (SPAdes, QUAST, BWA, samtools).

- [Assignment 1: SNP analysis, comparative genomics and protein structure of *ATP7B*](#assignment-1-snp-analysis-comparative-genomics-and-protein-structure-of-atp7b)
- [Assignment 2: Read mapping and genome assembly of a *Klebsiella pneumoniae* phage](#assignment-2-read-mapping-and-genome-assembly-of-a-klebsiella-pneumoniae-phage)

## Assignment 1: SNP analysis, comparative genomics and protein structure of *ATP7B*

**Gene:** *ATP7B*, a copper-transporting ATPase. Mutations in this gene cause Wilson's disease.

Folder: [`assignment-1-atp7b`](assignment-1-atp7b) | Full write-up: [`Bioinformatics assignment 2.pdf`](assignment-1-atp7b/report/assignment-1-report.pdf)

### 1. SNPs from ClinVar

Three single nucleotide variants, each on a different RefSeq transcript:

| ClinVar ID | rsID | HGVS (coding) | Protein change | Consequence | Clinical significance |
|---|---|---|---|---|---|
| 3067976 | rs749171049 | NM_000053.4:c.4103T>C | p.Leu1368Pro | Missense | Likely pathogenic |
| 617987 | rs965523214 | NM_001406511.1:c.-54-62G>A | none | Intron variant | Uncertain significance |
| 538871 | rs28362532 | NM_001406512.1:c.-143T>C | none | 5' UTR variant | Benign |

Transcripts analysed: NM_000053.4 (variant 1), NM_001406511.1 (variant 6) and NM_001406512.1 (variant 7).

### 2. Retrieving coding sequences

The CDS of each transcript was downloaded from NCBI Nucleotide. One transcript was also fetched programmatically with Biopython and the Entrez API ([`scripts/fetch_cds.py`](assignment-1-atp7b/scripts/fetch_cds.py)).

### 3. Multiple sequence alignment and phylogeny

- **Clustal Omega** (EMBL-EBI) on the three human CDS: all three are 4,398 bp and align without a single difference, which means the transcript variants differ only in their untranslated regions. The two non-coding SNPs above therefore fall outside the aligned region.
- **MUSCLE in R** ([`scripts/msa_muscle.R`](assignment-1-atp7b/scripts/msa_muscle.R)) on four sequences: human variants 6 and 7, dog (*Canis lupus familiaris*) and rat (*Rattus norvegicus*) *ATP7B* CDS. Pairwise nucleotide identity over the aligned CDS:

| | Human v7 | Dog | Rat |
|---|---|---|---|
| **Human v6** | 100% | 85.9% | 83.7% |
| **Human v7** | | 85.9% | 83.7% |
| **Dog** | | | 80.6% |

The phylogenetic tree places the human transcripts with dog, and rat on a separate branch.

### 4. Pairwise protein alignment (BLASTp)

Isoform a (NP_001393441.1) against isoform c (NP_001230111.1), with five scoring matrices (BLOSUM62, BLOSUM80, BLOSUM90, PAM30, PAM70). BLOSUM90 gave the highest maximum score, which fits two very similar sequences. In the BLAST hit table, the top alignment covers 1,332 residues with 86.4% identity and 89.0% positives (E-value 0.0). The alignment repeats a GMTCQSCV-type motif, the copper-binding motif of the metal-binding domains of ATP7B.

### 5. Protein structure

AlphaFold/ColabFold was not available, so structures were built with SWISS-MODEL (human isoform a, human isoform c, mouse Atp7b), and the matching AlphaFold DB models were downloaded for human (P35670) and mouse (Q64446). Ramachandran statistics came from MolProbity.

| Protein | Mean pLDDT | Residues with pLDDT > 90 | Ramachandran favored | Ramachandran outliers |
|---|---|---|---|---|
| Human isoform a | 71.7 | 10.8% | 87.6% | 5.95% |
| Mouse Atp7b | 71.5 | 10.4% | 88.4% | 6.71% |
| Human isoform c | 72.0 | 10.0% | 89.2% | 6.14% |

pLDDT and Ramachandran values are as recorded in the report. TM-align was also unavailable, so isoform c was superimposed on isoform a with PDBeFold (SSM): Q-score 0.409, RMSD 0.789 Å over 931 aligned residues, 98.5% sequence identity. The Q-score is PDBeFold's own measure and is not the same as the TM-score.

## Assignment 2: Read mapping and genome assembly of a *Klebsiella pneumoniae* phage

Short-read whole-genome sequencing data of a *Klebsiella pneumoniae* phage (SRA accession [SRR36584168](https://www.ncbi.nlm.nih.gov/sra/SRR36584168)) is checked for quality, assembled de novo, and mapped to a reference genome. Parts B and C were run in Google Colab.

Folder: [`assignment-2-read-mapping-assembly`](assignment-2-read-mapping-assembly)

### Part A: Data and quality control

The reads were downloaded with `fastq-dump` from the SRA Toolkit. The questions about the file were answered in R with ShortRead, Biostrings and ggplot2 ([`scripts/part_a_analysis.R`](assignment-2-read-mapping-assembly/scripts/part_a_analysis.R)), and the quality control plots were made in Python because FastQC was not available. Screenshots and plots are in [`results/part_a`](assignment-2-read-mapping-assembly/results/part_a).

| Question | Result |
|---|---|
| Number of reads | 2,421,219 |
| First read | `SRR36584168.1 A01757:294:HFYJGDRX7:2:2101:2826:1000 length=251` |
| Occurrences of `GCGATCAACTCGCTGG` | 7,499 |
| Read length | Most reads are 250 to 251 bp (mean 244 bp), with very few shorter than 200 bp |
| Quality per read | Mean Phred score peaks at Q36 to Q37, and very few reads fall below Q30 |
| Quality along the read | 36.6 at the start of the read, falling gently to about 35 at position 250 |

In the first 1,000 reads (saved in [`first_1000_reads.fastq`](assignment-2-read-mapping-assembly/results/part_a/first_1000_reads.fastq)) the base qualities take only three values (Q11, Q25 and Q37), which is why the box plot looks flat. Quality is high throughout, and the slight drop towards the 3' end is typical of Illumina chemistry. Trimming the read ends would be an option for sensitive analyses.

### Part B: De novo assembly

The reads were assembled with SPAdes 3.13.1 as single-end reads (`--careful`, automatic k-mer selection).

| | Contigs | Scaffolds |
|---|---|---|
| Sequences | 54 | 50 |
| Total length | 149,588 bp | 149,988 bp |
| Largest | 30,788 bp | 30,788 bp |
| N50 | 8,023 bp | 16,944 bp |
| Sequences of at least 1 kb | 19 | 15 |

The contigs fall into two groups by k-mer coverage: three large contigs have about 5,000 to 6,600x coverage, while the others have roughly 6 to 17x.

### Part C: Read mapping

The reads were mapped to the reference provided by the instructor (32 sequences, about 139 kb in total) with BWA-MEM, then the SAM file was converted to a sorted, indexed BAM file with samtools.

| Metric | Result |
|---|---|
| Primary reads | 2,421,219 |
| Primary reads mapped | 2,419,790 (99.94%) |
| Supplementary alignments | 7,897 |

- **IGV:** a 195 bp window on sequence `k141_9` ([screenshot](assignment-2-read-mapping-assembly/results/mapping/igv_k141_9_193-387.png)) shows the read stack with a few mismatches and short deletions.
- **CIGAR strings:** interpretations of `29S21=1X25=`, `20M2I1M1D10M` and `5M10N25M` are in [`cigar_interpretation.md`](assignment-2-read-mapping-assembly/results/mapping/cigar_interpretation.md).
- The FASTQ, SAM and BAM files are too large for GitHub and are not included, and the reference genome is course material and is not included either.

### Commands

```
fastq-dump SRR36584168
spades.py -s SRR36584168.fastq -o spades_output --careful -t 4 --memory 8
bwa index reference.fasta
bwa mem -t 4 -R '@RG\tID:phage\tSM:Klebsiella_phage\tLB:lib1\tPL:ILLUMINA' reference.fasta SRR36584168.fastq > aligned_reads.sam
samtools view -bS aligned_reads.sam > aligned_reads.bam
samtools sort aligned_reads.bam -o aligned_reads_sorted.bam
samtools index aligned_reads_sorted.bam
samtools flagstat aligned_reads_sorted.bam
```

## Repository layout

```
assignment-1-atp7b/
├── data/       RefSeq CDS FASTA files and the ortholog/isoform FASTA used for MSA
├── scripts/    fetch_cds.py (Python), msa_muscle.R (R)
├── results/
│   ├── msa/         Clustal Omega alignment, R/MUSCLE alignment (FASTA, PHYLIP, text)
│   ├── blast/       BLASTp hit table and subject protein sequence
│   └── structures/  SWISS-MODEL and AlphaFold PDB files, Ramachandran plot
└── report/     PDF report with screenshots and discussion

assignment-2-read-mapping-assembly/
├── scripts/    part_a_analysis.R
└── results/
    ├── part_a/     screenshots, quality plots, first 1000 reads
    ├── assembly/   SPAdes contigs, scaffolds and log
    └── mapping/    flagstat summary, IGV screenshot, CIGAR interpretation
```

## Run the scripts

```
pip install biopython
python assignment-1-atp7b/scripts/fetch_cds.py NM_000053.4
```

Set your own email in `Entrez.email` first. For the R script, install `msa`, `Biostrings` and `ape`, then from the `assignment-1-atp7b` folder run `Rscript scripts/msa_muscle.R`.

## Notes

Course assignments for Introduction to Bioinformatics, University of Tehran, 2025 to 2026. Assignment 2 was completed as a group project.
