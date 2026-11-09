
# Bioinformatics Coursework

Assignments from *Introduction to Bioinformatics* (University of Tehran). The work uses public databases and web tools (ClinVar, NCBI, Clustal Omega, BLAST, AlphaFold DB, SWISS-MODEL) together with short Python and R scripts.

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
```

## Run the scripts

```
pip install biopython
python assignment-1-atp7b/scripts/fetch_cds.py NM_000053.4
```

Set your own email in `Entrez.email` first. For the R script, install `msa`, `Biostrings` and `ape`, then from the `assignment-1-atp7b` folder run `Rscript scripts/msa_muscle.R`.

## Notes

Course assignments for Introduction to Bioinformatics, University of Tehran, Fall 2025.
