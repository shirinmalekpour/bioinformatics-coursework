# Interpreting CIGAR strings

A CIGAR string describes how a read lines up with the reference. Each part is a number followed by an operation letter:

| Letter | Meaning | Uses read bases | Uses reference bases |
|---|---|---|---|
| `M` | aligned base (match or mismatch) | yes | yes |
| `=` | exact match | yes | yes |
| `X` | mismatch | yes | yes |
| `I` | insertion (bases in the read that are not in the reference) | yes | no |
| `D` | deletion (reference bases missing from the read) | no | yes |
| `N` | skipped region of the reference (for example an intron in RNA-seq) | no | yes |
| `S` | soft clip (bases kept in the read but not aligned) | yes | no |

## `29S21=1X25=`

29 bases at the start of the read are soft-clipped (present in the read but not aligned, for example adapter sequence or poor-quality bases), then 21 bases match the reference exactly, 1 base is a mismatch (a substitution or SNP), and 25 bases match exactly. The read is 76 bases long and covers 47 bases of the reference.

## `20M2I1M1D10M`

20 aligned bases, then a 2-base insertion (two bases in the read that are not in the reference), 1 aligned base, a 1-base deletion (one reference base that is missing from the read), and 10 aligned bases. The read is 33 bases long and covers 32 bases of the reference.

## `5M10N25M`

5 aligned bases, then 10 reference bases are skipped, then 25 aligned bases. The `N` operation is typical of RNA-seq, where a read spans a splice junction and the skipped stretch is an intron. The read is 30 bases long and covers 40 bases of the reference.
