[![Build Status](https://travis-ci.org/etiennedanis/problem-set-2.svg?branch=master)](https://travis-ci.org/etiennedanis/problem-set-2)

# Problem Set 2

| **Due**: 5pm on Feb 16

## Workflow

See [here](https://github.com/MOLB7621/Discussion/issues/9) for the
general workflow.

Your answers should be in YAML format in a file called `answers.yml` at
the top level of the repository.

```
# answers.yml should look like this:
answer-1: 123
answer-2: 456
```

For this problem set you will use BEDtools to do some calculations. **Each
problem is worth 10 points**.

These files are in the `data/` directory.

- Fasta file with human genome sequence: `fasta/hg19.chr22.fa`
- BED file containing ChIP-seq peaks for H3K4me3 in Hela cells:
  `bed/encode.h3k4me3.hela.chr22.bed.gz`
- BED file with all genes in hg19: `bed/genes.hg19.bed.gz`.
- File containing peak calls for ENCODE transcription factor ChIP-seq
  experiements: `bed/encode.tfbs.chr22.bed.gz`.
- Bedgraph with CTCF ChIP-seq data in bedGraph format: `bedgraph/ctcf.hela.chr22.bg.gz`
- A "genome file" with chromosome size info: `genome/hg19.genome`
- A file containing transcription start sites (TSS) for `chr22`: `bed/tss.hg19.chr22.bed.gz`

## Question 1

Use BEDtools intersect to identify the size of the largest overlap between
CTCF and H3K4me3 locations.

## Question 2

Use BEDtools to calculate the GC content of nucleotides 19,000,000 to
19,000,500 on chr22 of `hg19` genome build. Report the GC content
as a fraction (e.g., 0.50).

## Question 3

Use BEDtools to identify the length of the CTCF ChIP-seq peak (i.e.,
interval) that has the largest mean signal in `ctcf.hela.chr22.bg.gz`.

## Question 4

Use BEDtools to identify the gene promoter (defined as 1000 bp upstream of
a TSS) with the highest median signal in `ctcf.hela.chr22.bg.gz`.  Report
the gene name (e.g., 'ABC123')

## Question 5

Use BEDtools to identify the longest interval on `chr22` that is not
covered by `genes.hg19.bed.gz`. Report the interval like `chr1:100-500`.

## Question 6 (extra credit)

Use one or more BEDtools that we haven't covered in class. Be creative.

