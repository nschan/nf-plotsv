# nf-syri

Here, whole genome alignments are created and then passed to syri.

# Usage

## Samplesheet

Samplesheet layout is 

```
name,fasta
```

Reference Name can be provided using `--reference`, reference genome path is `--ref_genome`

## Pairwise mode

To create a lot of consecutive pairwise alignments to create a plot across many genomes, use `--pairwise`.
This will create pairwise alignments from top to bottom of the samplesheet (i.e. align row2 on row1, row3 on row2, row4 on row3, etc)

# Steps

minimap2

 Align query to reference

syri

 Use alignments