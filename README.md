# nf-syri

Here, whole genome alignments are created and then passed to syri and then to plotsr.

# Usage

## Samplesheet

Samplesheet layout is 

```
name,fasta
```

Reference name can be provided using `--reference`, reference genome path is `--ref_genome`

## Pairwise mode

`--pairwise`: create consecutive pairwise alignments from the samplesheet to create a plot across many genomes.
This will create pairwise alignments from top to bottom of the samplesheet (i.e. align row2 on row1, row3 on row2, row4 on row3, etc) and then create a plot using plotsr.

# Steps

minimap2

 Align query to reference

syri

 Use alignments

plotsr
 
 plot syri output