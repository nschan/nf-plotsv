# nf-syri

Here, whole genome alignments are created and then passed to syri.

# Usage

## Samplesheet

Samplesheet layout is 

```
name,fasta
```

Reference Name can be provided using `--reference`, reference genome path is `--ref_genome`

# Steps

minimap2

 Align query to reference

syri

 Use alignments