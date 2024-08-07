# nf-plotsv

Here, whole genome alignments are created, oriented if needed, and then passed to syri and then to plotsr. The workflow is called `plotsv`, but the repo is named `nf-syri` and I am too lazy to make new one..

# Usage

After cloning into home:

`nextflow run ~/nf-syri --samplesheet samplesheet.csv -profile biohpc_gen,charliecloud`

## Samplesheet

Samplesheet layout is 

```
name,fasta
```

Reference name can be provided using `--reference`, reference genome path is `--ref_genome`.
> Note: The reference genome needs to be prefiltered to contain the same chromsomes as those selected using `params.subset_pattern` when **not** using pairwise mode

## Params

Default params are defined in [`nextflow.config`](nextflow.config):

| Parameter | Effect | Default |
|  ---  |  ---   |   ---   |
| samplesheet | Samplesheet to be used | `false` |
| reference  | Reference Name | `'Col-CEN_v1.2'` |
| ref_genome | Referemce genome fasta | `'/dss/dsslegfs01/pn73so/pn73so-dss-0000/becker_common/reference_genomes/Arabidopsis/Col-CEN/Col-CEN_v1.2.fasta'` |
| reorient | Reorient sequences to have them all go the same direction? | `false` |
| pairwise | Use pairwise mode (see below) | `true` |
| subset_pattern | Pattern used for subsetting genomes in samplesheet | `"Chr[1-5]"` |
| out | Directory for results | `'./'` |
| plotsr_conf | Config file for plotsr | `"$projectDir/assets/plotsr_config.conf"` |


## Pairwise mode

`--pairwise`: create consecutive pairwise alignments from the samplesheet to create a plot across many genomes.
This will create pairwise alignments from top to bottom of the samplesheet (i.e. align row2 on row1, row3 on row2, row4 on row3, etc) and then create a _single_ plot using plotsr.

## plotsr

`plotsr` is controlled via a config file. Defaults to the one in `assets/`

# Steps

minimap2 / samtools:
 -  Align queries to reference

fixchr:
 - Orient chromosomes

minimap2 / samtools:
 - Align correctly oriented chromosomes

syri:
 - Use alignments

plotsr:
 - plot syri output