# nf-plotsv

This [`nextflow`](https://nextflow.io) pipeline can be used to analyze genomes using [`syri`](https://schneebergerlab.github.io/syri/) and [`plotsr`](https://github.com/schneebergerlab/plotsr/). Whole genome alignments are created using minimap2, oriented and re-alinged if needed, and then passed to syri and then to plotsr.
This [`nextflow`](https://nextflow.io) pipeline can be used to analyze genomes using [`syri`](https://schneebergerlab.github.io/syri/) and [`plotsr`](https://github.com/schneebergerlab/plotsr/). Whole genome alignments are created using minimap2, oriented and re-alinged if needed, and then passed to syri and then to plotsr.

# Usage

After cloning into `$HOME` this pipeline can be run locally with `conda` environments like this:

```
nextflow run ~/nf-syri --samplesheet samplesheet.csv -profile local,conda
```

Or locally with docker

```
nextflow run ~/nf-syri --samplesheet samplesheet.csv -profile local,docker
```

Or on biohpc_gen@LRZ using charliecloud

```
nextflow run ~/nf-syri --samplesheet samplesheet.csv -profile biohpc_gen,charliecloud
```

To run on other infrastructure, a matching config needs to be created (or copied from [`nf-core/configs`](https://github.com/nf-core/configs/tree/master/conf))

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
| plotsr_args | Extra arguments for plotsr | `'-S 0.7 -W 5 -H 5 -f 12'` |
| plotsr_tracks | path to tracks file for plotsr, see [here](https://github.com/schneebergerlab/plotsr/blob/master/README.md#visualising-tracks). | `''` |
| plotsr_colors | Color palette to use. | [alphabet palette](https://kwstat.github.io/pals/reference/discrete.html) |


> plotsr_tracks needs to be a full path, e.g. something like "$PWD/tracks.txt" should work.
> plotsr_colors need to be provided as hex or something, see nextflow.config

### Pairwise mode

`--pairwise`: create consecutive pairwise alignments from the samplesheet to create a plot across many genomes.
This will create pairwise alignments from top to bottom of the samplesheet (i.e. align row2 on row1, row3 on row2, row4 on row3, etc) and then create a _single_ plot using plotsr.

### plotsr

`plotsr` is controlled via a config file. Defaults to the one in `assets/`

# Pipeline graph

```mermaid
%%{init: {'theme': 'dark',
          'themeVariables':{
            'commitLabelColor': '#cccccc',
            'commitLabelBackground': '#434244',
            'commitLabelFontSize': '12px',
            'tagLabelFontSize': '12px',
            'git0': '#8db7d2',
            'git1': '#58508d',
            'git2': '#bc5090',
            'git3': '#ff6361',
            'git4': '#ffa600',
            'git5': '#74a892',
            'git6': '#d69e49',
            'git7': '#00ffff'
            },
          'gitGraph': {
            'mainBranchName': "Prepare Genome",
             'parallelCommits': false
             } 
          }
}%%

gitGraph TB:
  commit id: "Subset by name [--subset_pattern]"
  commit id: "Alignment"
  commit id: "Identify / reorient misoriented [--reorient]"
  branch "Align"
  checkout "Align"
  commit id: "Align sequences"
  branch "SyRi / plotsr"
  checkout "SyRi / plotsr"
  commit id: "Run SyRi"
  commit id: "Plot results"
```

# Pipeline graph

```mermaid
%%{init: {'theme': 'dark',
          'themeVariables':{
            'commitLabelColor': '#cccccc',
            'commitLabelBackground': '#434244',
            'commitLabelFontSize': '12px',
            'tagLabelFontSize': '12px',
            'git0': '#8db7d2',
            'git1': '#58508d',
            'git2': '#bc5090',
            'git3': '#ff6361',
            'git4': '#ffa600',
            'git5': '#74a892',
            'git6': '#d69e49',
            'git7': '#00ffff'
            },
          'gitGraph': {
            'mainBranchName': "Prepare Genome",
             'parallelCommits': false
             } 
          }
}%%

gitGraph TB:
  commit id: "Subset by name"
  branch "Reorient"
  commit id: "Alignment"
  branch "Align"
  commit id: "Identify misoriented"
  checkout "Align"
  merge "Reorient" tag: "Oriented sequences"
  commit id: "Align sequences"
  branch "SyRi / plotsr"
  checkout "SyRi / plotsr"
  commit id: "Run SyRi"
  commit id: "Plot results"
```

# Software used

This pipeline uses the following packages:

 - [`minimap2`](https://github.com/lh3/minimap2)

 - [`seqtk`](https://github.com/lh3/seqtk)

 - [`seqkit`](https://bioinf.shenwei.me/seqkit/)

 - [`syri`](https://schneebergerlab.github.io/syri/)
 
 - [`fixchr`](https://github.com/schneebergerlab/fixchr)

 - [`plotsr`](https://github.com/schneebergerlab/plotsr)
