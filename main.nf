nextflow.enable.dsl = 2 

include { SEQTK_ORIENT } from './modules/seqtk/main'
include { SEQTK_SUBSET as SUBSET } from './modules/seqtk/main'
include { ALIGN_GENOMES } from './modules/align/main'
include { FIXCHR } from './modules/fixchr/main'
include { SYRI } from './modules/syri/main' 
include { PLOTSR } from './modules/plotsr/main' 

// Pairwise modules
include { ALIGN_PAIRWISE } from './modules/align/main'
include { SYRI_PAIRWISE } from './modules/syri/main' 
include { PLOTSR_PAIRWISE } from './modules/plotsr/main' 


/*
Samplesheet: 
name,fasta
*/
log.info """\
=======================================================================================================================
=======================================================================================================================
                                                    nf-plotsv
                                                    ---------
                       Plot structural variation across genomes using the 'Schneeberger tools'
-----------------------------------------------------------------------------------------------------------------------
Niklas Schandry                                  niklas@bio.lmu.de                      gitlab.lrz.de/beckerlab/nf-syri      
-----------------------------------------------------------------------------------------------------------------------
  Results directory  : ${params.out}

  Parameters:
     samplesheet     : ${params.samplesheet}
     reference       : ${params.reference}
     ref_genome      : ${params.ref_genome}
     reorient        : ${params.reorient}
     pairwise        : ${params.pairwise}
     subset_pattern  : ${params.subset_pattern}
     plotsr config   : ${params.plotsr_conf}
     plotsr args     : ${params.plotsr_args}
     plotsr tracks   : ${params.plotsr_tracks}
=======================================================================================================================
=======================================================================================================================
"""
    .stripIndent(false)
/*
  PREPARE GENOMES
  ---------------
    Subset to whatever is the pattern
    branch into reference and non-referene geomes
    align genomes
    fix orientation
Prepare genomes uses FIXCHR to find inverted chromsomes, and then passes those to SEQTK_ORIENT to rc the inverted chromosomes..
*/

workflow PREPARE_GENOMES {
  take:
    input

  main:

    if(params.reorient) {
      input
        | SUBSET
        | branch { row ->
            REF: row[0] == params.reference
            ASSEMBLY: row[0] != params.reference }
        | set { ch_branched }
      ALIGN_GENOMES(ch_branched.ASSEMBLY, tuple(params.reference, params.ref_genome))
        | FIXCHR
        | map { it -> [name: it[0], path: it[1]]}
        | set { fixed }
      ch_branched.REF
        .concat(fixed)
        .set { fixed }
    } else {
          input
          | SUBSET
          | map { it -> [name: it[0], path: it[1]]}
          | set { fixed }
    }
    
  emit:
    fixed
} 

/*
This workflow is simple, and largely an exercise in manipulating channels.
  In pairwise mode, the input channel is collated so that two consecutive rows become a tuple, which is then unnnested
  After passing everything through alingment and syri, the syri outputs are joined to input channel (to preserve the order)
  and then collected into a tuple which is deconstructed in PLOTSR_PAIRWISE into a list of arguments.
  The nextflow part is okay, but naturally the PLOTSR_PAIRWISE process is pile sad bash.
*/

workflow PLOTSV {
  ch_input = Channel.fromPath(params.samplesheet) 
              | splitCsv(header:true) 

  ch_input.map { it -> [name: it.name] }
   | set { ch_order }

  PREPARE_GENOMES(ch_input)
 
  if(params.pairwise) {
    ch_order
      .cross(PREPARE_GENOMES.out)
      .map { it -> it[1] }
      .collate(2, 1, false)
      .set { ch_chunked }
      

    ch_chunked.map { it -> [name_A: it[0].name, genome_A: it[0].path, name_B: it[1].name, genome_B: it[1].path]}
      .set { ch_chunked }

    /*
    ch_chunked.map { it -> [[it.name_A, it.name_B]] }
      //.toList()
      //.dump(tag: 'chunk_order')
      .set { ch_chunk_order }
    */

    ch_chunked
    | ALIGN_PAIRWISE
    | SYRI_PAIRWISE
    

    SYRI_PAIRWISE
      .out
      .syri_out
      .map { it -> it[2] }
      .collect()
      .dump(tag: 'SYRI_out')
      .set { plotsr_in }
      
    ch_order
      .cross(PREPARE_GENOMES.out)
      .map { it -> it[1] }
      .map { it -> it.path }
      .collect()
      .set { ch_prepared_files }

    ch_order
      .map { it -> it.name }
      .flatten()
      .collect()
      //.dump(tag: 'plotsr_names')
      .set { ch_names }

    PLOTSR_PAIRWISE(plotsr_in, ch_names, ch_prepared_files, params.plotsr_conf, params.plotsr_args, params.plotsr_tracks, params.plotsr_colors)
    
  } else {
    ALIGN_GENOMES(PREPARE_GENOMES.out, tuple(params.reference, params.ref_genome))
    SYRI(ALIGN_GENOMES.out)
    PLOTSR(SYRI.out.syri_out, params.reference, params.plotsr_conf, params.plotsr_args)
  }
}

workflow { PLOTSV() }