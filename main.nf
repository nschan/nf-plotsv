nextflow.enable.dsl = 2 

include { ALIGN_GENOMES } from './modules/align/main'
include { SYRI as RUN_SYRI } from './modules/syri/main' 
include { PLOTSR } from './modules/plotsr/main' 
include { FIXCHR } from './modules/fixchr/main'
include { SEQTK_ORIENT } from './modules/seqtk/main'
include { SEQTK_SUBSET as SUBSET } from './modules/seqtk/main'
include { ALIGN_PAIRWISE } from './modules/align/main'
include { SYRI_PAIRWISE } from './modules/syri/main' 
include { PLOTSR_PAIRWISE } from './modules/plotsr/main' 


/*
Samplesheet: 
name,fasta
*/
log.info """\
======================================================================================================================================================
======================================================================================================================================================
minimap2 -> syri -> plotsr
------------------------------------------------------------------------------------------------------------------------------------------------------
Niklas Schandry                                      niklas@bio.lmu.de                                                                    
------------------------------------------------------------------------------------------------------------------------------------------------------
  Results directory  : ${params.out}

  Parameters:
     samplesheet     : ${params.samplesheet}
     collect         : ${params.collect}
     reference       : ${params.reference}
     ref_genome      : ${params.ref_genome}
     pairwise        : ${params.pairwise}
     subset_pattern  : ${params.subset_pattern}
======================================================================================================================================================
======================================================================================================================================================
"""
    .stripIndent(false)
/*
  PREPARE GENOMES
  ---------------
    Subset to chromosomes 1-5
    branch into reference and non-referene geomes
    align genomes
    fix orientation
Prepare genomes uses FIXCHR to find inverted chromsomes, and then passes those to SEQTK_ORIENT to rc the inverted chromosomes..
*/

workflow PREPARE_GENOMES {
  take:
    input

  main:
    input
    | SUBSET
    | branch { row ->
        REF: row[0] == params.reference
        ASSEMBLY: row[0] != params.reference
      }
    | set { ch_branched }

    ALIGN_GENOMES(ch_branched.ASSEMBLY, tuple(params.reference, params.ref_genome))
    | FIXCHR
    //| SEQTK_ORIENT
    | map { it -> [name: it[0], path: it[1]]}
    | set { fixed }

    ch_branched.REF
    .concat(fixed)
    .set { fixed }


  emit:
    fixed
} 

/*
This workflow is simple, and largely an exercise in manipulating channels.
  In pairwise mode, the input channel is collated so that two consecutive rows become a tuple, which is then unnnested
  After passing everything through alingment and syri, the syri outputs are joined to input channel (to preserve the order)
  and then collected into a tuple which is deconstructed in PLOTSR_PAIRWISE into a list of arguments.
*/

workflow SYRI {
  ch_input = Channel.fromPath(params.samplesheet) 
              | splitCsv(header:true) 

  ch_input.map { it -> [name: it.name] }
   | set { ch_order }

  PREPARE_GENOMES(ch_input)
 
  if(params.pairwise) {
    ch_order
      .cross(PREPARE_GENOMES.out)
      .map { it -> it[1]}
      .collate(2, 1, false)
      .set { ch_chunked }
      

    ch_chunked.map { it -> [ref_name: it[0].name, ref_path: it[0].path, query_name: it[1].name, query_path: it[1].path]}
      .set { ch_chunked }


    ch_chunked.map { it -> [ref_name: it.ref_name, query_name: it.query_name] }
      .set { ch_chunk_order }

    ch_chunked
    | ALIGN_PAIRWISE
    | SYRI_PAIRWISE

    ch_chunk_order
      .cross(SYRI_PAIRWISE.out.syri_out)
      .map { it -> it[1][2] }
      .collect()
      .set { plotsr_in }
    
    PREPARE_GENOMES.out
      .map { it -> it.path }
      .collect()
      .set { ch_prepared_files }

    ch_order
      .map { it -> it.name }
      .flatten()
      .collect()
      .set { ch_names }

    PLOTSR_PAIRWISE(plotsr_in, ch_names, ch_prepared_files)
    
  } else {
    ALIGN_GENOMES(PREPARE_GENOMES.out, tuple(params.reference, params.ref_genome))
    RUN_SYRI(ALIGN_GENOMES.out, params.reference)
    PLOTSR(RUN_SYRI.out.syri_out, params.reference)
  }
}

workflow { SYRI() }