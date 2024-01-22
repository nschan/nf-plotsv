nextflow.enable.dsl = 2 
params.publish_dir_mode = 'copy'
params.samplesheet = false
params.reference = 'Col-CEN_v1.2'
params.ref_genome = '/dss/dsslegfs01/pn73so/pn73so-dss-0000/becker_common/reference_genomes/Arabidopsis/Col-CEN/Col-CEN_v1.2.fasta'
params.pairwise = false
params.out = './results'

include { ALIGN_GENOMES} from './modules/align/main'
include { SYRI } from './modules/syri/main' 
include { PLOTSR } from './modules/plotsr/main' 
include { ALIGN_PAIRWISE } from './modules/align/main'
include { SYRI_PAIRWISE } from './modules/syri/main' 
include { PLOTSR_PAIRWISE } from './modules/plotsr/main' 

/*
Samplesheet: 
name,fasta
*/

workflow WGA {
  take: ch_input
  
  main:
    ALIGN_GENOMES(ch_input, params.ref_genome)
    alignment = ALIGN_GENOMES.out
  emit:
    alignment
}

workflow VARIATION {
  take: alignment
  
  main:
    SYRI(alignment, params.reference)
    PLOTSR(SYRI.out.syri_out, params.reference)
}

workflow RUN_SYRI {
    if(params.samplesheet) {
        ch_input = Channel.fromPath(params.samplesheet) 
                        .splitCsv(header:true) 
    }
    else {
    exit 1, 'Input samplesheet not specified!'
    }
    WGA(ch_input) 
    VARIATION(WGA.out)
}

workflow {
  ch_input = Channel.fromPath(params.samplesheet) 
                        .splitCsv(header:true) 
  if(params.pairwise) {
    if(params.samplesheet) {
       ch_chunked = ch_input | collate(2, 1, false)
    }
    ch_chunked.map { row -> [ref_name = row[0].name, ref_path = row[0].path, query_name = row[1].name, query_path = row[1].path]}
    .set {ch_chunked}
    ALIGN_PAIRWISE(ch_chunked) 
    SYRI_PAIRWISE(ALIGN_PAIRWISE.out)
    ch_order = ch_chunked.map {it -> [reference = it.ref_name, query = it.query_name]}
    ch_syri = SYRI_PAIRWISE.out.syri_out
              .collect()
    ch_order.join(ch_syri, by: [0,1])
      .set {ch_syri}
    PLOTSR_PAIRWISE(ch_syri)
  } else {
  RUN_SYRI()
  }
}