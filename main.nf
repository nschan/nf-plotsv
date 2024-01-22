nextflow.enable.dsl = 2 
params.publish_dir_mode = 'copy'
params.samplesheet = false
params.reference = 'Col-CEN_v1.2'
params.ref_genome = '/dss/dsslegfs01/pn73so/pn73so-dss-0000/becker_common/reference_genomes/Arabidopsis/Col-CEN/'

include { ALIGN_GENOMES} from './modules/align/main'
include { SYRI } from './modules/syri/main' 
include { PLOTSR } from './modules/plotsr/main' 

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
  RUN_SYRI()
}