include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process PLOTSR {
    tag "$meta"
    label 'process_low'
    publishDir "${params.out}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename,
                                        options:params.options, 
                                        publish_dir:"${task.process}".replace(':','/').toLowerCase(), 
                                        publish_id:meta) }
    input:
        tuple val(meta), path(syri_out)
        val(reference)

    output:
        tuple val(meta), path("*plotsr.pdf"), emit: figure

    script:
        """
        plotsr ${meta}_on_${reference}.syri.out $reference $meta -H 8 -W 5 -o ${meta}_on_${reference}.plotsr.pdf
        """
}