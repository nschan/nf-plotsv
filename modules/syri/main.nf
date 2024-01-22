include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process SYRI {
    tag "$meta"
    label 'process_low'
    publishDir "${params.out}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename,
                                        options:params.options, 
                                        publish_dir:"${task.process}".replace(':','/').toLowerCase(), 
                                        publish_id:meta) }
    input:
        tuple val(meta), path(alignment)
        val(reference)

    output:
        tuple val(meta), path("*syri.out"), emit: syri_out
        tuple val(meta), path("*syri.vcf"), emit: syri_vcf

    script:
        """
        syri -c -r ${reference} -q ${meta} -k -F B -ncores $task.cpus
        mv syri.out ${meta}_on_${reference}.syri.out
        mv syri.vcf ${meta}_on_${reference}.syri.vcf
        """
}

process SYRI_PAIRWISE {
    tag "${query}_${reference}"
    label 'process_low'
    publishDir "${params.out}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename,
                                        options:params.options, 
                                        publish_dir:"${task.process}".replace(':','/').toLowerCase(), 
                                        publish_id:meta) }
    input:
        tuple val(reference), val(query), path(alignment)

    output:
        tuple val("${query}_${reference}"), path("*syri.out"), emit: syri_out
        tuple val("${query}_${reference}"), path("*syri.vcf"), emit: syri_vcf

    script:
        """
        syri -c  ${alignment} -r ${reference} -q ${query} -k -F B -ncores $task.cpus
        mv syri.out ${query}_on_${reference}.syri.out
        mv syri.vcf ${query}_on_${reference}.syri.vcf
        """
}