include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process SYRI {
    fair true
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
        micromamba run -n base syri -c -r ${reference} -q ${meta} -k -F B -ncores $task.cpus
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
                                        publish_id:"${query}_${reference}") }
    input:
        tuple val(reference), path(reference_genome), val(query), path(query_genome), path(alignment), path(index)

    output:
        tuple val(reference), val(query), path("*syri.out"), emit: syri_out
        tuple val(reference), val(query), path("*syri.vcf"), emit: syri_vcf

    script:
        """
        micromamba run -n base syri -c ${alignment} -r ${reference_genome} -q ${query_genome} -k -F B --nc $task.cpus
        mv syri.out ${query}_on_${reference}.syri.out
        mv syri.vcf ${query}_on_${reference}.syri.vcf
        """
}