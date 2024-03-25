include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process ALIGN_GENOMES {
    fair true
    tag "$query"
    label 'process_low'
    publishDir "${params.out}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename,
                                        options:params.options, 
                                        publish_dir:"${task.process}".replace(':','/').toLowerCase(), 
                                        publish_id:meta) }
    input:
        tuple val(query), path(query_genome)
        tuple val(reference), path(reference_genome)

    output:
        tuple val(reference), path(reference_genome), val(query), path(query_genome), path("*.bam"), emit: alignment

    script:
        """
        minimap2 -t $task.cpus \\
            -ax asm5 \\
            --eqx \\
            ${reference_genome} ${query_genome} \\
            | samtools sort -O BAM -@ $task.cpus  > ${query}_on_${reference}.bam
        """
}

process ALIGN_PAIRWISE {
    fair true
    tag "${query}_${reference}"
    label 'process_low'
    publishDir "${params.out}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename,
                                        options:params.options, 
                                        publish_dir:"${task.process}".replace(':','/').toLowerCase(), 
                                        publish_id:"${query}_${reference}") }
    input:
        tuple val(reference), path(reference_genome), val(query), path(query_genome)

    output:
        tuple val(reference), path(reference_genome), val(query), path(query_genome), path("*.bam"), path("*.bai"), emit: alignment

    script:
        """
        minimap2 -t $task.cpus \\
            -ax asm5 \\
            --eqx \\
            ${reference_genome} ${query_genome} \\
            | samtools sort -O BAM -@ $task.cpus > ${query}_on_${reference}.bam
            samtools index ${query}_on_${reference}.bam
        """
}