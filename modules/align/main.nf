process ALIGN_GENOMES {
    fair true
    tag "$query"
    label 'process_low'
    publishDir(
      path: { "${params.out}/${task.process}".replace(':','/').toLowerCase() }, 
      mode: 'copy',
      overwrite: true,
      saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) }
    ) 
    input:
        tuple val(query), path(query_genome)
        tuple val(reference), path(reference_genome)

    output:
        tuple val(reference), path(reference_genome), val(query), path(query_genome), path("*.bam"),path("*.bai"), emit: alignment

    script:
        """
        minimap2 -t $task.cpus \\
            -ax asm5 \\
            --eqx \\
            ${reference_genome} ${query_genome} \\
            | samtools sort -O BAM -@ $task.cpus  > ${query}_on_${reference}.bam
        samtools index ${query}_on_${reference}.bam
        """
}

process ALIGN_PAIRWISE {
    fair true
    tag "${name_A}_on_${name_B}"
    label 'process_low'
    publishDir(
      path: { "${params.out}/${task.process}".replace(':','/').toLowerCase() }, 
      mode: 'copy',
      overwrite: true,
      saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) }
    ) 
    input:
        tuple val(name_A), path(genome_A), val(name_B), path(genome_B)

    output:
        tuple val(name_A), path(genome_A), val(name_B), path(genome_B), path("*.bam"), path("*.bai"), emit: alignment

    script:
        """
        minimap2 -t $task.cpus \\
            -ax asm5 \\
            --eqx \\
            ${genome_A} ${genome_B} \\
            | samtools sort -O BAM -@ $task.cpus > ${name_B}_on_${name_A}.bam
            samtools index ${name_B}_on_${name_A}.bam
        """
}