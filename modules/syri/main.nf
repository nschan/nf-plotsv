process SYRI {
    fair true
    tag "${query}_on_${reference}"
    label 'process_low'
    publishDir(
      path: { "${params.out}/${task.process}".replace(':','/').toLowerCase() }, 
      mode: 'copy',
      overwrite: true,
      saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) }
    ) 
    input:
        tuple val(reference), path(reference_genome), val(query), path(query_genome), path(alignment), path(index)

    output:
        tuple val(query), path("*syri.out"), emit: syri_out
        tuple val(query), path("*syri.vcf"), emit: syri_vcf

    script:
        """
        micromamba run -n base syri -c ${alignment} -r ${reference_genome} -q ${query_genome} -k -F B --nc $task.cpus
        mv syri.out ${query}_on_${reference}.syri.out
        mv syri.vcf ${query}_on_${reference}.syri.vcf
        """
}

process SYRI_PAIRWISE {
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
        tuple val(name_A), path(genome_A), val(name_B), path(genome_B), path(alignment), path(index)

    output:
        tuple val(name_A), val(name_B), path("*syri.out"), emit: syri_out
        tuple val(name_A), val(name_B), path("*syri.vcf"), emit: syri_vcf

    script:
        """
        micromamba run -n base syri -c ${alignment} -r ${genome_B} -q ${genome_A} -k -F B --nc $task.cpus
        mv syri.out ${name_A}_on_${name_B}.syri.out
        mv syri.vcf ${name_A}_on_${name_B}.syri.vcf
        """
}