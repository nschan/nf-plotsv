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
        tuple val(reference), val(query), path("*syri.out"), emit: syri_out
        tuple val(reference), val(query), path("*syri.vcf"), emit: syri_vcf

    script:
        """
        micromamba run -n base syri -c ${alignment} -r ${reference_genome} -q ${query_genome} -k -F B --nc $task.cpus
        mv syri.out ${query}_on_${reference}.syri.out
        mv syri.vcf ${query}_on_${reference}.syri.vcf
        """
}