process SYRI {
    fair true
    tag "$meta"
    label 'process_low'
    publishDir(
      path: { "${params.out}/${task.process}".replace(':','/').toLowerCase() }, 
      mode: 'copy',
      overwrite: true,
      saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) }
    ) 
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