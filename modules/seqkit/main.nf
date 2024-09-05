process SEQKIT_GET_LENGTH {
    tag "$meta"
    label 'process_medium'
    fair true
    publishDir(
      path: { "${params.out}/${task.process}".replace(':','/').toLowerCase() }, 
      mode: 'copy',
      overwrite: true,
      saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) }
    ) 
    input:
      tuple val(meta), path(genome_fasta)
  
    output:
      tuple val(meta), path("*length.txt"), emit: length

    script:
      def prefix = task.ext.prefix ?: "${meta}"
      
  """
  seqkit fx2tab --length --name ${genome_fasta} > ${meta}_length.txt
  """
}