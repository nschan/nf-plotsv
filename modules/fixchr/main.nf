process FIXCHR {
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
        tuple val(reference), path(reference_genome), val(query), path(query_genome), path(alignment), path(index)

    output:
        tuple val(query), path("*qry.filtered.fa"), emit: fixed_query

    script:
        """
        micromamba run -n base fixchr -c ${alignment} -F B -r ${reference_genome} -q ${query_genome} --prefix ${query} --log INFO --contig_size 5000 
        """
}