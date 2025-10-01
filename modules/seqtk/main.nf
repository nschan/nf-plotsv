process SEQTK_ORIENT {
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
        tuple val(reference), val(query), path(query_genome), path(alignment_info)

    output:
        tuple val(query), path("*oriented.fa"), emit: oriented

    script:
        """
        cat ${alignment_info} | awk '{ if (\$9 == -1) {print \$11}}' > inverted_seqs.txt
        cat ${alignment_info} | awk '{ if (\$9 == 1) {print \$11}}' > forward_seqs.txt
        seqtk subseq ${query_genome} inverted_seqs.txt > ${query}_rev.fa
        seqtk subseq ${query_genome} forward_seqs.txt > ${query}_fwd.fa
        seqtk seq -r ${query}_rev.fa > ${query}_rev_rc.fa
        cat ${query}_fwd.fa ${query}_rev_rc.fa > ${query}_oriented.fa
        """
}
process SEQTK_SUBSET {
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
        tuple val(meta), path(genome)

    output:
        tuple val(meta), path("*_subset.fa"), emit: subset

    script:
    def pattern = params.subset_pattern
        """
        if [[ ${genome} == *.gz ]];
            then
                gzip -dc ${genome} > ${meta}_genome.fa
        else
            ln -s ${genome} ${meta}_genome.fa
        fi
        grep $pattern ${meta}_genome.fa | sed 's/>//' > names.lst
        seqtk subseq ${meta}_genome.fa names.lst > ${meta}_subset.fa
        """
}