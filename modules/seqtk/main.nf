include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process SEQTK_ORIENT {
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