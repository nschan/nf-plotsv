include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process FIXCHR {
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
        tuple val(reference), path(reference_genome), val(query), path(query_genome), path(alignment)

    output:
        tuple val(reference), val(query), path(query_genome), path("*input_alignments.txt"), emit: alignment_info

    script:
        """
        micromamba run -n base fixchr -c ${alignment} -F B -r ${reference_genome} -q ${query_genome} --prefix ${query} --log INFO --contig_size 5000 
        """
}