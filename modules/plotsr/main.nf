include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process PLOTSR {
    tag "$meta"
    label 'process_low'
    publishDir "${params.out}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename,
                                        options:params.options, 
                                        publish_dir:"${task.process}".replace(':','/').toLowerCase(), 
                                        publish_id:meta) }
    input:
        tuple val(meta), path(syri_out)
        val(reference)

    output:
        tuple val(meta), path("*plotsr.pdf"), emit: figure

    script:
        """
        micromamba run -n base plotsr ${meta}_on_${reference}.syri.out $reference $meta -H 8 -W 5 -o ${meta}_on_${reference}.plotsr.pdf
        """
}
process PLOTSR_PAIRWISE {
    tag "BIGPLOT"
    label 'process_low'
    publishDir "${params.out}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename,
                                        options:params.options, 
                                        publish_dir:"${task.process}".replace(':','/').toLowerCase(), 
                                        publish_id:"plot") }
    input:
        path(in_files)
        val(names)
        val(genomes)

    output:
        path("*.pdf"), emit: figure

    script:
    """
    files_array=( ${in_files} )
    files="\${files_array[@]/#/--sr }"
    echo \$files > files.txt

    names_array=(${names})
    echo \$names_array > names.txt

    genomes_array=(${genomes})
    echo \$genomes_array > genomes.txt

    sed -i 's/\\[//g' genomes.txt 
    sed -i 's/\\]//g' genomes.txt 
    sed -i 's/,//g' genomes.txt 
    for x in `cat genomes.txt`
    do
        echo \$x
    done > genomes.col

    sed -i 's/\\[//g' names.txt 
    sed -i 's/\\]//g' names.txt 
    sed -i 's/,//g' names.txt 
    for x in `cat names.txt`
    do
        echo \$x
    done > names.col

    paste names.col genomes.col >> plotsr_infile.tsv

    micromamba run -n base plotsr --genomes plotsr_infile.tsv \\
    \$files \\
    -o plot.pdf
    """
} 