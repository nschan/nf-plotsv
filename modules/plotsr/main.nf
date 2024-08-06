process PLOTSR {
    tag "$meta"
    label 'process_low'
    publishDir(
      path: { "${params.out}/${task.process}".replace(':','/').toLowerCase() }, 
      mode: 'copy',
      overwrite: true,
      saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) }
    ) 
    input:
        tuple val(meta), path(syri_out)
        val(reference)
        path(plotsr_conf)

    output:
        tuple val(meta), path("*plotsr.pdf"), emit: figure

    def plotsr_conf = file("$projectDir/assets/plotsr_config.conf", checkIfExists: true)
    script:
        """
        micromamba run -n base plotsr \\
            ${meta}_on_${reference}.syri.out $reference $meta \\
            -H 8 -W 5 \\
            --cfg $plotsr_conf \\
            -o ${meta}_on_${reference}.plotsr.pdf 
        """
}
process PLOTSR_PAIRWISE {
    tag "BIGPLOT"
    label 'process_low'
    publishDir(
      path: { "${params.out}/${task.process}".replace(':','/').toLowerCase() }, 
      mode: 'copy',
      overwrite: true,
      saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) }
    ) 
    input:
        path(in_files)
        val(names)
        val(genomes)
        path(plotsr_conf)

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

    paste genomes.col names.col >> plotsr_infile.tsv

    micromamba run -n base plotsr --genomes plotsr_infile.tsv \\
    \$files \\
    --cfg ${plotsr_conf} \\
    -o plot.pdf
    """
} 