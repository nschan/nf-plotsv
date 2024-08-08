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
        val reference 
        path plotsr_conf 
        path extra_args
    output:
        tuple val(meta), path("*plotsr.pdf"), emit: figure

    def plotsr_conf = file("$projectDir/assets/plotsr_config.conf", checkIfExists: true)
    def plotsr_args = extra_args ?: ''
    script:
        """
        micromamba run -n base plotsr \\
            ${meta}_on_${reference}.syri.out $reference $meta \\
            -H 8 -W 5 \\
            --cfg $plotsr_conf \\
            $plotsr_args \\
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
        path in_files
        val names 
        val genomes
        path plotsr_conf 
        val extra_args
        val tracks
        val palette

    output:
        path("*.pdf"), emit: figure

    script:

    def plotsr_tracks = tracks.equals('') || tracks == null ? '' : "--tracks ${tracks}"
    def plotsr_palette = palette.equals('') || palette == null ? '\\#8F7C00' : "${palette.toString().replace("#","\\#")}"

    """
    files_array=( ${in_files} )
    files="\${files_array[@]/#/--sr }"
    echo \$files > files.txt

    names_array=( ${names} )
    echo \$names_array > names.txt
    sed -i 's/\\[//g' names.txt 
    sed -i 's/\\]//g' names.txt 
    sed -i 's/,//g' names.txt 
    for x in `cat names.txt`
    do
        echo \$x
    done > names.col
    len_names=\$(cat names.col | wc -l)

    genomes_array=( ${genomes} )
    echo \$genomes_array > genomes.txt
    sed -i 's/\\[//g' genomes.txt 
    sed -i 's/\\]//g' genomes.txt 
    sed -i 's/,//g' genomes.txt 
    for x in `cat genomes.txt`
    do
        echo \$x
    done > genomes.col
    len_genomes=\$(cat genomes.col | wc -l)

    color_array=( ${plotsr_palette} )
    len_colors=\${#color_array[@]}

    # Here is a bunch of crap to make sure that colors are the same length as the samplesheet

    ## If they are equal its nice 
    if [ \$len_colors -eq \$len_names ]; then
        color_array2=( "\${color_array[@]}" )
    
    ## If there are more colors than genomes, subset colors
    elif [ \$len_colors -gt \$len_names ]; then
        color_array2=( "\${color_array[@]:0:\$((\$len_names-1))}" )

    ## If there are less colors than genomes, repeat them
    elif [ \$len_colors -lt \$len_names ]; then
        len_fac=\$((\$len_names / \$len_colors)) # Take full divisions
        len_mod=\$((\$len_names % \$len_colors)) # Take mod
        for i in \$(seq 1 \$len_fac); do color_array2+=(\${ color_array[@]}); done
        color_array2+=(\${color_array[@]:0:\$((\$len_mod))}) 
    fi

    # prefix the color hexcode with lc:
    colors="\${color_array2[@]/#/lc:}"
    echo \$colors > colors.txt
    
    # Turn txt into col
    for x in `cat colors.txt`
    do
        echo \$x
    done > colors.col

    paste genomes.col names.col colors.col >> plotsr_infile.tsv

    micromamba run -n base plotsr --genomes plotsr_infile.tsv \\
    \$files \\
    --cfg ${plotsr_conf} \\
    $extra_args \\
    ${plotsr_tracks} \\
    -o plot.pdf
    """
} 