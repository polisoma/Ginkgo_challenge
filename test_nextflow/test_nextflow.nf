params.input_folder = '/path/to/input'

process checkInputFolder {
    output:
    file 'sample_names.txt' into sampleNamesCh

    script:
    """
    if [ ! -d "${params.input_folder}" ]; then
        echo "Error: Input folder does not exist or is not a directory."
        exit 1
    fi

    if ! ls "\${params.input_folder}"/*R1.paired.fq.gz > /dev/null 2>&1; then
        echo "Error: No files found in input folder matching the pattern *R1.paired.fq.gz."
        exit 1
    fi

    ls "\${params.input_folder}"/*R1.paired.fq.gz | xargs -n 1 basename | sed 's/_R1.paired.fq.gz//' > sample_names.txt
    """
}

process bwaIndex {
    input:
    file referenceFasta from 'annotation/MN908947.3.fasta'

    script:
    """
    bwa index annotation/MN908947.3.fasta
    """
}

process fastqc {
    input:
    file 'sample_names.txt' from sampleNamesCh
    set val(sample_name)

    script:
    """
    mkdir -p fastqc
    fastqc -o fastqc/ "\${params.input_folder}/\${sample_name}_R1.paired.fq.gz" "\${params.input_folder}/\${sample_name}_R2.paired.fq.gz"
    """
}



