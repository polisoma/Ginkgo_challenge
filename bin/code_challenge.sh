#!/bin/bash

# Function to display help
display_help() {
    echo "Usage: $0 input_folder"
    echo "Example: $0 /path/to/input"
    echo
    echo "Required Tools:"
    echo "- bwa 0.7.17 (and /annotation/MN908947.3.fasta; it will create the index if it doesn't exist)"
    echo "- samtools 1.19.2"
    echo "- trimmomatic 0.39 (and /adapters/NexteraPE.fa)"
    echo "- freebayes 1.3.6"
    echo "- picard 3.1.1"
    echo "- bcftools 1.19"
    echo "- snpEff 5.2 (it will download NC_045512.2)"
    exit 1
}

# Check for correct number of arguments
if [ "$#" -ne 1 ]; then
    display_help
fi

# Save the input folder path as a variable
INPUT_folder=$1

# Define the reference filename without the .fasta extension
REFERENCE="MN908947.3"

# Define the path to the adapters file
ADAPTERS="adapters/NexteraPE.fa"

# folder for logs
mkdir -p logs

# Check if the input folder exists and contains the required files
if [ ! -d "$INPUT_folder" ]; then
    echo "Error: Input folder does not exist or is not a directory."
    exit 1
fi

if ! ls "${INPUT_folder}"/*R1.paired.fq.gz > /dev/null 2>&1; then
    echo "Error: No files found in input folder matching the pattern *R1.paired.fq.gz."
    exit 1
fi

# Create an array to store the names of the sample files
sample_names=()

# Read the names of the sample files from the fastq folder
for file in "$INPUT_folder"/*R1.paired.fq.gz; do
    # Extract the base name of the file without the path and suffix
    sample_name=$(basename "$file" .R1.paired.fq.gz)
    # Add the sample name to the array
    sample_names+=("$sample_name")
done

# Check if the reference index exists, otherwise create it
file="annotation/${REFERENCE}.fasta.bwt"
if [ ! -e "$file" ]; then
    echo "File $file does not exist. Making bwa index..."
    bwa index "annotation/${REFERENCE}.fasta" "annotation/${REFERENCE}"
    if [ $? -ne 0 ]; then
        echo "Error: Failed to create bwa index."
        exit 1
    fi
fi

# Run fastQC to assess the quality before mapping
mkdir -p fastqc
fastqc -o fastqc/ "$INPUT_folder"/*R1.paired.fq.gz "$INPUT_folder"/*R2.paired.fq.gz
if [ $? -ne 0 ]; then
    echo "Error: Failed to run fastQC."
    exit 1
fi

# Loop over each sample
for sample_name in "${sample_names[@]}"; do
    echo "Processing sample: $sample_name"
    # quality is good --> light trimming and adaptor trimming
    mkdir -p trimming

    # Run trimmomatic for each sample
    trimmomatic PE -phred33 -trimlog "logs/${sample_name}_trimming.log" \
    "$INPUT_folder/${sample_name}.R1.paired.fq.gz" "$INPUT_folder/${sample_name}.R2.paired.fq.gz" \
    "trimming/${sample_name}_R1.trimmed.fq.gz" "trimming/${sample_name}_R1.unpaired.fq.gz" \
    "trimming/${sample_name}_R2.trimmed.fq.gz" "trimming/${sample_name}_R2.unpaired.fq.gz" \
    "ILLUMINACLIP:$ADAPTERS:2:30:10" \
    "HEADCROP:10" "LEADING:20" "TRAILING:20" "SLIDINGWINDOW:4:20" "MINLEN:36"
    if [ $? -ne 0 ]; then
        echo "Error: Failed to run trimmomatic for sample $sample_name."
        exit 1
    fi

        # Mapping to the genome with bwa
    mkdir -p "alignment/$sample_name"
    bwa mem -M \
    "annotation/${REFERENCE}.fasta" \
    "trimming/${sample_name}_R1.trimmed.fq.gz" "trimming/${sample_name}_R2.trimmed.fq.gz" \
    2> "logs/${sample_name}_bwa.err" \
    > "alignment/$sample_name/${sample_name}.sam"
    if [ $? -ne 0 ]; then
        echo "Error: Failed to map reads for sample $sample_name."
        exit 1
    fi

    # Samtools
    samtools view -SbhF 4 -q 10 "alignment/$sample_name/${sample_name}.sam" | samtools sort > "alignment/$sample_name/${sample_name}.bam"
    if [ $? -ne 0 ]; then
        echo "Error: Failed to convert SAM to BAM for sample $sample_name."
        exit 1
    fi

    # Mark duplicates and make sample_stat with some quality and sanity checks
    picard MarkDuplicates \
    INPUT="alignment/$sample_name/${sample_name}.bam" \
    OUTPUT="alignment/$sample_name/${sample_name}_rmdup.bam" \
    METRICS_FILE="alignment/$sample_name/${sample_name}_stat" \
    REMOVE_DUPLICATES=True
    if [ $? -ne 0 ]; then
        echo "Error: Failed to mark duplicates for sample $sample_name."
        exit 1
    fi

    # Index BAM file
    samtools index "alignment/$sample_name/${sample_name}_rmdup.bam"
    if [ $? -ne 0 ]; then
        echo "Error: Failed to index BAM file for sample $sample_name."
        exit 1
    fi

    # Remove intermediate files including sam (too heavy for storage)
    rm "alignment/$sample_name/${sample_name}.sam" "alignment/$sample_name/${sample_name}.bam"
    if [ $? -ne 0 ]; then
        echo "Error: Failed to remove intermediate files for sample $sample_name."
        exit 1
    fi

    # Calling variants without any filter criteria (will do filtering on vcf file)
    # Compress vcf file since they are usually heavy
    mkdir -p "vcf_results/$sample_name"
    freebayes -f "annotation/${REFERENCE}.fasta" "alignment/$sample_name/${sample_name}_rmdup.bam" > "vcf_results/$sample_name/var.vcf"
    if [ $? -ne 0 ]; then
        echo "Error: Failed to call variants for sample $sample_name."
        exit 1
    fi

    # Filtering variants (QUAL and DP for specific sample; we have only one)
    bcftools=../software/bcftools/bcftools
    $bcftools filter -o "vcf_results/$sample_name/var_filtered.vcf" -i 'QUAL>20 && FORMAT/DP[0]>10' "vcf_results/$sample_name/var.vcf"
    if [ $? -ne 0 ]; then
        echo "Error: Failed to filter variants for sample $sample_name."
        exit 1
    fi

    # Variant annotation (the chromosome name will be not the same as in fasta)
    snpEff download NC_045512.2
    # Making the chromosome names the same in MN908947.3 and NC_045512.2
    perl -pi -e 's/MN908947.3/NC_045512.2/g' "vcf_results/$sample_name/var_filtered.vcf"
    snpEff -v -stats "vcf_results/$sample_name/${sample_name}.html" NC_045512.2 "vcf_results/$sample_name/var_filtered.vcf" > "vcf_results/$sample_name/${sample_name}.ann.vcf"
    if [ $? -ne 0 ]; then
        echo "Error: Failed to annotate variants for sample $sample_name."
        exit 1
    fi
done


