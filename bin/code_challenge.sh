#!/bin/bash
# required tools to be installed
# fasta should be in the annotation folder

set -euo pipefail

# options
INPUT_R1=$1
INPUT_R2=$1
INPUT_folder=$1
OUTPUT_folder=$2
REFERENCE_folder=annotation #$3

# the folder structure:
mkdir -p logs

# creating the bwa index (if doesn't exist in the annotation folder)
file=$REFERENCE_folder/MN908947.3.fasta.bwt
if [ -e "$file" ]; then
    echo "File $file exists -> OK"
else
    echo "File $file does not exist. Making bwa index..."
    bwa index annotation/MN908947.3.fasta annotation/MN908947.3
fi

# runnign fastQC to assess the quality before mapping
mkdir -p fastqc
fastqc -o fastqc/ fastq/*.fq.gz

# quality is good --> light trimming and adaptor trimming
mkdir -p trimming 
# adapters were guessed based in the infor from SRA: NEB_ARTICv1 design
adapters=adapters/NexteraPE.fa

trimmomatic PE -phred33 -trimlog logs/trimming.log \
fastq/sample.R1.paired.fq.gz fastq/sample.R2.paired.fq.gz \
trimming/sample.R1.trimmed.fq.gz trimming/sample.R1.unpaired.fq.gz \
trimming/sample.R2.trimmed.fq.gz trimming/sample.R2.unpaired.fq.gz \
ILLUMINACLIP:$adapters:2:30:10 \
HEADCROP:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:36

# mapping to the genome with bwa
# -M 
mkdir -p alignment
# bwa -t parameter is ommited, but would be useful to specify the number of cores
bwa mem -M \
$REFERENCE_folder/MN908947.3.fasta \
trimming/sample.R1.trimmed.fq.gz trimming/sample.R2.trimmed.fq.gz \
2> logs/bwa.err \
> alignment/sample.sam

# samtools -@ ommited, but good to specify the number of cores
# filtering out unmapped reads, including only reads with a quality of 10 or higher
samtools view -SbhF 4 -q 10 alignment/sample.sam | samtools sort > alignment/sample.bam
	
# marking duplicates and making sample_stat with some quality and sanity checks
picard MarkDuplicates \
INPUT=alignment/sample.bam \
OUTPUT=alignment/sample.rmdup.bam \
METRICS_FILE=alignment/sample_stat \ 
REMOVE_DUPLICATES=True

# making bam index
samtools index alignment/sample.rmdup.bam
# removing intermediate files incuding sam (too heavy for storage)
rm alignment/sample.sam alignment/sample.bam

# calling variants without any filter criteria (will do filtering on vcf file)
# compress vcf file since they are usually heavy
mkdir -p vcf_results
freebayes -f  $REFERENCE_folder/MN908947.3.fasta alignment/sample.rmdup.bam > vcf_results/var.vcf

# filtering variants (QUAL and DP for specific sample; we have only one)
bcftools=../software/bcftools/bcftools
$bcftools filter -o vcf_results/var_filtered.vcf -i 'QUAL>20 && FORMAT/DP[0]>10' vcf_results/var.vcf

# variant annotation (the chomosome name will be not the same as in fasta)
snpEff download NC_045512.2
# making the chromosome names the same in MN908947.3 and NC_045512.2
perl -pi -e 's/MN908947.3/NC_045512.2/g' vcf_results/var_filtered.vcf
snpEff -v -stats vcf_results/sample.html NC_045512.2 vcf_results/var_filtered.vcf > vcf_results/sample.ann.vcf


