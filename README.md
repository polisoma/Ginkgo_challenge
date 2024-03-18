## Introduction
This pipeline performs the alignment of sequencing reads (after light pre-processing) and variant calling.

## How to Run and Requirements
- All code is in the wrapper `/bin/code_challenge.sh`. All tools and dependencies need to be installed locally.
- To run `/bin/code_challenge.sh`, you need to provide the location of fastq files (e.g., `/fastq`) and the genome fasta file (`/annotation`).
- For read trimming, you need to have `/adapters/NexteraPE.fa`.
- The wrapper will create all necessary folders and files (and will delete heavy intermediate files, like `.sam`).
- All errors will be logged in `/logs`.
- The final files will be stored in `/vcf_results`: `var_filtered.vcf` (filtered for QUAL and DP) and `sample.html` (the annotation of variants).
- Tools required locally:
  - bwa 0.7.17 (and `/annotation/MN908947.3.fasta`; it will create the index if it doesn't exist)
  - samtools 1.19.2
  - trimmomatic 0.39 (and `/adapters/NexteraPE.fa`)
  - freebayes 1.3.6
  - picard 3.1.1
  - bcftools 1.19
  - snpEff 5.2 (it will download NC_045512.2)
- The number of cores is not specified (for samtools and so on), but should be included for scaling up.
- Structure of the folders and required files to run the wrapper:
.
 * [your_input_folder](./your_input_folder)
   * [*.R1.paired.fq.gz](./your_input_folder/*.R1.paired.fq.gz)
   * [*.R1.paired.fq.gz](./your_input_folder/*.R1.paired.fq.gz)
 * [annotation](./annotation)
   * [MN908947.3.fasta](./annotation/MN908947.3.fasta)
 * [adapters](./adapters)
   * [NexteraPE.fa](./adapters/NexteraPE.fa)
 * [bin](./bin)
   * [code_challenge.sh](./bin/code_challenge.sh)


## Results and interpretations

