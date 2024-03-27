## Introduction
This pipeline performs the alignment of sequencing reads (after light pre-processing) and variant calling.

## How to Run and Requirements
- All code is in the wrapper `/bin/code_challenge.sh`. All tools and dependencies need to be installed locally.
- To run `/bin/code_challenge.sh`, you need to provide the location of fastq files (e.g., `/fastq`) 

```bash
  bin/code_challenge.sh your_input_folder
 ```
- The genome fasta file needs to be in `/annotation/MN908947.3.fasta`.
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

 * [your_input_folder]
   * [*.R1.paired.fq.gz]
   * [*.R1.paired.fq.gz]
 * [annotation]
   * [MN908947.3.fasta]
 * [adapters]
   * [NexteraPE.fa]
 * [bin]
   * [code_challenge.sh]

## Results and interpretations
- The quality of the sequencing is ok, no major drops in reads numbers durign pre-processing. 
- The quality of variants and depth are also ok (even before filtering): 
```bash
bcftools query -f '%DP\n' vcf_results/sample/var.vcf | awk '{sum += $1} END {print "Average depth of coverage:", sum/NR}'
 ```
Average depth of coverage: 44.0816
```bash
bcftools query -f '%QUAL\n' vcf_results/sample/var.vcf | awk '{ sum += $1 } END { if (NR > 0) print "Average QUAL score:", sum / NR }'
 ```
Average QUAL score: 1184.97
```bash
bcftools filter -i 'QUAL < 20' vcf_results/sample/var.vcf | bcftools view -H | wc -l
 ```
Only 8 low quality sites, easy to filter out
```bash
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' vcf_results/sample/var.vcf | awk '$5 > 0.95 {print $1, $2, $3, $4, $5}'
 ```
- There are 40 positions with AF=1 (>0.95): it could be a contamination or population-specific variant
- To better assess contaminations in general, I'd include Fastq Screen that would align the sequences against a panel of different genomes.
- There are 6 variants with AF < 0.05: may be considered potential mosaic candidates.
- All types of variants: mosaic or potential contaminations should be assessed and validated further, for example by amplicon-sequencing with UMIs to account for PCR duplications and confidently detect rare variants

