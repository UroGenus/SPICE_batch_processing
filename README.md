# SPICE batch processing

This repository contains scripts for batch processing of raw .bam files using the [SPICE pipeline](https://github.com/demichelislab/SPICE-pipeline-CWL). SPICE is designed to generate a comprehensive view of the genomic landscape of matched tumor and normal samples by leveraging allele-specific information from high quality next generation sequencing data. Specifically, using targeted DNA data (e.g. WES), the pipeline first applies tumor purity and ploidy correction, generates a quantitative measure of aneuploidy (asP), and then calls a series of genomic aberrations, including allele-specific copy number aberrations, SNVs with copy number corrected allelic fractions and indels.

### Installation on UBELIX
```
git clone https://github.com/demichelislab/SPICE-pipeline-CWL
conda install -c conda-forge cwltool
module load Workspace_Home
mkdir $SCRATCH
chmod 700 $SCRATCH
```

### Process one file 
```
conda activate spice-env
module load nodejs/14.17.0-GCCcore-10.3.0
cwltool --singularity PATH_TO/SPICE-pipeline-CWL/cwl/workflows/pipeline.cwl PATH_TO/parameters.yaml

```

### Parameters for IonTorrent (one example file)
```
bam_file_normal:
  class: File
  path: /storage/research/dbmr_urology/Prostate_PDX/IonTorrent/bam/IonXpress-026.bam
bam_file_tumor:
  class: File
  path: /storage/research/dbmr_urology/Prostate_PDX/IonTorrent/bam/IonXpress-025.bam
reference_genome_fasta_file:
  class: File
  path: /storage/research/dbmr_urology/Prostate_PDO/hg19.fa
kit_target_bed_file:
  class: File
  path: /storage/research/dbmr_urology/Prostate_PDO/WG_IAD127899.20170720.designed.sorted.merged.bed
kit_bait_bed_file:
  class: File
  path: /storage/research/dbmr_urology/Prostate_PDO/WG_IAD127899.20170720.designed.sorted.merged.bed
kit_target_interval_file:
  class: File
  path: /storage/research/dbmr_urology/Prostate_PDO/SPICE_data/WG_IAD127899.20170720.interval_list
kit_bait_interval_file:
  class: File
  path: /storage/research/dbmr_urology/Prostate_PDO/SPICE_data/WG_IAD127899.20170720.interval_list
snps_in_kit_vcf_file:
  class: File
  path: /storage/research/dbmr_urology/Prostate_PDO/SPICE_data/WG_IAD127899.20170720.snp.vcf
ethseq_snps_vcf_file:
  class: File
  path: /storage/research/dbmr_urology/Prostate_PDO/SPICE_data/ethseq-universal_exonic_model-hg19.vcf
ethseq_snps_gds_file:
  class: File
  path: /storage/research/dbmr_urology/Prostate_PDO/SPICE_data/ethseq-universal_exonic_model-hg19.gds
spia_snps_vcf_file:
  class: File
  path: /storage/research/dbmr_urology/Prostate_PDO/SPICE_data/snps_spia_default-hg19.vcf
sample_sex: m
vep_reference_genome_version: GRCh37
vep_data_directory:
  class: Directory
  path: /storage/research/dbmr_urology/Prostate_PDO/SPICE_data/homo_sapiens
threads: 5
create_reports: true
log_to_file: true
```

### Parameters for WES (one example file)
```
bam_file_normal:
  class: File
  path: /storage/research/dbmr_urology/SequencingData/PNPCa_seq_cornell/cornell_wes/PM1548-SK1N1.bam
bam_file_tumor:
  class: File
  path: /storage/research/dbmr_urology/SequencingData/PNPCa_seq_cornell/cornell_wes/PM1548-SK1T1.bam
reference_genome_fasta_file:
  class: File
  path: /storage/research/dbmr_urology/Prostate_PDO/SPICE_data/hg19_WES.fa
kit_target_bed_file:
  class: File
  path: /storage/research/dbmr_urology/Prostate_PDO/SPICE_data/Haloplex_Regions_b37_Covered_noHeader_with_genenames.bed
kit_bait_bed_file:
  class: File
  path: /storage/research/dbmr_urology/Prostate_PDO/SPICE_data/Haloplex_Regions_b37_Covered_noHeader_with_genenames.bed
kit_target_interval_file:
  class: File
  path: /storage/research/dbmr_urology/Prostate_PDO/SPICE_data/Haloplex_Regions_b37_Covered_noHeader_with_genenames.interval_list
kit_bait_interval_file:
  class: File
  path: /storage/research/dbmr_urology/Prostate_PDO/SPICE_data/Haloplex_Regions_b37_Covered_noHeader_with_genenames.interval_list
snps_in_kit_vcf_file:
  class: File
  path: /storage/research/dbmr_urology/Prostate_PDO/SPICE_data/Haloplex_Regions_b37_Covered_noHeader_with_genenames.snp.vcf
ethseq_snps_vcf_file:
  class: File
  path: /storage/research/dbmr_urology/Prostate_PDO/SPICE_data/ethseq-universal_exonic_model-hg19_nochr.vcf
ethseq_snps_gds_file:
  class: File
  path: /storage/research/dbmr_urology/Prostate_PDO/SPICE_data/ethseq-universal_exonic_model-hg19_nochr.gds
spia_snps_vcf_file:
  class: File
  path: /storage/research/dbmr_urology/Prostate_PDO/SPICE_data/snps_spia_default-hg19_nochr.vcf
sample_sex: m
vep_reference_genome_version: GRCh37
vep_data_directory:
  class: Directory
  path: /storage/research/dbmr_urology/Prostate_PDO/SPICE_data/homo_sapiens
threads: 5
create_reports: true
log_to_file: true
```

### Obtaining files required for processing

#### kit_target_bed_file, kit_bait_bed_file

##### IonTorrent

```
sort -k 1,1 -k2,2n WG_IAD127899.20170720.designed.bed > WG_IAD127899.20170720.designed.sorted.bed
module load UHTS/Analysis/BEDTools/2.29.2
bedtools merge -c 4,5,6 -o distinct,distinct,distinct -i WG_IAD127899.20170720.designed.sorted.bed > WG_IAD127899.20170720.designed.sorted.merged.bed
```

#### reference_genome_fasta_file

##### IonTorrent
```
cd DIR_WITH_PipeIT.img
singularity exec --cleanenv -B /path/to/current/directory/ PipeIT.img cp /PipeIT_resources/usr/local/res/hg19.fasta .
module load vital-it/7
module load UHTS/Analysis/samtools/1.10
samtools faidx hg19_IonTorrent.fa
cp hg19_IonTorrent.fa.fai hg19_IonTorrent.fai
cp hg19_IonTorrent.dict hg19_IonTorrent.fa.dict
```

#### kit_target_interval_file, kit_bait_interval_file 

##### IonTorrent
```
wget https://github.com/broadinstitute/picard/releases/download/2.26.5/picard.jar
module load Java/11.0.2
java -jar PATHTO/picard.jar CreateSequenceDictionary R=/storage/research/dbmr_urology/Prostate_PDO/hg19_IonTorrent.fa O=/storage/research/dbmr_urology/Prostate_PDO/hg19_IonTorrent.dict
java -jar PATHTO/picard.jar BedToIntervalList  I=/storage/research/dbmr_urology/Prostate_PDO/WG_IAD127899.20170720.designed.sorted.merged.bed O=/storage/research/dbmr_urology/Prostate_PDO/SPICE_data/WG_IAD127899.20170720.interval_list SD=/storage/research/dbmr_urology/Prostate_PDO/hg19_IonTorrent.fa
```

##### WES
```
awk -F'\t' -v OFS="\t" '{ print $2, $3, $4, $1 }' S06588914_Agilent_Clinical_Exome_TargetGeneSymbols.txt > S06588914_Agilent_Clinical_Exome_TargetGeneSymbols.bed
module load UHTS/Analysis/BEDTools/2.29.2
bedtools intersect -wa -wb -a Haloplex_Regions_b37_Covered_noHeader_noAnn.bed -b S06588914_Agilent_Clinical_Exome_TargetGeneSymbols.bed > tmp.bed
awk -F'\t' -v OFS="\t" '{ print $1, $2, $3, $7 }' tmp.bed > Haloplex_Regions_b37_Covered_noHeader_with_genenames.bed
java -jar PATHTO/picard.jar BedToIntervalList I=/storage/research/dbmr_urology/Prostate_PDO/SPICE_data/Haloplex_Regions_b37_Covered_noHeader_with_genenames.bed O=/storage/research/dbmr_urology/Prostate_PDO/SPICE_data/Haloplex_Regions_b37_Covered_noHeader_with_genenames.interval_list SD=/storage/research/dbmr_urology/Prostate_PDO/SPICE_data/hg19_WES.fa
```

#### snps_in_kit_vcf_file

One needs to install GATK. Failed to do it using [UBELIX EasyBuild instructions](https://hpc-unibe-ch.github.io/software/EasyBuild.html), it was then installed locally by UBELIX admins. After the installation, do 

```
download Homo_sapiens_assembly19.dbsnp.vcf, Homo_sapiens_assembly19.dict, Homo_sapiens_assembly19.fasta, Homo_sapiens_assembly19.fasta.fai from https://console.cloud.google.com/storage/browser/gcp-public-data--broad-references/hg19/v0
module load Workspace_Home/1.1
module load Workspace/home GATK
gatk IndexFeatureFile -I Homo_sapiens_assembly19.dbsnp.vcf
```

##### IonTorrent
```
./SPICE.preprocessing.sh
cat WG_IAD127899.20170720.snp.vcf | gawk '{if($0 ~ /^#/){print $0}else{if($5 !~ /,/){print $0}}}' > WG_IAD127899.20170720.snp.wout.comma.vcf
```

##### WES
```
gatk SelectVariants -R /storage/research/dbmr_urology/Prostate_PDO/SPICE_data/Homo_sapiens_assembly19.fasta -V /storage/research/dbmr_urology/Prostate_PDO/SPICE_data/Homo_sapiens_assembly19.dbsnp.vcf -L Haloplex_Regions_b37_Covered_noHeader_with_genenames.interval_list -select-type-to-include SNP -O Haloplex_Regions_b37_Covered_noHeader_with_genenames.snp.vcf
```

#### CURRENT ISSUES

SPICE pipeline breaks on the following command
```
'ethseq' '--enable_plot' '--num_threads' '5' '/var/lib/cwl/stg0c4dda1c-37d5-42a9-9499-22d03456c9ac/ethseq-universal_exonic_model-hg19.gds' '/var/lib/cwl/stg9646c714-3e33-4ff0-8dfb-47ca0577276b/genotype.vcf' 'ethseq'
[2021-12-10 09:19:14] Running EthSEQ
[2021-12-10 09:19:14] Working directory: ethseq
[2021-12-10 09:19:14] Create ethseq folder
[2021-12-10 09:19:14] Create target model from VCF
[2021-12-10 09:19:19] Create aggregated model
[2021-12-10 09:19:21] ERROR: Target and reference models are not compatible.
```
Hypothesis: [ethseq](https://github.com/cibiobcg/EthSEQ) default parameter is `bam.chr.encoding = FALSE`. We might need to change that.

To run ethseq on Ubelix:
```
R
library('EthSEQ')
## Run the analysis
ethseq.Analysis(
  bam.list = file.path(data.dir,"BAMs_List.txt"),
  out.dir = out.dir,
  model.gds = system.file("extdata","Reference_SS2_10000SNPs.gds",
     package="EthSEQ"),
  verbose=TRUE,
  aseq.path = out.dir,
  mbq=20,
  mrq=20,
  mdc=10,
  run.genotype = TRUE,
  composite.model.call.rate = 1,
  cores=1,
  bam.chr.encoding = FALSE) # chromosome names encoded without "chr" prefix in BAM files
```
