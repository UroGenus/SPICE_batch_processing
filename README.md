# SPICE batch processing

This repository contains scripts for batch processing of raw .bam files using the [SPICE pipeline](https://github.com/demichelislab/SPICE-pipeline-CWL). SPICE is designed to generate a comprehensive view of the genomic landscape of matched tumor and normal samples by leveraging allele-specific information from high quality next generation sequencing data. Specifically, using targeted DNA data (e.g. WES), the pipeline first applies tumor purity and ploidy correction, generates a quantitative measure of aneuploidy (asP), and then calls a series of genomic aberrations, including allele-specific copy number aberrations, SNVs with copy number corrected allelic fractions and indels.

### Installation on UBELIX
```
git clone https://github.com/demichelislab/SPICE-pipeline-CWL
conda install -c conda-forge cwltool
```

### Process one file 
```
conda activate spice-env
module load nodejs/14.17.0-GCCcore-10.3.0
cwltool --no-container PATH_TO/SPICE-pipeline-CWL/cwl/workflows/pipeline.cwl PATH_TO/parameters.yaml

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
  path: /storage/research/dbmr_urology/Prostate_PDO/WG_IAD127899.20170720.designed.bed
kit_bait_bed_file:
  class: File
  path: /storage/research/dbmr_urology/Prostate_PDO/WG_IAD127899.20170720.designed.bed
kit_target_interval_file:
  class: File
  path: /storage/research/dbmr_urology/Prostate_PDO/SPICE_data/WG_IAD127899.20170720.interval_list
kit_bait_interval_file:
  class: File
  path: /storage/research/dbmr_urology/Prostate_PDO/SPICE_data/WG_IAD127899.20170720.interval_list
snps_in_kit_vcf_file:
  class: File
  path: /storage/research/dbmr_urology/Prostate_PDO/SPICE_data/WG_IAD127899.20170720.w.chr.snp
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
  path: /storage/research/dbmr_urology/Prostate_PDO/SPICE_data/homo_sapiens_vep_104_GRCh37
threads: 5
create_reports: true
log_to_file: true
```

### Parameters for WES (one example file)
```
TBC
```

### Obtaining files required for processing

#### kit_target_interval_file, kit_bait_interval_file 

##### IonTorrent
```
cd DIR_WITH_PipeIT.img
singularity exec --cleanenv -B /path/to/current/directory/ PipeIT.img cp /PipeIT_resources/usr/local/res/hg19.fasta .
wget https://github.com/broadinstitute/picard/releases/download/2.26.5/picard.jar
module load Java/11.0.2
java -jar picard.jar CreateSequenceDictionary R=/storage/research/dbmr_urology/Prostate_PDO/hg19.fa O=/storage/research/dbmr_urology/Prostate_PDO/hg19.dict
java -jar picard.jar BedToIntervalList  I=/storage/research/dbmr_urology/Prostate_PDO/WG_IAD127899.20170720.designed.bed O=/storage/research/dbmr_urology/Prostate_PDO/WG_IAD127899.20170720.interval_list SD=/storage/research/dbmr_urology/Prostate_PDO/hg19.fa
module load vital-it/7
module load UHTS/Analysis/samtools/1.10
samtools faidx hg19.fa
cp hg19.fa.fai hg19.fai
cp hg19.dict hg19.fa.dict
```

##### WES
```
???
```

#### snps_in_kit_vcf_file

One needs to install GATK. Failed to do it using [UBELIX EasyBuild instructions](https://hpc-unibe-ch.github.io/software/EasyBuild.html), it was then installed locally by UBELIX admins. After the installation, do 

```
download Homo_sapiens_assembly19.dbsnp.vcf, Homo_sapiens_assembly19.dict, Homo_sapiens_assembly19.fasta, Homo_sapiens_assembly19.fasta.fai from https://console.cloud.google.com/storage/browser/gcp-public-data--broad-references/hg19/v0
module load Workspace/home GATK
gatk IndexFeatureFile -I Homo_sapiens_assembly19.dbsnp.vcf
gatk SelectVariants -R /storage/research/dbmr_urology/Prostate_PDO/SPICE_data/Homo_sapiens_assembly19.fasta -V Homo_sapiens_assembly19.dbsnp.vcf -L /storage/research/dbmr_urology/Prostate_PDO/SPICE_data/WG_IAD127899.20170720.wout.chr.interval_list -select-type-to-include SNP -O WG_IAD127899.20170720.snp

```

#### CURRENT ISSUES

SPICE pipeline breaks on the following command
```
'gatk4' 'CollectHsMetrics' '--INPUT' '/tmp/xg6t3p0f/stga4de0ed6-3df9-4950-87e4-db415f469bba/IonXpress-025.bam' '--BAIT_INTERVALS' '/tmp/xg6t3p0f/stg1b18cac0-af3f-4143-8bab-ad2519a13013/WG_IAD127899.20170720.interval_list' '--TARGET_INTERVALS' '/tmp/xg6t3p0f/stg1b18cac0-af3f-4143-8bab-ad2519a13013/WG_IAD127899.20170720.interval_list' '--VALIDATION_STRINGENCY' 'LENIENT' '--OUTPUT' 'hsmetrics_tumor/hsmetrics_tumor.txt' '--PER_TARGET_COVERAGE' 'hsmetrics_tumor/hsmetrics_per_target_coverage_tumor.txt' '--REFERENCE_SEQUENCE' '/tmp/xg6t3p0f/stgcf87b360-d190-47d7-af60-77185726fc13/hg19.fa' 2> picard.log 1>&2
```

Tried to call the command explicitly like this:
```
gatk CollectHsMetrics --INPUT /storage/research/dbmr_urology/Prostate_PDX/IonTorrent/bam/IonXpress-026.bam --BAIT_INTERVALS /storage/research/dbmr_urology/Prostate_PDO/SPICE_data/WG_IAD127899.20170720.interval_list --TARGET_INTERVALS /storage/research/dbmr_urology/Prostate_PDO/SPICE_data/WG_IAD127899.20170720.interval_list --VALIDATION_STRINGENCY LENIENT --OUTPUT tmp/hsmetrics_tumor.txt --PER_TARGET_COVERAGE tmp/hsmetrics_per_target_coverage_tumor.txt --REFERENCE_SEQUENCE /storage/research/dbmr_urology/Prostate_PDO/hg19.fa
```

**Error message**: *htsjdk.samtools.util.SequenceUtil$SequenceListsDifferException: Sequence dictionaries are not the same size (25, 93)* 

**Hypothesis**: Alignment was done using different hg19.fasta file as compared to the one downloaded (https://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/hg19.fa.gz). Now trying to find out, which fasta file has been originally used.


