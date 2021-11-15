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
cwltool --singularity PATH_TO/SPICE-pipeline-CWL/cwl/workflows/pipeline.cwl PATH_TO/parameters.yaml

```

### Parameters for IonTorrent (one example file)
```
bam_file_normal:
  class: File
  path: /storage/research/dbmr_urology/Prostate_PDX/IonTorrent/IonXpress-026.bam
bam_file_tumor:
  class: File
  path: /storage/research/dbmr_urology/Prostate_PDX/IonTorrent/IonXpress-025.bam
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
  path: path/to/snps_in_kit.vcf
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
wget https://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/hg19.fa.gz
gunzip hg19.fa.gz
java -jar picard.jar CreateSequenceDictionary R=/storage/research/dbmr_urology/Prostate_PDO/hg19.fa O=/storage/research/dbmr_urology/Prostate_PDO/hg19.dict
java -jar picard.jar BedToIntervalList  I=/storage/research/dbmr_urology/Prostate_PDO/WG_IAD127899.20170720.designed.bed O=/storage/research/dbmr_urology/Prostate_PDO/WG_IAD127899.20170720.interval_list SD=/storage/research/dbmr_urology/Prostate_PDO/hg19.fa
```

##### WES
```
???
```

