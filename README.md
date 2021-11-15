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

#### snps_in_kit_vcf_file

One needs to install GATK. Failed to do it using [UBELIX EasyBuild instructions](https://hpc-unibe-ch.github.io/software/EasyBuild.html), it was installed by UBELIX admins. After the installation, run 

```
module load Workspace/home GATK
gatk SelectVariants -R /storage/research/dbmr_urology/Prostate_PDO/hg19.fa -V Homo_sapiens_assembly19.dbsnp.vcf -L WG_IAD127899.20170720.interval_list -select-type-to-include SNP -O WG_IAD127899.20170720.snp

```

**ERROR MESSAGE BECAUSE OF CHROMOSOME FORMAT MISMATCH**: 
*A USER ERROR has occurred: Badly formed genome unclippedLoc: Contig chr1 given as location, but this contig isn't present in the Fasta sequence dictionary*

Format of IonTorrent panel located at `/storage/research/dbmr_urology/Prostate_PDO/WG_IAD127899.20170720.designed.bed`
```
chr1	11174374	11174494	REGION_1_1.15817	.	GENE_ID=REGION_1;Pool=1
chr1	11188005	11188099	REGION_4_1.18508	.	GENE_ID=REGION_4;Pool=1
...
chrX	123229137	123229260	STAG2_167026_1.40045	.	GENE_ID=STAG2_167026;Pool=1
chrX	123229236	123229348	STAG2_167026_1.63590	.	GENE_ID=STAG2_167026;Pool=2
```
Format of IonTorrent panel intervals located at `/storage/research/dbmr_urology/Prostate_PDO/SPICE_data/WG_IAD127899.20170720.interval_list` that was generated with `picard.jar BedToIntervalList`
```
@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:249250621	M5:1b22b98cdeb4a9304cb5d48026a85128	UR:file:/storage/research/dbmr_urology/Prostate_PDO/hg19.fa
@SQ	SN:chr2	LN:243199373	M5:a0d9851da00400dec1098a9255ac712e	UR:file:/storage/research/dbmr_urology/Prostate_PDO/hg19.fa
@SQ	SN:chr3	LN:198022430	M5:641e4338fa8d52a5b781bd2a2c08d3c3	UR:file:/storage/research/dbmr_urology/Prostate_PDO/hg19.fa
@SQ	SN:chr4	LN:191154276	M5:23dccd106897542ad87d2765d28a19a1	UR:file:/storage/research/dbmr_urology/Prostate_PDO/hg19.fa
```
Format of `/storage/research/dbmr_urology/Prostate_PDO/SPICE_data/Homo_sapiens_assembly19.dbsnp.vcf` downloaded from `https://console.cloud.google.com/storage/browser/gcp-public-data--broad-references/hg19/v0`
```
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    
1       10144   rs144773400     TA      T       .       PASS    ASP;RSPOS=10145;SAO=0;SSR=0;VC=DIV;VP=050000000004000000000200;WGT=0;dbSNPBuildID=134
1       10228   rs143255646     TA      T       .       PASS    ASP;RSPOS=10229;SAO=0;SSR=0;VC=DIV;VP=050000000004000000000200;WGT=0;dbSNPBuildID=134
1       10234   rs145599635     C       T       .       PASS    ASP;RSPOS=10234;SAO=0;SSR=0;VC=SNV;VP=050000000004000000000100;WGT=0;dbSNPBuildID=134
```
Format of `/storage/research/dbmr_urology/Prostate_PDO/SPICE_data/Homo_sapiens_assembly19.dict` downloaded from `https://console.cloud.google.com/storage/browser/gcp-public-data--broad-references/hg19/v0`

```
@HD     VN:1.0  SO:unsorted
@SQ     SN:1    LN:249250621    UR:http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta   AS:GRCh37       M5:1
b22b98cdeb4a9304cb5d48026a85128     SP:Homo Sapiens
@SQ     SN:2    LN:243199373    UR:http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta   AS:GRCh37       M5:a
0d9851da00400dec1098a9255ac712e     SP:Homo Sapiens
@SQ     SN:3    LN:198022430    UR:http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta   AS:GRCh37       M5:f
dfd811849cc2fadebc929bb925902e5     SP:Homo Sapiens
```

**ISSUE**: *How to rename 1 to chr1 in Homo_sapiens_assembly19.dbsnp.vcf and Homo_sapiens_assembly19.dict OR rename chr1 to 1 in WG_IAD127899.20170720.interval_list?* OR maybe there is dbsnp with format chr1?
