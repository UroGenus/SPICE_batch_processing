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
wget https://github.com/broadinstitute/picard/releases/download/2.26.5/picard.jar
module load Java/11.0.2
java -jar picard.jar CreateSequenceDictionary R=/storage/research/dbmr_urology/Prostate_PDO/hg19.fa O=/storage/research/dbmr_urology/Prostate_PDO/hg19.dict
java -jar picard.jar BedToIntervalList  I=/storage/research/dbmr_urology/Prostate_PDO/WG_IAD127899.20170720.designed.bed O=/storage/research/dbmr_urology/Prostate_PDO/WG_IAD127899.20170720.interval_list SD=/storage/research/dbmr_urology/Prostate_PDO/hg19.fa
```

##### WES
```
???
```

#### snps_in_kit_vcf_file

One needs to installu GATK. Failed to do it using [UBELIX EasyBuild instructions](https://hpc-unibe-ch.github.io/software/EasyBuild.html), it was installed by UBELIX admins. After the installation, run 

```
download Homo_sapiens_assembly19.dbsnp.vcf, Homo_sapiens_assembly19.dict, Homo_sapiens_assembly19.fasta, Homo_sapiens_assembly19.fasta.fai from https://console.cloud.google.com/storage/browser/gcp-public-data--broad-references/hg19/v0
module load Workspace/home GATK
gatk IndexFeatureFile -I Homo_sapiens_assembly19.dbsnp.vcf
gatk SelectVariants -R /storage/research/dbmr_urology/Prostate_PDO/SPICE_data/Homo_sapiens_assembly19.fasta -V Homo_sapiens_assembly19.dbsnp.vcf -L /storage/research/dbmr_urology/Prostate_PDO/Panos_test_dir/WG_IAD127899.20170720.wout.chr.interval_list -select-type-to-include SNP -O WG_IAD127899.20170720.snp

```

**ERROR MESSAGE BECAUSE OF CHROMOSOME FORMAT MISMATCH**: 
*A USER ERROR has occurred: Input files reference and features have incompatible contigs: No overlapping contigs found.
  reference contigs = [chr1, chr2, chr3, chr4, chr5, chr6, chr7, chrX, chr8, chr9, chr10, chr11, chr12, chr13, chr14, chr15, chr16, chr17, chr18, chr20, chrY, chr19, chr22, chr21, chr6_ssto_hap7, chr6_mcf_hap5, chr6_cox_hap2, chr6_mann_hap4, chr6_apd_hap1, chr6_qbl_hap6, chr6_dbb_hap3, chr17_ctg5_hap1, chr4_ctg9_hap1, chr1_gl000192_random, chrUn_gl000225, chr4_gl000194_random, chr4_gl000193_random, chr9_gl000200_random, chrUn_gl000222, chrUn_gl000212, chr7_gl000195_random, chrUn_gl000223, chrUn_gl000224, chrUn_gl000219, chr17_gl000205_random, chrUn_gl000215, chrUn_gl000216, chrUn_gl000217, chr9_gl000199_random, chrUn_gl000211, chrUn_gl000213, chrUn_gl000220, chrUn_gl000218, chr19_gl000209_random, chrUn_gl000221, chrUn_gl000214, chrUn_gl000228, chrUn_gl000227, chr1_gl000191_random, chr19_gl000208_random, chr9_gl000198_random, chr17_gl000204_random, chrUn_gl000233, chrUn_gl000237, chrUn_gl000230, chrUn_gl000242, chrUn_gl000243, chrUn_gl000241, chrUn_gl000236, chrUn_gl000240, chr17_gl000206_random, chrUn_gl000232, chrUn_gl000234, chr11_gl000202_random, chrUn_gl000238, chrUn_gl000244, chrUn_gl000248, chr8_gl000196_random, chrUn_gl000249, chrUn_gl000246, chr17_gl000203_random, chr8_gl000197_random, chrUn_gl000245, chrUn_gl000247, chr9_gl000201_random, chrUn_gl000235, chrUn_gl000239, chr21_gl000210_random, chrUn_gl000231, chrUn_gl000229, chrM, chrUn_gl000226, chr18_gl000207_random]
  features contigs = [1, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 2, 20, 21, 22, 3, 4, 5, 6, 7, 8, 9, GL000191.1, GL000192.1, GL000193.1, GL000194.1, GL000195.1, GL000196.1, GL000197.1, GL000198.1, GL000199.1, GL000200.1, GL000201.1, GL000202.1, GL000203.1, GL000204.1, GL000205.1, GL000206.1, GL000207.1, GL000208.1, GL000209.1, GL000210.1, GL000211.1, GL000212.1, GL000213.1, GL000214.1, GL000215.1, GL000216.1, GL000217.1, GL000218.1, GL000219.1, GL000220.1, GL000221.1, GL000222.1, GL000223.1, GL000224.1, GL000225.1, GL000226.1, GL000227.1, GL000228.1, GL000229.1, GL000230.1, GL000231.1, GL000232.1, GL000233.1, GL000234.1, GL000235.1, GL000236.1, GL000237.1, GL000238.1, GL000239.1, GL000240.1, GL000241.1, GL000242.1, GL000243.1, GL000244.1, GL000245.1, GL000246.1, GL000247.1, GL000248.1, GL000249.1, MT, X, Y]*

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

**ISSUE**: *How to rename 1 to chr1 in Homo_sapiens_assembly19.dbsnp.vcf and Homo_sapiens_assembly19.dict OR rename chr1 to 1 in WG_IAD127899.20170720.interval_list? OR maybe there is dbsnp with format chr1 (e.g. there is something [here](https://genome.ucsc.edu/cgi-bin/hgTables?hgta_table=cytoBand&hgta_doSchema=describe%20table%20schema), but I'm not sure)?*

**ALSO TRIED**: Following [this instructions](https://www.biostars.org/p/410789/), but GCF_000001405.25_GRCh37.p13_assembly_report.txt does not seem to have all chromosome names that we have in our panel.
