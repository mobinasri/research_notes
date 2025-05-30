## Create assembly-based truth sets using HPRC-R2 assemblies

The samples HG001-007 are commonly being used for training models for variant calling and the models might overfit to these samples therefore testing on the same samples
cannot reveal such bias. To measure how much the models are overfit to these samples we should create some high-quality truth sets using samples other than HG001-007. 
I follow the steps mentioned below to create some new truth sets using HPRC-Release2 data:
- Select 5 HPRC-R2 samples with phased diploid assemblies polished with DeepPolisher (using PacBio-HiFi and ONT data).
  Those samples should be absent from HPRC-R1 since we want to use graph-v1.1 for testing pangenome-aware DeepVariant.
- Map assemblies to GRCh38 or CHM13-v2.0 and perform variant calling with dipcall (I prefer CHM13-v2.0 but I have to check if DRAGEN works with it)
- Copy output VCF and high-confidence BED file along with related illumina read data to the gs bucket `gs://pepper-deepvariant/mobinasri`


### Select five samples
```
cd /private/groups/patenlab/masri/internship/assembly_truth_sets
# Download list of files for HPRC-R2 assemblies
wget https://raw.githubusercontent.com/human-pangenomics/hprc_intermediate_assembly/refs/heads/main/data_tables/assemblies_pre_release_v0.6.1.index.csv

# Download list of HPRC-R1 samples
wget https://raw.githubusercontent.com/human-pangenomics/HPP_Year1_Data_Freeze_v1.0/refs/heads/main/sample_metadata/hprc_year1_sample_metadata.txt
```

```
# select 5 samples randomly
cat assemblies_pre_release_v0.6.1.index.csv | head -n1 > assemblies_pre_release_v0.6.1.index.selected_five.csv
cat hprc_year1_sample_metadata.txt | \
    cut -f1 | \
    grep -v -F -f - assemblies_pre_release_v0.6.1.index.csv | \
    shuf -n5 >> assemblies_pre_release_v0.6.1.index.selected_five.csv

# list of R2 samples
cat assemblies_pre_release_v0.6.1.index.selected_five.csv | awk -F',' '{print $1}'
sample_id
HG01975
NA20809
NA20346
HG03209
NA18879
```
### Download references
In the dipcall github page it was mentioned that:

> The above is applied to autosomes and female chrX. For a male sample, parent1 is assumed to be the father and parent2 the mother. Dipcall treats PARs the same way as autosomes. However, outside PARs, dipcall filters out chrX regions covered by father contigs and filters out chrY regions covered by mother contigs. To make proper calls on sex chromosomes, users should hard mask PARs on the reference chrY.

#### CHM13-v2.0
```
# chm13 reference
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz

# PAR regions
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0_PAR.bed

# remove chrY from PAR Bed
cat chm13v2.0_PAR.bed | grep chrX > chm13v2.0_PAR.only_X.bed
cat chm13v2.0_PAR.bed | grep chrY > chm13v2.0_PAR.only_Y.bed


# hard-mask PAR-chrY
gunzip chm13v2.0.fa.gz
bedtools maskfasta -fi chm13v2.0.fa -bed chm13v2.0_PAR.only_Y.bed -fo chm13v2.0.PAR_Y_masked.fa

# index
samtools faidx chm13v2.0.PAR_Y_masked.fa
```

#### GRCH38
```
# GRCh38 reference
```
### Run Dipcall
### Convert CRAM to FASTQ for short read data
### Copy to gs bucket

