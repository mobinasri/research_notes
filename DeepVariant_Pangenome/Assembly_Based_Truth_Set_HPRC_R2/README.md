## Create assembly-based truth sets using HPRC-R2 assemblies

The samples HG001-007 are commonly being used for training models for variant calling and the models might overfit to these samples therefore testing on the same samples
cannot reveal such bias. To measure how much the models are overfit to these samples we should create some high-quality truth sets using samples other than HG001-007. 
I follow the steps mentioned below to create some new truth sets using HPRC-Release2 data:
- Select 5 HPRC-R2 samples with phased diploid assemblies polished with DeepPolisher (using PacBio-HiFi and ONT data).
  Those samples should be absent from HPRC-R1 since we want to use graph-v1.1 for testing pangenome-aware DeepVariant.
- Map assemblies to GRCh38 or CHM13-v2.0 and perform variant calling with dipcall (I prefer CHM13-v2.0 but I have to check if DRAGEN works with it)
- Copy output VCF and high-confidence BED file along with related illumina read data in the gs bucket `gs://pepper-deepvariant/mobinasri`


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
### Run Dipcall
### Convert CRAM to FASTQ for short read data
### Copy to gs bucket

