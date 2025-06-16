## Create assembly-based truth sets using HPRC-R2 assemblies

The HG001–HG007 samples are commonly used for training variant calling models. As a result, models may become overfit to these samples, and evaluating performance on the same data cannot effectively reveal such bias. To assess the extent of overfitting, it is essential to generate high-quality truth sets using samples outside of the HG001–HG007 set.
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
# select 5 samples randomly from trio assemblies (grep mat_hprc to consider only trio)
cat assemblies_pre_release_v0.6.1.index.csv | head -n1 > assemblies_pre_release_v0.6.1.index.selected_five.csv
cat hprc_year1_sample_metadata.txt | \
    cut -f1 | \
    grep -v -F -f - assemblies_pre_release_v0.6.1.index.csv | \
    grep mat_hprc | \
    shuf -n5 >> assemblies_pre_release_v0.6.1.index.selected_five.csv

# list of R2 samples
cat assemblies_pre_release_v0.6.1.index.selected_five.csv | awk -F',' '{print $1}'
sample_id
HG02984
HG02280
HG04184
HG01255
HG03831
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
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/GRCh38/assemblies/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz

# PAR regions
wget https://raw.githubusercontent.com/lh3/dipcall/refs/heads/master/data/hs38.PAR.bed

# remove chrY from PAR Bed
cat hs38.PAR.bed | grep chrX > hs38.PAR.only_X.bed
cat hs38.PAR.bed | grep chrY > hs38.PAR.only_Y.bed

# hard-mask PAR-chrY
gunzip -c GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz > hs38.fa
bedtools maskfasta -fi hs38.fa -bed hs38.PAR.only_Y.bed -fo hs38.PAR_Y_masked.fa

# index
samtools faidx hs38.PAR_Y_masked.fa
```
### Run Dipcall
I wrote a bash script for running Dipcall, which is available in `dipcall_files/dipcall.bash`. I prepared input json for the 5 samples selected above. Dipcall was run twice for each sample; once using GRCh38 and once using chm13-v2.0 as the reference. The json files are available in `dipcall_files/GRCh38_based` and `dipcall_files/CHM13_based`. I submitted dipcall jobs on the UCSC Pheonix cluster using the following commands:
```
# The input json files for the CHM13 reference are copied to this directory
cd /private/groups/patenlab/masri/internship/assembly_truth_sets/dipcall_results_CHM13
for i in $(ls | grep json);do sbatch --cpus-per-task=32 ../../../apps/bash_scripts/dipcall.bash $i  ;done

# The input json files for the GRCh38 reference are copied to this directory
cd /private/groups/patenlab/masri/internship/assembly_truth_sets/dipcall_results_GRCh38
for i in $(ls | grep json);do sbatch --cpus-per-task=32 ../../../apps/bash_scripts/dipcall.bash $i  ;done
```

To make the Dipcall call sets more conservative and easier to use for benchmarking I restricted them to autosomes. Therefore I created a new BED file per Dipcall run by removing the sex chromosomes from the Dipcall high-confidence BED file.
```
# for CHM13
cd /private/groups/patenlab/masri/internship/assembly_truth_sets/dipcall_results_CHM13
for i in $(find . | grep "dip.bed");do cat $i | grep -v -e"X" -e"Y" > ${i%%.bed}.no_XY.bed;done

# for GRCh38
cd /private/groups/patenlab/masri/internship/assembly_truth_sets/dipcall_results_GRCh38
for i in $(find . | grep "dip.bed");do cat $i | grep -v -e"X" -e"Y" > ${i%%.bed}.no_XY.bed;done
```

Upload Dipcall outputs:
```
# for GRCh38
cd /private/groups/patenlab/masri/internship/assembly_truth_sets/dipcall_results_GRCh38
for i in $(find . | grep -e "dip.vcf.gz" -e "dip.no_XY.bed" -e "dip.bed");do gsutil cp $i gs://pepper-deepvariant/mobinasri/pangenome_paper/hprc_r2_dipcall_call_sets/GRCh38_based/${i##./};done

# for CHM13
cd /private/groups/patenlab/masri/internship/assembly_truth_sets/dipcall_results_CHM13
for i in $(find . | grep -e "dip.vcf.gz" -e "dip.no_XY.bed" -e "dip.bed");do gsutil cp $i gs://pepper-deepvariant/mobinasri/pangenome_paper/hprc_r2_dipcall_call_sets/CHM13_based/${i##./};done
```

List of files uploaded to gcp
```
################
#### GRCh38 ####
################

gs://pepper-deepvariant/mobinasri/pangenome_paper/hprc_r2_dipcall_call_sets/GRCh38_based/HG01255/:
gs://pepper-deepvariant/mobinasri/pangenome_paper/hprc_r2_dipcall_call_sets/GRCh38_based/HG01255/HG01255_hprc_r2_v1.0.1_dipcall_GRCh38.dip.bed
gs://pepper-deepvariant/mobinasri/pangenome_paper/hprc_r2_dipcall_call_sets/GRCh38_based/HG01255/HG01255_hprc_r2_v1.0.1_dipcall_GRCh38.dip.no_XY.bed
gs://pepper-deepvariant/mobinasri/pangenome_paper/hprc_r2_dipcall_call_sets/GRCh38_based/HG01255/HG01255_hprc_r2_v1.0.1_dipcall_GRCh38.dip.vcf.gz

gs://pepper-deepvariant/mobinasri/pangenome_paper/hprc_r2_dipcall_call_sets/GRCh38_based/HG02280/:
gs://pepper-deepvariant/mobinasri/pangenome_paper/hprc_r2_dipcall_call_sets/GRCh38_based/HG02280/HG02280_hprc_r2_v1.0.1_dipcall_GRCh38.dip.bed
gs://pepper-deepvariant/mobinasri/pangenome_paper/hprc_r2_dipcall_call_sets/GRCh38_based/HG02280/HG02280_hprc_r2_v1.0.1_dipcall_GRCh38.dip.no_XY.bed
gs://pepper-deepvariant/mobinasri/pangenome_paper/hprc_r2_dipcall_call_sets/GRCh38_based/HG02280/HG02280_hprc_r2_v1.0.1_dipcall_GRCh38.dip.vcf.gz

gs://pepper-deepvariant/mobinasri/pangenome_paper/hprc_r2_dipcall_call_sets/GRCh38_based/HG02984/:
gs://pepper-deepvariant/mobinasri/pangenome_paper/hprc_r2_dipcall_call_sets/GRCh38_based/HG02984/HG02984_hprc_r2_v1.0.1_dipcall_GRCh38.dip.bed
gs://pepper-deepvariant/mobinasri/pangenome_paper/hprc_r2_dipcall_call_sets/GRCh38_based/HG02984/HG02984_hprc_r2_v1.0.1_dipcall_GRCh38.dip.no_XY.bed
gs://pepper-deepvariant/mobinasri/pangenome_paper/hprc_r2_dipcall_call_sets/GRCh38_based/HG02984/HG02984_hprc_r2_v1.0.1_dipcall_GRCh38.dip.vcf.gz

gs://pepper-deepvariant/mobinasri/pangenome_paper/hprc_r2_dipcall_call_sets/GRCh38_based/HG03831/:
gs://pepper-deepvariant/mobinasri/pangenome_paper/hprc_r2_dipcall_call_sets/GRCh38_based/HG03831/HG03831_hprc_r2_v1.0.1_dipcall_GRCh38.dip.bed
gs://pepper-deepvariant/mobinasri/pangenome_paper/hprc_r2_dipcall_call_sets/GRCh38_based/HG03831/HG03831_hprc_r2_v1.0.1_dipcall_GRCh38.dip.no_XY.bed
gs://pepper-deepvariant/mobinasri/pangenome_paper/hprc_r2_dipcall_call_sets/GRCh38_based/HG03831/HG03831_hprc_r2_v1.0.1_dipcall_GRCh38.dip.vcf.gz

gs://pepper-deepvariant/mobinasri/pangenome_paper/hprc_r2_dipcall_call_sets/GRCh38_based/HG04184/:
gs://pepper-deepvariant/mobinasri/pangenome_paper/hprc_r2_dipcall_call_sets/GRCh38_based/HG04184/HG04184_hprc_r2_v1.0.1_dipcall_GRCh38.dip.bed
gs://pepper-deepvariant/mobinasri/pangenome_paper/hprc_r2_dipcall_call_sets/GRCh38_based/HG04184/HG04184_hprc_r2_v1.0.1_dipcall_GRCh38.dip.no_XY.bed
gs://pepper-deepvariant/mobinasri/pangenome_paper/hprc_r2_dipcall_call_sets/GRCh38_based/HG04184/HG04184_hprc_r2_v1.0.1_dipcall_GRCh38.dip.vcf.gz

##################
#### CHM13-v2 ####
##################

gs://pepper-deepvariant/mobinasri/pangenome_paper/hprc_r2_dipcall_call_sets/CHM13_based/HG01255/:
gs://pepper-deepvariant/mobinasri/pangenome_paper/hprc_r2_dipcall_call_sets/CHM13_based/HG01255/HG01255_hprc_r2_v1.0.1_dipcall_CHM13.dip.bed
gs://pepper-deepvariant/mobinasri/pangenome_paper/hprc_r2_dipcall_call_sets/CHM13_based/HG01255/HG01255_hprc_r2_v1.0.1_dipcall_CHM13.dip.no_XY.bed
gs://pepper-deepvariant/mobinasri/pangenome_paper/hprc_r2_dipcall_call_sets/CHM13_based/HG01255/HG01255_hprc_r2_v1.0.1_dipcall_CHM13.dip.vcf.gz

gs://pepper-deepvariant/mobinasri/pangenome_paper/hprc_r2_dipcall_call_sets/CHM13_based/HG02280/:
gs://pepper-deepvariant/mobinasri/pangenome_paper/hprc_r2_dipcall_call_sets/CHM13_based/HG02280/HG02280_hprc_r2_v1.0.1_dipcall_CHM13.dip.bed
gs://pepper-deepvariant/mobinasri/pangenome_paper/hprc_r2_dipcall_call_sets/CHM13_based/HG02280/HG02280_hprc_r2_v1.0.1_dipcall_CHM13.dip.no_XY.bed
gs://pepper-deepvariant/mobinasri/pangenome_paper/hprc_r2_dipcall_call_sets/CHM13_based/HG02280/HG02280_hprc_r2_v1.0.1_dipcall_CHM13.dip.vcf.gz

gs://pepper-deepvariant/mobinasri/pangenome_paper/hprc_r2_dipcall_call_sets/CHM13_based/HG02984/:
gs://pepper-deepvariant/mobinasri/pangenome_paper/hprc_r2_dipcall_call_sets/CHM13_based/HG02984/HG02984_hprc_r2_v1.0.1_dipcall_CHM13.dip.bed
gs://pepper-deepvariant/mobinasri/pangenome_paper/hprc_r2_dipcall_call_sets/CHM13_based/HG02984/HG02984_hprc_r2_v1.0.1_dipcall_CHM13.dip.no_XY.bed
gs://pepper-deepvariant/mobinasri/pangenome_paper/hprc_r2_dipcall_call_sets/CHM13_based/HG02984/HG02984_hprc_r2_v1.0.1_dipcall_CHM13.dip.vcf.gz

gs://pepper-deepvariant/mobinasri/pangenome_paper/hprc_r2_dipcall_call_sets/CHM13_based/HG03831/:
gs://pepper-deepvariant/mobinasri/pangenome_paper/hprc_r2_dipcall_call_sets/CHM13_based/HG03831/HG03831_hprc_r2_v1.0.1_dipcall_CHM13.dip.bed
gs://pepper-deepvariant/mobinasri/pangenome_paper/hprc_r2_dipcall_call_sets/CHM13_based/HG03831/HG03831_hprc_r2_v1.0.1_dipcall_CHM13.dip.no_XY.bed
gs://pepper-deepvariant/mobinasri/pangenome_paper/hprc_r2_dipcall_call_sets/CHM13_based/HG03831/HG03831_hprc_r2_v1.0.1_dipcall_CHM13.dip.vcf.gz

gs://pepper-deepvariant/mobinasri/pangenome_paper/hprc_r2_dipcall_call_sets/CHM13_based/HG04184/:
gs://pepper-deepvariant/mobinasri/pangenome_paper/hprc_r2_dipcall_call_sets/CHM13_based/HG04184/HG04184_hprc_r2_v1.0.1_dipcall_CHM13.dip.bed
gs://pepper-deepvariant/mobinasri/pangenome_paper/hprc_r2_dipcall_call_sets/CHM13_based/HG04184/HG04184_hprc_r2_v1.0.1_dipcall_CHM13.dip.no_XY.bed
gs://pepper-deepvariant/mobinasri/pangenome_paper/hprc_r2_dipcall_call_sets/CHM13_based/HG04184/HG04184_hprc_r2_v1.0.1_dipcall_CHM13.dip.vcf.gz
```

### Convert CRAM to FASTQ for short read data
Since most tools work with fastq I convert CRAM to paired-end fastq files to make downstream analyses easier. The scripts and input json files for this conversion is available in `cram2fq/` directory.
```
# submit jobs on SLURM-based cluster
sbatch cram2fq.bash ${INPUT_SAMPLE_JSON}
```
The fastq files are uploaded to the gcp bucket:
```
gsutil ls -R gs://pepper-deepvariant/mobinasri/pangenome_paper/hprc_r2_short_reads | grep fastq.gz | sort
gs://pepper-deepvariant/mobinasri/pangenome_paper/hprc_r2_short_reads/HG01255/HG01255.illumina.R1.fastq.gz
gs://pepper-deepvariant/mobinasri/pangenome_paper/hprc_r2_short_reads/HG01255/HG01255.illumina.R2.fastq.gz
gs://pepper-deepvariant/mobinasri/pangenome_paper/hprc_r2_short_reads/HG02280/HG02280.illumina.R1.fastq.gz
gs://pepper-deepvariant/mobinasri/pangenome_paper/hprc_r2_short_reads/HG02280/HG02280.illumina.R2.fastq.gz
gs://pepper-deepvariant/mobinasri/pangenome_paper/hprc_r2_short_reads/HG02984/HG02984.illumina.R1.fastq.gz
gs://pepper-deepvariant/mobinasri/pangenome_paper/hprc_r2_short_reads/HG02984/HG02984.illumina.R2.fastq.gz
gs://pepper-deepvariant/mobinasri/pangenome_paper/hprc_r2_short_reads/HG03831/HG03831.illumina.R1.fastq.gz
gs://pepper-deepvariant/mobinasri/pangenome_paper/hprc_r2_short_reads/HG03831/HG03831.illumina.R2.fastq.gz
gs://pepper-deepvariant/mobinasri/pangenome_paper/hprc_r2_short_reads/HG04184/HG04184.illumina.R1.fastq.gz
gs://pepper-deepvariant/mobinasri/pangenome_paper/hprc_r2_short_reads/HG04184/HG04184.illumina.R2.fastq.gz
```

### Create haplotype-sampled graphs
To be able to run pangenome-aware DV with personalized pangenome I created haplotype-sampled graphs with 16 and 32 haplotypes. The scripts and input json files for running the haplotype sampling workflow is available in `haplotype_sampling/` directory. The WDL is available here: https://github.com/mobinasri/vg_wdl/blob/master/workflows/haplotype_sampling_customized.wdl

The gbz files are uploaded to the gcp bucket:
```
gsutil ls -R gs://pepper-deepvariant/mobinasri/pangenome_paper/hprc_r2_sampled_graphs | grep gbz | sort
gs://pepper-deepvariant/mobinasri/pangenome_paper/hprc_r2_sampled_graphs/HG01255.illumina/analysis/haplotype_sampling_customized_outputs/HG01255.illumina.hap_num_16.CHM13_removed.gbz
gs://pepper-deepvariant/mobinasri/pangenome_paper/hprc_r2_sampled_graphs/HG01255.illumina/analysis/haplotype_sampling_customized_outputs/HG01255.illumina.hap_num_32.CHM13_removed.gbz
gs://pepper-deepvariant/mobinasri/pangenome_paper/hprc_r2_sampled_graphs/HG02280.illumina/analysis/haplotype_sampling_customized_outputs/HG02280.illumina.hap_num_16.CHM13_removed.gbz
gs://pepper-deepvariant/mobinasri/pangenome_paper/hprc_r2_sampled_graphs/HG02280.illumina/analysis/haplotype_sampling_customized_outputs/HG02280.illumina.hap_num_32.CHM13_removed.gbz
gs://pepper-deepvariant/mobinasri/pangenome_paper/hprc_r2_sampled_graphs/HG02984.illumina/analysis/haplotype_sampling_customized_outputs/HG02984.illumina.hap_num_16.CHM13_removed.gbz
gs://pepper-deepvariant/mobinasri/pangenome_paper/hprc_r2_sampled_graphs/HG02984.illumina/analysis/haplotype_sampling_customized_outputs/HG02984.illumina.hap_num_32.CHM13_removed.gbz
gs://pepper-deepvariant/mobinasri/pangenome_paper/hprc_r2_sampled_graphs/HG03831.illumina/analysis/haplotype_sampling_customized_outputs/HG03831.illumina.hap_num_16.CHM13_removed.gbz
gs://pepper-deepvariant/mobinasri/pangenome_paper/hprc_r2_sampled_graphs/HG03831.illumina/analysis/haplotype_sampling_customized_outputs/HG03831.illumina.hap_num_32.CHM13_removed.gbz
gs://pepper-deepvariant/mobinasri/pangenome_paper/hprc_r2_sampled_graphs/HG04184.illumina/analysis/haplotype_sampling_customized_outputs/HG04184.illumina.hap_num_16.CHM13_removed.gbz
gs://pepper-deepvariant/mobinasri/pangenome_paper/hprc_r2_sampled_graphs/HG04184.illumina/analysis/haplotype_sampling_customized_outputs/HG04184.illumina.hap_num_32.CHM13_removed.gbz
```

### Map short reads to GRCh38 with vg-1.55
Short reads were mapped to graphv1.1 with vg-1.55 (I used this version to match the previous DeepVariant experiments). They were then surjected onto GRCh38 reference.
The scripts and input json files for mapping with vg giraffe are available in `map_short_reads/` directory.
```
# submit jobs on SLURM-based cluster for each sample
sbatch map_vg_bam.sh ${INPUT_SAMPLE_JSON}
```

The gcp urls:
```
gs://pepper-deepvariant/mobinasri/pangenome_paper/hprc_r2_short_read_mappings/graph_v1.1_grch38_vg_1.55/HG01255/HG01255.illumina.graph_v1.1_grch38.vg_1.55.bam
gs://pepper-deepvariant/mobinasri/pangenome_paper/hprc_r2_short_read_mappings/graph_v1.1_grch38_vg_1.55/HG01255/HG01255.illumina.graph_v1.1_grch38.vg_1.55.bam.bai
gs://pepper-deepvariant/mobinasri/pangenome_paper/hprc_r2_short_read_mappings/graph_v1.1_grch38_vg_1.55/HG02280/HG02280.illumina.graph_v1.1_grch38.vg_1.55.bam
gs://pepper-deepvariant/mobinasri/pangenome_paper/hprc_r2_short_read_mappings/graph_v1.1_grch38_vg_1.55/HG02280/HG02280.illumina.graph_v1.1_grch38.vg_1.55.bam.bai
gs://pepper-deepvariant/mobinasri/pangenome_paper/hprc_r2_short_read_mappings/graph_v1.1_grch38_vg_1.55/HG02984/HG02984.illumina.graph_v1.1_grch38.vg_1.55.bam
gs://pepper-deepvariant/mobinasri/pangenome_paper/hprc_r2_short_read_mappings/graph_v1.1_grch38_vg_1.55/HG02984/HG02984.illumina.graph_v1.1_grch38.vg_1.55.bam.bai
gs://pepper-deepvariant/mobinasri/pangenome_paper/hprc_r2_short_read_mappings/graph_v1.1_grch38_vg_1.55/HG03831/HG03831.illumina.graph_v1.1_grch38.vg_1.55.bam
gs://pepper-deepvariant/mobinasri/pangenome_paper/hprc_r2_short_read_mappings/graph_v1.1_grch38_vg_1.55/HG03831/HG03831.illumina.graph_v1.1_grch38.vg_1.55.bam.bai
gs://pepper-deepvariant/mobinasri/pangenome_paper/hprc_r2_short_read_mappings/graph_v1.1_grch38_vg_1.55/HG04184/HG04184.illumina.graph_v1.1_grch38.vg_1.55.bam
gs://pepper-deepvariant/mobinasri/pangenome_paper/hprc_r2_short_read_mappings/graph_v1.1_grch38_vg_1.55/HG04184/HG04184.illumina.graph_v1.1_grch38.vg_1.55.bam.bai
```



