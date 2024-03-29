## Preparing workflows for haplotype sampling and making GBZ files with different number of haplotypes
### Comment 1 : 01/29/2023

#### Modifying the haplotype sampling workflow
I created a fork of https://github.com/vgteam/vg_wdl in here https://github.com/mobinasri/vg_wdl to be able to make the changes I want in the workflows.
I made a customized version of Parsa's workflow for haplotype sampling. The original version is here: https://github.com/vgteam/vg_wdl/blob/master/workflows/haplotype_sampling.wdl and 
my version is here: https://github.com/mobinasri/vg_wdl/blob/master/workflows/haplotype_sampling_customized.wdl

The customized version can take an array of read files of any type. I imported `extract_reads.wdl` for extracting each read file, which can be either BAM,CRAM,FASTQ or FASTA.GZ
and receiving a fastq file. This haplotype sampling workflow will take a list of haplotype numbers and will create an array of gbz files; each one is a graph created with 
a different haplotype number. My main motivation for writing this WDL was that we are going to train multiple models with different number of haplotypes and then investigate which one acheives
a higher accuracy. For this aim we need to create multiple gbz files for different haplotype numbers. With this WDL we can run the workflow only once and create all the graphs 
that we need. 

#### Preparing csv files for tesing the workflow

I'm going to use one of the bam files Parsa generated by mapping short reads to GRCh38 with Giraffe. I will use this one:
```
gs://pepper-deepvariant/seeskand/dv_training/HG003_35x_merged.positionsorted.bam
gs://pepper-deepvariant/seeskand/dv_training/HG003_35x_merged.positionsorted.bam.bai
```

I made these two csv files for my run on Phoenix:
[data_table.csv](https://github.com/mobinasri/research_notes/blob/main/DeepVariant_Pangenome/Prepare_workflows_and_tables_for_haplotype_sampling/files/data_table.csv)
and
[input_mapping.csv](https://github.com/mobinasri/research_notes/blob/main/DeepVariant_Pangenome/Prepare_workflows_and_tables_for_haplotype_sampling/files/input_mapping.csv)

I ran the commands below
```
# set env variables
USER_NAME="masri"
EMAIL="masri@ucsc.edu"
WDL_PATH="/private/groups/patenlab/masri/apps/vg_wdl/workflows/haplotype_sampling_customized.wdl"
WDL_NAME=$(basename ${WDL_PATH%%.wdl})
INPUT_MAPPING_CSV="/private/groups/patenlab/masri/haplotype_sampling/HG003/input_mapping.csv"
INPUT_DATA_TABLE_CSV="/private/groups/patenlab/masri/haplotype_sampling/HG003/data_table.csv"
WORKING_DIR="/private/groups/patenlab/masri/haplotype_sampling/HG003"
```
```
# get the scripts for running my job
cd /private/groups/patenlab/masri/apps
git clone https://github.com/human-pangenomics/hprc_intermediate_assembly
LAUNCH_FROM_TABLE_PY="/private/groups/patenlab/masri/apps/hprc_intermediate_assembly/hpc/launch_from_table.py"
LAUNCH_WORKFLOW_JOB_ARRAY_BASH="/private/groups/patenlab/masri/apps/hprc_intermediate_assembly/hpc/launch_workflow_job_array_single_machine.sh"
```
```
# create directories for putting input and output json files
mkdir -p ${WORKING_DIR}
cd ${WORKING_DIR}

mkdir -p ${WDL_NAME}_output_jsons
mkdir -p ${WDL_NAME}_input_jsons

# create input json files
cd ${WORKING_DIR}/${WDL_NAME}_input_jsons
python3  ${LAUNCH_FROM_TABLE_PY} \
    --data_table ${INPUT_DATA_TABLE_CSV} \
    --field_mapping ${INPUT_MAPPING_CSV} \
    --workflow_name ${WDL_NAME}

```
```
# Run the job
cd ${WORKING_DIR}
cd ${WDL_NAME}_output_jsons
INPUT_JSON_DIR=${WORKING_DIR}/${WDL_NAME}_input_jsons
mkdir -p ${WDL_NAME}_logs

sbatch      --job-name=${WDL_NAME}_${USERNAME} \
            --cpus-per-task=16 \
            --mem=180G \
            --mail-user=${EMAIL} \
            --output=${WDL_NAME}_logs/${WDL_NAME}_%A_%a.log \
            --array=1-1%1  \
            --partition=long  \
            --time=3-00:00:00 \
            ${LAUNCH_WORKFLOW_JOB_ARRAY_BASH} \
            --wdl ${WDL_PATH} \
            --sample_csv  ${INPUT_DATA_TABLE_CSV} \
            --input_json_dir ${INPUT_JSON_DIR}
```
