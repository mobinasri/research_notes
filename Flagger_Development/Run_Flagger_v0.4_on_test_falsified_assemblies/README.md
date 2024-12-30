# Run Flagger v0.4 on test falsified assemblies
Here I run non-HMM Flagger (v0.4) on test falsified assemblies (two misassembly rates; 0.5% and 2%) for 3 different platforms; HiFi_DC_1.2, ONT_R9 and ONT_R10. I'm running Flagger-v0.4 on each bam file with two downsampling rates; 1.0 and 0.5. All reads were mapped with minimap2 with appropriate presets (`lr:hqae` for HiFi and ONT-R10, `map-ont` for ONT-R9). 

## Create json files
```
cd /private/groups/patenlab/masri/t2t/HG002_v1.1/falsifier_runs_oct_2024/test/non_hmm_flagger
WORKING_DIR=${PWD}

## Make sure you are in the working directory
cd ${WORKING_DIR}

## Get the script for creating input json files.
wget https://raw.githubusercontent.com/human-pangenomics/hprc_intermediate_assembly/1f61ff0043442d8350a282ef3533def588bee8dc/hpc/launch_from_table.py

WDL_PATH=/private/groups/patenlab/masri/apps/flagger_v0.4.0/flagger/wdls/workflows/flagger_end_to_end.wdl
WDL_FILENAME=$(basename ${WDL_PATH})
WDL_NAME=${WDL_FILENAME%%.wdl}

## Make a directory for saving input json files
mkdir -p ${WDL_NAME}_input_jsons
cd ${WDL_NAME}_input_jsons

## Make input json files
## One json will be created per row
python3  ${WORKING_DIR}/launch_from_table.py \
            --data_table ${WORKING_DIR}/data_table.csv \
            --field_mapping ${WORKING_DIR}/input_mapping.csv \
            --workflow_name ${WDL_NAME}
```

## Run jobs (12 jobs in total)
```
## Make sure you are in the working directory
cd ${WORKING_DIR}

## Get the bash script for running WDLs on Slurm using Toil
wget https://raw.githubusercontent.com/human-pangenomics/hprc_intermediate_assembly/b81bbb9540eaf5632a53faba43be71a0974f14f6/hpc/toil_sbatch_single_machine.sh

## Set environment variables for sbatch
USERNAME="masri"
EMAIL="masri@ucsc.edu"
TIME_LIMIT="70:00:00"

## Partition should be modifed based on the available partitions on the server
PARTITION="long"

## Go to the execution directory
mkdir -p ${WDL_NAME}_logs

## Run jobs arrays
sbatch      --job-name=${WDL_NAME}_${USERNAME} \
            --cpus-per-task=32 \
            --mem=128G \
            --mail-user=${EMAIL} \
            --output=${WDL_NAME}_logs/${WDL_NAME}_%A_%a.log \
            --array=1-12%12  \
            --time=${TIME_LIMIT} \
            --partition=${PARTITION} \
            ${WORKING_DIR}/toil_sbatch_single_machine.sh \
            --wdl ${WDL_PATH} \
            --sample_csv  ${WORKING_DIR}/data_table.csv \
            --input_json_path ${WORKING_DIR}/${WDL_NAME}_input_jsons/\${SAMPLE_ID}_${WDL_NAME}.json
```
