
```
WORKING_DIR="/private/groups/patenlab/masri/t2t/HG002_v1.0.1/falsifier_runs/run_July_2024/tune_alpha_with_workflow_start_only"

cd ${WORKING_DIR}

## Get the script for creating input json files.
wget https://raw.githubusercontent.com/human-pangenomics/hprc_intermediate_assembly/1f61ff0043442d8350a282ef3533def588bee8dc/hpc/launch_from_table.py

## Save WDL path and name in environment variables
WDL_PATH="/private/groups/patenlab/masri/flagger/start_only_mode/flagger/wdls/workflows/tune_hyperparameter_alpha.wdl"
WDL_FILENAME=$(basename ${WDL_PATH})
WDL_NAME=${WDL_FILENAME%%.wdl}

## Make a folder for saving files related to run e.g. input and output jsons
mkdir -p runs_toil_slurm
cd runs_toil_slurm

## Make a directory for saving input json files
mkdir -p ${WDL_NAME}_input_jsons
cd ${WDL_NAME}_input_jsons

## Make input json files
## One json will be created per row
python3  ${WORKING_DIR}/launch_from_table.py \
            --data_table ${WORKING_DIR}/data_table.csv \
            --field_mapping ${WORKING_DIR}/input_mapping.csv \
            --workflow_name ${WDL_NAME}


## Make sure you are in the working directory. Check step 1 for setting ${WORKING_DIR} if it's not set already
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
mkdir -p runs_toil_slurm/${WDL_NAME}_logs
cd runs_toil_slurm

## Run jobs arrays
sbatch      --job-name=${WDL_NAME}_${USERNAME} \
            --cpus-per-task=64 \
            --mem=256G \
            --mail-user=${EMAIL} \
            --output=${WDL_NAME}_logs/${WDL_NAME}_%A_%a.log \
            --array=1-1%1  \
            --time=${TIME_LIMIT} \
            --partition=${PARTITION} \
            ${WORKING_DIR}/toil_sbatch_single_machine.sh \
            --wdl ${WDL_PATH} \
            --sample_csv  ${WORKING_DIR}/data_table.csv \
            --input_json_path ${WORKING_DIR}/runs_toil_slurm/${WDL_NAME}_input_jsons/\${SAMPLE_ID}_${WDL_NAME}.json
```
