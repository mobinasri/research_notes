# set working directory
WORKING_DIR="/private/groups/patenlab/masri/haplotype_sampling/graph_v2.0/remove_chm13"
cd ${WORKING_DIR}

## Save WDL path and name in environment variables
WDL_PATH="/private/groups/patenlab/masri/apps/vg_wdl/workflows/remove_sample_from_gbz.wdl"
WDL_FILENAME=$(basename ${WDL_PATH})
WDL_NAME=${WDL_FILENAME%%.wdl}


## Make a folder for saving files related to run e.g. input and output jsons
cd ${WORKING_DIR}
mkdir -p runs_toil_slurm
cd runs_toil_slurm

## Make a directory for saving input json files
mkdir -p ${WDL_NAME}_input_jsons
cd ${WDL_NAME}_input_jsons

python3 /private/groups/patenlab/masri/apps/hprc_intermediate_assembly/hpc/launch_from_table.py \
     --data_table ${WORKING_DIR}/data_table.csv \
     --field_mapping ${WORKING_DIR}/input_mapping.csv \
     --workflow_name ${WDL_NAME}


## Make sure you are in the working directory
cd ${WORKING_DIR}

## Set environment variables for sbatch
USERNAME="masri"
EMAIL="masri@ucsc.edu"
TIME_LIMIT="48:00:00"

## Partition should be modifed based on the available partitions on the server
PARTITION="long"


## Go to the execution directory
mkdir -p runs_toil_slurm/${WDL_NAME}_logs
cd runs_toil_slurm

## Run jobs arrays
sbatch      --job-name=${WDL_NAME}_${USERNAME} \
            --cpus-per-task=32 \
            --mem=128G \
            --mail-user=${EMAIL} \
            --mail-type=FAIL,END \
            --output=${WDL_NAME}_logs/${WDL_NAME}_%A_%a.log \
            --array=[1-2]%2  \
            --time=${TIME_LIMIT} \
            --partition=${PARTITION} \
            /private/groups/hprc/qc_hmm_flagger/hprc_intermediate_assembly/hpc/toil_sbatch_single_machine.sh \
            --wdl ${WDL_PATH} \
            --sample_csv  ${WORKING_DIR}/data_table.csv \
            --input_json_path ${WORKING_DIR}/runs_toil_slurm/${WDL_NAME}_input_jsons/\${SAMPLE_ID}_${WDL_NAME}.json
