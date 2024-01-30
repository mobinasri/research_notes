## Run Falsifier for HG002-T2T-v1.0.1

### Comment 1 : 01/30/2023

#### Map maternal haplotype against paternal and a create a PAF that can be used by Falsifier

```
USER_NAME="masri"
EMAIL="masri@ucsc.edu"

WDL_PATH="/private/groups/patenlab/masri/apps/flagger/wdls/tasks/alignment/asm2asm_aligner.wdl"
WDL_NAME=$(basename ${WDL_PATH%%.wdl})

INPUT_MAPPING_CSV="/private/groups/patenlab/masri/flagger/T2T/v1.0.1/mat_to_pat_alignment/input_mapping.csv"
INPUT_DATA_TABLE_CSV="/private/groups/patenlab/masri/flagger/T2T/v1.0.1/mat_to_pat_alignment/data_table.csv"

WORKING_DIR="/private/groups/patenlab/masri/flagger/T2T/v1.0.1/mat_to_pat_alignment"

# create the working directory
mkdir -p ${WORKING_DIR}
cd ${WORKING_DIR}

mkdir ${WDL_NAME}_output_jsons
mkdir ${WDL_NAME}_input_jsons


LAUNCH_FROM_TABLE_PY="/private/groups/patenlab/masri/apps/hprc_intermediate_assembly/hpc/launch_from_table.py"
LAUNCH_WORKFLOW_JOB_ARRAY_BASH="/private/groups/patenlab/masri/apps/hprc_intermediate_assembly/hpc/launch_workflow_job_array_single_machine.sh"

# create input json files
cd ${WORKING_DIR}/${WDL_NAME}_input_jsons
python3  ${LAUNCH_FROM_TABLE_PY} \
    --data_table ${INPUT_DATA_TABLE_CSV} \
    --field_mapping ${INPUT_MAPPING_CSV} \
    --workflow_name ${WDL_NAME}

cd ${WORKING_DIR}
cd ${WDL_NAME}_output_jsons
INPUT_JSON_DIR=${WORKING_DIR}/${WDL_NAME}_input_jsons
mkdir -p ${WDL_NAME}_logs

sbatch      --job-name=${WDL_NAME}_${USERNAME} \
            --cpus-per-task=32 \
            --mem=64G \
            --mail-user=${EMAIL} \
            --output=${WDL_NAME}_logs/${WDL_NAME}_%A_%a.log \
            --array=1-1%1  \
            --partition=medium  \
            --time=12:00:00 \
            ${LAUNCH_WORKFLOW_JOB_ARRAY_BASH} \
            --wdl ${WDL_PATH} \
            --sample_csv  ${INPUT_DATA_TABLE_CSV} \
            --input_json_dir ${INPUT_JSON_DIR}

```
