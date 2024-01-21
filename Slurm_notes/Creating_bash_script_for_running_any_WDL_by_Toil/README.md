## Creating bash script for running any WDL by Toil
### Comment 1: 01/21/2024

I created a new bash script for running any WDL by Toil when we have a data table and we want run the related workflow on each row of the table.
Julian has already written a wdl for hifiasm [here](https://github.com/human-pangenomics/hprc_intermediate_assembly/blob/main/assembly/batch1/launch_hifiasm_array.sh) 
I just modified that script to make it work with any WDL and take the path to WDL as an argument. The file is temporarily uploaded in https://github.com/mobinasri/research_notes/blob/main/Slurm_notes/Creating_bash_script_for_running_any_WDL_by_Toil/files/launch_workflow_job_array.sh.

One can use this script using the commands below:
```
# Set some environment variables

USER_NAME="your username"
EMAIL="your email"

WDL_PATH="/path/to/workflow.wdl"
WDL_NAME=$(basename ${WDL_PATH%%.wdl})

INPUT_MAPPING_CSV="/path/to/input_mapping.csv"
INPUT_DATA_TABLE_CSV="/path/to/data_table.csv"

WORKING_DIR="/path/to/working/dir"

# create the working directory
mkdir -p ${WORKING_DIR}
cd ${WORKING_DIR}

mkdir ${WDL_NAME}_output_jsons
mkdir ${WDL_NAME}_input_jsons

# get the python script for making input json files
# git clone https://github.com/human-pangenomics/hprc_intermediate_assembly
# available in this path hprc_intermediate_assembly/hpc/launch_from_table.py

LAUNCH_FROM_TABLE_PY="path/to/hprc_intermediate_assembly/hpc/launch_from_table.py"

LAUNCH_WORKFLOW_JOB_ARRAY_BASH="/private/groups/patenlab/masri/hprc/hprc_intermediate_assembly/hpc/launch_workflow_job_array.sh"

```
```

cd ${WORKING_DIR}/${WDL_NAME}_input_jsons

python3  ${LAUNCH_FROM_TABLE_PY} \
    --data_table ${INPUT_DATA_TABLE_CSV} \
    --field_mapping ${INPUT_MAPPING_CSV} \
    --workflow_name ${WDL_NAME}

```

```
cd ${WORKING_DIR}
cd ${WDL_NAME}_output_jsons

TOIL_SLURM_ARGS="--time=3-0:00 --partition=long"
INPUT_JSON_DIR=${WORKING_DIR}/${WDL_NAME}_input_jsons

mkdir ${WDL_NAME}_logs

# modify --array and --partition based on what you need
sbatch \
     --job-name=${WDL_NAME}_${USERNAME} \
     --mail-user=${EMAIL} \
     --output=${WDL_NAME}_logs/${WDL_NAME}_%x_%j_%A_%a.log \
     --array=1-10%10 \
     --partition=long \
     ${LAUNCH_WORKFLOW_JOB_ARRAY_BASH} \
     ${INPUT_DATA_TABLE_CSV} \
     ${WDL_PATH} \
     "${TOIL_SLURM_ARGS}" \
     ${INPUT_JSON_DIR}
```

