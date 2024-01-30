## Run Falsifier for HG002-T2T-v1.0.1

### Comment 1 : 01/30/2023

#### Map maternal haplotype against paternal and a create a PAF that can be used by Falsifier
I made two csv files for running this alignment ([data_table.csv](https://github.com/mobinasri/research_notes/blob/main/Flagger_HG002_T2T/Run_Falsifier_for_v1.0.1/files/data_table.csv) and [input_mapping.csv](https://github.com/mobinasri/research_notes/blob/main/Flagger_HG002_T2T/Run_Falsifier_for_v1.0.1/files/input_mapping.csv))
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

### Comment 2 : 01/30/2023

After the alignment job is done I made a PAF file from BAM file:

```
# ssh into a node after being allowed by Slurm
salloc  -c 8 --mem 32G bash -c 'ssh -Y $(scontrol show hostnames | head -n 1)'

# get the docker with paftools
docker run -v$PWD:$PWD -it -u$(id -u):$(id -g) mobinasri/long_read_aligner:v0.3.0

# convert bam to paf
k8 /home/apps/minimap2-2.24/misc/paftools.js sam2paf <(samtools view -h /private/groups/patenlab/masri/flagger/T2T/v1.0.1/mat_to_pat_alignment/asm2asm_aligner_output_jsons/HG002_v1.0.1/asm2asm_aligner_outputs/hg002v1.0.1.mat.mat_to_pat.sorted.bam ) >  /private/groups/patenlab/masri/flagger/T2T/v1.0.1/mat_to_pat_alignment/asm2asm_aligner_output_jsons/HG002_v1.0.1/asm2asm_aligner_outputs/hg002v1.0.1.mat.mat_to_pat.sorted.paf
```
